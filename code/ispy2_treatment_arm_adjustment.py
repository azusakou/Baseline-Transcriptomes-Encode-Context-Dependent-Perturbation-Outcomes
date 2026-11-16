"""Matched I-SPY2 treatment-arm adjustment, without modifying existing results.

The existing 987 patient identifiers, outcome labels, order and held-out folds
are locked to clinical_predictions.tsv. Reproduce the corrected clinical()
pipeline in revision_analysis.py, including its training-fold ranked five
genes and its 8-dimensional, independently refitted PFE (not OOF probability
stacking). The only new predictor is treatment arm. No hyperparameter tuning.

Usage (from SpringerV1):
  python code/ispy2_treatment_arm_adjustment.py --help

Raw expression stays in the source archive. All new outputs use the prefix
ispy2_arm_adjustment in the separate output directory specified below.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy
from scipy import sparse
import sklearn
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (average_precision_score, balanced_accuracy_score,
                             brier_score_loss, roc_auc_score)
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import OneHotEncoder, StandardScaler
import torch
from torch import nn

HERE = Path(__file__).resolve().parent
SPRINGER = HERE.parent
SEED = 20261116
GENES = ['BAX', 'BCL2', 'MKI67', 'PCNA', 'ATF3']
PREFIX = 'ispy2_arm_adjustment'
# GEO GSE194040 describes one AMG-386 arm in the HER2-negative subset.
# The single hyphenated spelling is therefore treated as a typography alias.
# Retain a literal-label sensitivity analysis rather than silently collapsing it.
ARM_ALIASES = {'Paclitaxel + AMG-386': 'Paclitaxel + AMG 386'}


def strings(a):
    return np.asarray([x.decode() if isinstance(x, bytes) else x for x in a])


def hread(v):
    if isinstance(v, h5py.Dataset):
        a = v[:]
        return strings(a) if a.dtype.kind in 'OSU' else a
    kind = v.attrs.get('encoding-type', '')
    if kind == 'categorical':
        return pd.Categorical.from_codes(v['codes'][:], strings(v['categories'][:]))
    if 'values' in v:
        a = v['values'][:].astype(float)
        a[v['mask'][:]] = np.nan
        return a
    raise ValueError((v.name, kind))


def obs(path):
    with h5py.File(path) as h:
        g = h['obs']
        key = g.attrs.get('_index', '_index')
        return pd.DataFrame({k: hread(g[k]) for k in g if k != key},
                            index=pd.Index(hread(g[key]), name='cell_id'))


def gene_matrix(path, cells):
    with h5py.File(path) as h:
        g = h['obs']
        ids = pd.Index(hread(g[g.attrs.get('_index', '_index')]))
        idx = ids.get_indexer(cells)
        var = h['var']
        genes = pd.Index(hread(var[var.attrs.get('_index', '_index')]))
        gi = genes.get_indexer(GENES)
        assert (idx >= 0).all() and (gi >= 0).all()
        x = h['X']
        if isinstance(x, h5py.Dataset):
            return x[:][idx][:, gi]
        full = sparse.csr_matrix((x['data'][:], x['indices'][:], x['indptr'][:]),
                                 shape=tuple(x.attrs['shape']))
        # Preserve original matrix() float32 conversion for sparse inputs.
        return full[idx].astype(np.float32)[:, gi].toarray()


def rank_fit(a, b):
    """Exact original training-patient empirical CDF with midranks."""
    aa = np.empty_like(a, dtype=float)
    bb = np.empty_like(b, dtype=float)
    for j in range(a.shape[1]):
        ref = np.sort(a[np.isfinite(a[:, j]), j])
        n = len(ref)
        for v, dest in [(a[:, j], aa[:, j]), (b[:, j], bb[:, j])]:
            dest[:] = (np.searchsorted(ref, v, 'left') +
                       np.searchsorted(ref, v, 'right') + 1) / (2 * max(n, 1))
            dest[~np.isfinite(v)] = .5
    return aa, bb


class Encoder(nn.Module):
    def __init__(self, n):
        super().__init__()
        self.network = nn.Sequential(nn.Linear(n, 64), nn.ReLU(), nn.Dropout(.15), nn.Linear(64, 8))

    def forward(self, x):
        return nn.functional.normalize(self.network(x), dim=1)


def fitenc(x, y, seed):
    """Exact architecture, loss, random seeds and optimization of clinical()."""
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    net = Encoder(x.shape[1])
    opt = torch.optim.AdamW(net.parameters(), lr=.002, weight_decay=.001)
    for _ in range(80):
        order = rng.permutation(len(x))
        for start in range(0, len(x), 128):
            ii = order[start:start + 128]
            batch = torch.tensor(x[ii], dtype=torch.float32)
            lab = torch.tensor(np.tile(y[ii], 2))
            v = torch.cat([batch + .05 * torch.randn_like(batch), batch + .05 * torch.randn_like(batch)])
            z = net(v)
            sim = z @ z.T / .15
            sim = sim - sim.max(1, keepdim=True).values.detach()
            eye = torch.eye(len(z), dtype=torch.bool)
            positive = lab[:, None].eq(lab[None, :]) & ~eye
            lp = sim - torch.exp(sim).masked_fill(eye, 0).sum(1).clamp_min(1e-12).log()[:, None]
            loss = -(lp * positive).sum(1).div(positive.sum(1).clamp_min(1)).mean()
            opt.zero_grad()
            loss.backward()
            opt.step()
    net.eval()
    return net


def sha256(path):
    dig = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            dig.update(block)
    return dig.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project-root', type=Path, default=SPRINGER.parent)
    ap.add_argument('--data-root', type=Path, default=Path('source_data'))
    ap.add_argument('--reference-dir', type=Path, default=HERE.parent / 'results')
    ap.add_argument('--output-dir', type=Path, default=HERE.parent / 'results')
    ap.add_argument('--bootstrap-replicates', type=int, default=5000)
    args = ap.parse_args()
    assert args.bootstrap_replicates >= 1000
    args.output_dir.mkdir(parents=True, exist_ok=True)
    torch.set_num_threads(1)

    def save(rows, suffix):
        pd.DataFrame(rows).to_csv(args.output_dir / f'{PREFIX}_{suffix}.tsv', sep='\t', index=False)

    reference_path = args.reference_dir / 'clinical_predictions.tsv'
    fullref = pd.read_csv(reference_path, sep='\t', dtype={'patient': str, 'cell_id': str})
    ref = fullref.query("cohort == 'ispy2' and model == 'PFE'").copy().reset_index(drop=True)
    assert len(ref) == 987 and ref.patient.is_unique and ref.cell_id.is_unique
    y = ref.outcome.to_numpy(int)
    foldid = ref.fold.to_numpy(int)
    groups = ref.patient.to_numpy()
    assert set(np.unique(foldid)) == set(range(1, 6))

    expression_path = args.data_root / 'staged/gse194040/ispy2_pretreatment_expression.h5ad'
    embedding_path = args.data_root / 'results/patient_fate_embedding/ispy2/patient_fate_embedding.parquet'
    axes_path = args.data_root / 'processed/axes/ispy2_axes.parquet'
    weights_path = args.project_root / 'results/intrinsic_fate_axis_mapper_stability/intrinsic_fate_axis_weights.tsv'
    source_meta = obs(expression_path)
    meta = source_meta.loc[ref.cell_id].copy()
    old = pd.read_parquet(embedding_path).set_index('cell_id')
    # The original source order is important for stochastic encoder mini-batches.
    assert old.index.equals(pd.Index(ref.cell_id, name='cell_id'))
    assert np.array_equal(old.outcome.to_numpy(int), y)
    assert np.array_equal(old.patient.astype(str).to_numpy(), groups)
    assert np.array_equal(meta.patient_id.astype(str).to_numpy(), groups)
    assert np.array_equal(pd.to_numeric(meta.pcr).to_numpy(int), y)

    checks = []
    for col in ['arm', 'hr', 'her2']:
        saved = pd.read_csv(args.reference_dir / f'ispy2_{col}_metadata.tsv', sep='\t', dtype={'patient': str})
        assert saved.patient.is_unique
        saved = saved.set_index('patient').loc[groups, col]
        raw = meta[col]
        if col == 'arm':
            assert np.array_equal(saved.astype(str).to_numpy(), raw.astype(str).to_numpy())
        else:
            assert np.array_equal(pd.to_numeric(saved).to_numpy(float), pd.to_numeric(raw).to_numpy(float))
        checks.append({'check': f'original_{col}_metadata_matches_h5ad', 'passed': True})

    # Audit all source records for conflicting clinical labels within patients.
    conflict_rows = []
    for patient, frame in source_meta.groupby('patient_id', observed=True):
        for col in ['arm', 'hr', 'her2', 'pcr']:
            vals = frame[col].dropna().astype(str).unique()
            if len(vals) > 1:
                conflict_rows.append({'patient': patient, 'variable': col, 'values': '|'.join(vals)})
    save(pd.DataFrame(conflict_rows, columns=['patient', 'variable', 'values']), 'source_patient_conflicts')
    assert not conflict_rows, 'Conflicting source patient metadata requires explicit resolution'

    clinical = pd.DataFrame(index=meta.index)
    for col in ['hr', 'her2']:
        clinical[col] = meta[col].astype(str).replace({'nan': 'unknown', '<NA>': 'unknown'})
    arm = meta.arm.astype(str).str.strip().replace({'nan': 'unknown', '<NA>': 'unknown', '': 'unknown'})
    clinical['arm_literal'] = arm
    clinical['arm'] = arm.replace(ARM_ALIASES)
    audit = ref[['patient', 'cell_id', 'outcome', 'fold']].copy()
    for col in clinical:
        audit[col] = clinical[col].to_numpy()
    audit['arm_alias_changed'] = audit.arm.ne(audit.arm_literal)
    audit['geo_accession'] = meta.geo_accession.to_numpy()
    save(audit, 'patients')

    summary_rows = []
    for arm_column in ['arm', 'arm_literal']:
        for name, frame in audit.groupby(arm_column, dropna=False, sort=True):
            summary_rows.append({'arm_encoding': arm_column, 'arm': name,
                                 'n': len(frame), 'pcr_n': int(frame.outcome.sum()),
                                 'pcr_fraction': frame.outcome.mean(),
                                 'hr_positive_n': int(frame.hr.eq('1.0').sum()),
                                 'her2_positive_n': int(frame.her2.eq('1.0').sum()),
                                 'folds_present': frame.fold.nunique()})
    save(summary_rows, 'treatment_composition')
    save(audit.groupby(['arm', 'hr', 'her2', 'outcome'], observed=True).size().reset_index(name='n'),
         'treatment_receptor_outcome_composition')
    save(audit.groupby(['fold', 'arm', 'outcome'], observed=True).size().reset_index(name='n'),
         'fold_treatment_composition')

    weights = pd.read_csv(weights_path, sep='\t')
    selected = weights.assign(ab=weights.coefficient.abs()).sort_values('ab', ascending=False).groupby('target').head(20).feature.unique().tolist()
    allfeat = weights.feature.unique().tolist()
    aa = pd.read_parquet(axes_path).set_index('cell_id').loc[ref.cell_id]
    common = [f for f in allfeat if f in aa]
    selected_here = [f for f in selected if f in common]
    raw = aa[common].to_numpy(float)
    gx = gene_matrix(expression_path, ref.cell_id)
    save([{'gene': gene, 'missing_n': int((~np.isfinite(gx[:, j])).sum()),
           'n': len(gx), 'handling': 'rank_fit maps nonfinite values to 0.5, as in original clinical()'}
          for j, gene in enumerate(GENES)], 'five_gene_missingness')
    recomputed = np.zeros(len(y), int)
    for fold, (_, te) in enumerate(StratifiedGroupKFold(5, shuffle=True, random_state=20260713).split(raw, y, groups), 1):
        recomputed[te] = fold
    assert np.array_equal(recomputed, foldid), 'Saved folds differ from original splitter'
    checks.append({'check': 'stored_fold_ids_equal_original_patient_splitter', 'passed': True})

    model_names = ['arm_only', 'hr_her2_only', 'hr_her2_plus_arm', 'five_gene',
                   'hr_her2_plus_five_gene', 'hr_her2_plus_arm_plus_five_gene',
                   'PFE', 'hr_her2_plus_PFE', 'hr_her2_plus_arm_plus_PFE',
                   'arm_only_literal', 'hr_her2_plus_arm_literal',
                   'hr_her2_plus_arm_plus_five_gene_literal', 'hr_her2_plus_arm_plus_PFE_literal']
    pp = {name: np.full(len(y), np.nan) for name in model_names}
    foldmetrics, category_diagnostics = [], []
    embedding_rows = []
    for fold in range(1, 6):
        tr = np.flatnonzero(foldid != fold)
        te = np.flatnonzero(foldid == fold)
        assert not set(groups[tr]) & set(groups[te])
        trainrank, testrank = rank_fit(raw[tr], raw[te])
        features = []
        for rr in [trainrank, testrank]:
            potential = []
            for _, ww in weights.groupby('target', sort=True):
                ww = ww[ww.feature.isin(common)]
                ii = [common.index(v) for v in ww.feature]
                potential.append(((rr[:, ii] - ww.rank_reference_mean.to_numpy()) /
                                  ww.rank_reference_scale.to_numpy()) @ ww.coefficient.to_numpy() + ww.intercept.iloc[0])
            features.append(np.column_stack([rr[:, [common.index(v) for v in selected_here]], *potential]))
        tx, vx = rank_fit(*features)
        sc = StandardScaler().fit(tx)
        tx, vx = sc.transform(tx), sc.transform(vx)
        encoder = fitenc(tx, y[tr], 20260713 + fold)
        with torch.no_grad():
            tz = encoder(torch.tensor(tx, dtype=torch.float32)).numpy()
            vz = encoder(torch.tensor(vx, dtype=torch.float32)).numpy()
        for j, pos in enumerate(te):
            embedding_rows.append({'patient': groups[pos], 'fold': fold,
                                   **{f'PFE_{k + 1}': vz[j, k] for k in range(8)}})

        def onehot(columns):
            enc = OneHotEncoder(handle_unknown='ignore', sparse_output=False)
            tc = enc.fit_transform(clinical.iloc[tr][columns])
            vc = enc.transform(clinical.iloc[te][columns])
            for col, levels in zip(columns, enc.categories_):
                unseen = sorted(set(clinical.iloc[te][col]) - set(levels))
                category_diagnostics.append({'fold': fold, 'design': '+'.join(columns), 'column': col,
                                             'training_categories': '|'.join(map(str, levels)),
                                             'unseen_test_categories': '|'.join(unseen),
                                             'unseen_test_n': int(clinical.iloc[te][col].isin(unseen).sum())})
            return tc, vc

        c, cv = onehot(['hr', 'her2'])
        a, av = onehot(['arm'])
        ca, cav = onehot(['hr', 'her2', 'arm'])
        al, alv = onehot(['arm_literal'])
        cal, calv = onehot(['hr', 'her2', 'arm_literal'])
        gt, gv = rank_fit(gx[tr], gx[te])
        sg = StandardScaler().fit(gt)
        gt, gv = sg.transform(gt), sg.transform(gv)
        designs = {
            'arm_only': (a, av),
            'hr_her2_only': (c, cv),
            'hr_her2_plus_arm': (ca, cav),
            'five_gene': (gt, gv),
            'hr_her2_plus_five_gene': (np.column_stack([c, gt]), np.column_stack([cv, gv])),
            'hr_her2_plus_arm_plus_five_gene': (np.column_stack([ca, gt]), np.column_stack([cav, gv])),
            'PFE': (tz, vz),
            'hr_her2_plus_PFE': (np.column_stack([c, tz]), np.column_stack([cv, vz])),
            'hr_her2_plus_arm_plus_PFE': (np.column_stack([ca, tz]), np.column_stack([cav, vz])),
            'arm_only_literal': (al, alv),
            'hr_her2_plus_arm_literal': (cal, calv),
            'hr_her2_plus_arm_plus_five_gene_literal': (np.column_stack([cal, gt]), np.column_stack([calv, gv])),
            'hr_her2_plus_arm_plus_PFE_literal': (np.column_stack([cal, tz]), np.column_stack([calv, vz])),
        }
        for name, (xt, xv) in designs.items():
            pred = LogisticRegression(C=.2, max_iter=2000).fit(xt, y[tr]).predict_proba(xv)[:, 1]
            pp[name][te] = pred
            foldmetrics.append({'model': name, 'fold': fold, 'n_train': len(tr), 'n_test': len(te),
                                'pcr_n': int(y[te].sum()), 'auroc': roc_auc_score(y[te], pred),
                                'average_precision': average_precision_score(y[te], pred),
                                'brier': brier_score_loss(y[te], pred),
                                'balanced_accuracy': balanced_accuracy_score(y[te], pred >= .5)})
        print(f'Fitted fold {fold}: {len(tr)} training, {len(te)} held-out patients', flush=True)
    assert all(np.isfinite(v).all() for v in pp.values())
    save(foldmetrics, 'fold_metrics')
    save(category_diagnostics, 'category_diagnostics')
    save(embedding_rows, 'heldout_embedding')
    prediction_rows = []
    for model, pred in pp.items():
        frame = audit.copy()
        frame['model'] = model
        frame['predicted'] = pred
        prediction_rows.append(frame)
    save(pd.concat(prediction_rows, ignore_index=True), 'predictions')

    fidelity = []
    for original, new in [('PFE', 'PFE'), ('clinical', 'hr_her2_only'),
                          ('clinical_plus_PFE', 'hr_her2_plus_PFE'), ('five_gene', 'five_gene'),
                          ('clinical_plus_five_gene', 'hr_her2_plus_five_gene')]:
        prior = fullref.query("cohort == 'ispy2' and model == @original").set_index('patient').loc[groups]
        discrepancy = float(np.max(np.abs(prior.predicted.to_numpy() - pp[new])))
        fidelity.append({'original_model': original, 'rerun_model': new,
                         'max_abs_probability_difference': discrepancy,
                         'auroc_difference': roc_auc_score(y, pp[new]) - roc_auc_score(y, prior.predicted),
                         'average_precision_difference': average_precision_score(y, pp[new]) - average_precision_score(y, prior.predicted)})
        assert discrepancy < 1e-10, f'Original {original} predictions not reproduced ({discrepancy})'
    save(fidelity, 'original_prediction_fidelity')
    print('All five unchanged models reproduce existing predictions within 1e-10.', flush=True)

    # One row per unique patient; the identical resample is applied to all models.
    # These intervals condition on the fixed OOF fits, as in prior analyses.
    rng = np.random.default_rng(SEED)
    boot = np.empty((args.bootstrap_replicates, len(model_names), 2))
    accepted = 0
    for rep in range(args.bootstrap_replicates):
        ii = rng.integers(len(y), size=len(y))
        if len(np.unique(y[ii])) != 2:
            continue
        for j, model in enumerate(model_names):
            boot[accepted, j, 0] = roc_auc_score(y[ii], pp[model][ii])
            boot[accepted, j, 1] = average_precision_score(y[ii], pp[model][ii])
        accepted += 1
        if (rep + 1) % 1000 == 0:
            print(f'Paired patient bootstrap: {rep + 1}/{args.bootstrap_replicates}', flush=True)
    boot = boot[:accepted]
    pooled = []
    for j, model in enumerate(model_names):
        for k, (metric, fn) in enumerate([('auroc', roc_auc_score), ('average_precision', average_precision_score)]):
            lo, hi = np.quantile(boot[:, j, k], [.025, .975])
            pooled.append({'model': model, 'metric': metric, 'estimate': fn(y, pp[model]),
                           'ci_low': lo, 'ci_high': hi, 'n': len(y), 'pcr_n': int(y.sum()),
                           'bootstrap_replicates': accepted, 'arm_encoding': 'literal' if model.endswith('_literal') else 'harmonized'})
    save(pooled, 'pooled_metrics')
    contrasts = [
        ('hr_her2_plus_five_gene', 'hr_her2_only'),
        ('hr_her2_plus_arm_plus_five_gene', 'hr_her2_plus_arm'),
        ('hr_her2_plus_PFE', 'hr_her2_only'),
        ('hr_her2_plus_arm_plus_PFE', 'hr_her2_plus_arm'),
        ('hr_her2_plus_arm', 'hr_her2_only'),
        ('hr_her2_plus_arm_plus_five_gene', 'hr_her2_plus_five_gene'),
        ('hr_her2_plus_arm_plus_PFE', 'hr_her2_plus_PFE'),
        ('hr_her2_plus_arm_plus_five_gene', 'hr_her2_plus_arm_plus_PFE'),
        ('hr_her2_plus_arm_plus_five_gene_literal', 'hr_her2_plus_arm_literal'),
        ('hr_her2_plus_arm_plus_PFE_literal', 'hr_her2_plus_arm_literal'),
    ]
    deltas = []
    for model, baseline in contrasts:
        a, b = model_names.index(model), model_names.index(baseline)
        for k, (metric, fn) in enumerate([('auroc', roc_auc_score), ('average_precision', average_precision_score)]):
            lo, hi = np.quantile(boot[:, a, k] - boot[:, b, k], [.025, .975])
            deltas.append({'model': model, 'baseline': baseline, 'metric': metric,
                           'delta': fn(y, pp[model]) - fn(y, pp[baseline]),
                           'ci_low': lo, 'ci_high': hi, 'n': len(y),
                           'bootstrap_replicates': accepted})
    save(deltas, 'paired_bootstrap_increments')
    foldframe = pd.DataFrame(foldmetrics)
    fold_deltas = []
    for model, baseline in contrasts:
        a = foldframe[foldframe.model.eq(model)].set_index('fold')
        b = foldframe[foldframe.model.eq(baseline)].set_index('fold')
        for fold in a.index:
            fold_deltas.append({'model': model, 'baseline': baseline, 'fold': fold,
                                'auroc_delta': a.loc[fold, 'auroc'] - b.loc[fold, 'auroc'],
                                'average_precision_delta': a.loc[fold, 'average_precision'] - b.loc[fold, 'average_precision']})
    save(fold_deltas, 'fold_increments')
    checks += [{'check': 'patient_identifiers_unique', 'passed': True},
               {'check': 'all_model_predictions_finite', 'passed': True},
               {'check': 'no_heldout_patient_in_training', 'passed': True},
               {'check': 'all_original_model_predictions_reproduced', 'passed': True}]
    save(checks, 'validation_checks')
    paths = [Path(__file__), reference_path, expression_path, embedding_path, axes_path, weights_path,
             HERE / 'revision_analysis.py'] + [args.reference_dir / f'ispy2_{col}_metadata.tsv' for col in ['arm', 'hr', 'her2']]
    manifest = {
        'analysis': 'I-SPY2 held-out-patient treatment-arm adjustment',
        'n_patients': len(y), 'pcr_n': int(y.sum()), 'n_source_rows': len(source_meta),
        'n_source_patients': source_meta.patient_id.nunique(),
        'duplicate_source_patient': '629606: two platform records; preserve pre-existing first-record selection',
        'missing_arm_n': int(audit.arm.eq('unknown').sum()),
        'missing_hr_n': int(audit.hr.eq('unknown').sum()), 'missing_her2_n': int(audit.her2.eq('unknown').sum()),
        'arm_aliases': ARM_ALIASES, 'arm_alias_changed_n': int(audit.arm_alias_changed.sum()),
        'raw_arm_count': audit.arm_literal.nunique(), 'harmonized_arm_count': audit.arm.nunique(),
        'arm_source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194040',
        'alias_patient_source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5826057',
        'predictors': {'five_gene': GENES, 'PFE_input_features': selected_here,
                       'PFE_potentials': sorted(weights.target.unique().tolist())},
        'validation': 'Existing stored patient-held-out 5 folds; verified against StratifiedGroupKFold(shuffle=True, random_state=20260713)',
        'preprocessing': 'Clinical one-hot categories, empirical ranks and scalers fit only on outer training patients; unknown categories ignored; no interaction terms',
        'classifier': 'LogisticRegression(C=0.2, max_iter=2000); identical default solver and L2 penalty to original',
        'PFE': 'Full 8-dimensional encoder refit in each outer training fold, exact clinical() architecture, optimization, 80 epochs and seeds; no stacking of saved OOF probabilities',
        'bootstrap': {'replicates': accepted, 'seed': SEED, 'unit': 'patient',
                      'ci': 'paired percentile 95%; predictions held fixed; not full retraining uncertainty'},
        'versions': {'python': platform.python_version(), 'numpy': np.__version__, 'pandas': pd.__version__,
                     'scipy': scipy.__version__, 'sklearn': sklearn.__version__, 'h5py': h5py.__version__, 'torch': torch.__version__},
        'inputs': [{'path': str(p.resolve()), 'bytes': p.stat().st_size, 'sha256': sha256(p)} for p in paths],
    }
    (args.output_dir / f'{PREFIX}_provenance.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(pd.DataFrame(pooled).query("metric == 'auroc'").to_string(index=False), flush=True)
    print(pd.DataFrame(deltas).to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
