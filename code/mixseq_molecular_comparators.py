#!/usr/bin/env python3
"""Fixed molecular descriptor controls on preserved MIX-seq chemical pairs.

Offline reproduction (no chemical toolkit/network required):
  OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python \
    reproducibility/mixseq_molecular_comparators.py \
    --modeling-table /path/to/modeling_table.parquet

Descriptor table/provenance are frozen alongside the predictions. The molecular
models never receive perturbation identity. Identity controls are explicitly
separate. Output files are new and the script refuses to overwrite them.
"""
from pathlib import Path
import argparse, base64, hashlib, json, platform, warnings
import numpy as np
import pandas as pd
import scipy, sklearn
from scipy.stats import spearmanr, pearsonr
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler, OneHotEncoder

HERE=Path(__file__).resolve().parent
TARGETS=['p_apoptosis','p_proliferation','p_stress','external_resistance_quantile']
PHYS=['MolecularWeight','XLogP','TPSA','Complexity','HBondDonorCount',
      'HBondAcceptorCount','RotatableBondCount','HeavyAtomCount','Charge']
SEED=20261116
warnings.filterwarnings('ignore', message='An input array is constant')

def fingerprint(value):
    if pd.isna(value):return np.full(881,np.nan)
    b=base64.b64decode(value,validate=True)
    n=int.from_bytes(b[:4],'big')
    assert n==881 and len(b)==115, (n,len(b))
    return np.unpackbits(np.frombuffer(b[4:],dtype=np.uint8),bitorder='big')[:n].astype(float)

def design(num,cat):
    parts=[('num',make_pipeline(SimpleImputer(strategy='mean'),StandardScaler()),num)]
    if cat:parts.append(('cat',OneHotEncoder(handle_unknown='ignore',sparse_output=False),cat))
    return make_pipeline(ColumnTransformer(parts),Ridge(alpha=10))

def correlation(y,p,fn=spearmanr):
    return float(fn(y,p).statistic) if len(y)>2 and np.std(y)>1e-12 and np.std(p)>1e-12 else np.nan

def write_new(df,path):
    if path.exists():raise FileExistsError(f'Refusing to overwrite {path}')
    if path.suffix=='.parquet':df.to_parquet(path,index=False)
    else:df.to_csv(path,sep='\t',index=False)

def bootstrap_weights(group,replicates,seed):
    unique,inv=np.unique(group,return_inverse=True)
    rng=np.random.default_rng(seed)
    weights=rng.multinomial(len(unique),np.full(len(unique),1/len(unique)),size=replicates)
    counts=np.bincount(inv)
    return inv,weights,weights@counts

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--modeling-table',type=Path,default=HERE.parent/'inputs/mixseq_modeling_table.parquet')
    ap.add_argument('--descriptor-table',type=Path,default=HERE.parent/'results/mixseq_molecular_descriptors.tsv')
    ap.add_argument('--output-dir',type=Path,default=HERE.parent/'reproduced/molecular')
    ap.add_argument('--prefix',default='mixseq_molecular')
    ap.add_argument('--bootstrap-replicates',type=int,default=2000)
    ap.add_argument('--reference-predictions',type=Path,default=HERE.parent/'results/mixseq_predictions.parquet')
    args=ap.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
    def output(kind,suffix='tsv'):return args.output_dir/f'{args.prefix}_{kind}.{suffix}'
    kinds=['fold_metrics','summary','paired_comparisons','model_correlations','coverage','fold_assignments']
    for p in [output(k) for k in kinds]+[output('predictions','parquet'),output('run','json')]:
        if p.exists():raise FileExistsError(f'Refusing to overwrite {p}')
    rows=pd.read_parquet(args.modeling_table).reset_index(drop=True)
    assert rows.pair_id.is_unique and len(rows)==1044
    desc=pd.read_csv(args.descriptor_table,sep='\t').set_index('perturbation')
    assert set(rows.perturbation)==set(desc.index)
    fp=np.stack([fingerprint(v) for v in desc.Fingerprint2D])
    fpcols=[f'fp_{i:03d}' for i in range(881)]
    mol=pd.concat([desc[['structure_available']+PHYS],pd.DataFrame(fp,index=desc.index,columns=fpcols)],axis=1)
    rows=rows.join(mol,on='perturbation',validate='many_to_one')
    assert rows.structure_available.dtype==bool
    means=[c for c in rows if c.startswith('z') and c.endswith('_mean')]
    state=[c for c in rows if c.startswith('z')]
    assert len(means)==30 and len(state)==120
    time=['time_hours']
    models={
      'training_mean':([],[]),
      'time_only':(time,[]),
      'categorical_identity_time':(time,['perturbation']),
      'categorical_identity_mean_spread_time':(time+state,['perturbation']),
      'structure_fp_time':(time+fpcols,[]),
      'structure_fp_mean_time':(time+fpcols+means,[]),
      'structure_fp_mean_spread_time':(time+fpcols+state,[]),
      'physicochemical_time':(time+PHYS,[]),
      'physicochemical_mean_spread_time':(time+PHYS+state,[])}
    fold_records=[];folds={}
    for mode in ['leave_cell_line','leave_perturbation']:
        groups=rows.cell_line if mode=='leave_cell_line' else rows.perturbation
        splitter=GroupKFold(5) if mode=='leave_cell_line' else LeaveOneGroupOut()
        folds[mode]=list(splitter.split(rows,groups=groups))
        for f,(_,te) in enumerate(folds[mode],1):
            for i in te:fold_records.append(dict(validation=mode,fold=f,pair_id=rows.iloc[i].pair_id,cell_line=rows.iloc[i].cell_line,perturbation=rows.iloc[i].perturbation))
    fold_table=pd.DataFrame(fold_records)
    reference_match=None
    if args.reference_predictions.is_file():
        reference=pd.read_parquet(args.reference_predictions,filters=[('analysis','==','original')],columns=['validation','fold','pair_id']).drop_duplicates()
        a=fold_table[['validation','fold','pair_id']].sort_values(['validation','pair_id']).reset_index(drop=True)
        b=reference.sort_values(['validation','pair_id']).reset_index(drop=True)
        reference_match=a.equals(b)
        assert reference_match,'Outer folds differ from original revision_analysis.py predictions'
    covered=rows.structure_available.to_numpy()
    coverage=rows.groupby('perturbation',sort=True).agg(n_pairs=('pair_id','size'),structure_available=('structure_available','first'),n_cell_lines=('cell_line','nunique'),n_resistance=('external_resistance_quantile','count'),n_time_missing=('time_hours',lambda x:x.isna().sum())).reset_index()
    coverage=coverage.merge(desc.reset_index()[['perturbation','CID','lookup_name','full_source_compound_id','source_url']],on='perturbation',validate='one_to_one')
    predictions=[];fm=[]
    # The covered-only sensitivity keeps the original all-pair fold assignment;
    # no GroupKFold is regenerated after dropping the unmapped compound.
    for scope,eligible in [('all_pairs_partial_structure_coverage',np.ones(len(rows),bool)),('structure_covered_pairs',covered)]:
        for mode,split in folds.items():
            for f,(train,test) in enumerate(split,1):
                for target in TARGETS:
                    valid=rows[target].notna().to_numpy()&eligible
                    tr=train[valid[train]];te=test[valid[test]]
                    if len(tr)<20 or not len(te):continue
                    yy=rows.iloc[te][target].to_numpy(float)
                    base=np.repeat(rows.iloc[tr][target].mean(),len(te))
                    for model,(num,cat) in models.items():
                        pred=base if model=='training_mean' else np.clip(design(num,cat).fit(rows.iloc[tr],rows.iloc[tr][target]).predict(rows.iloc[te]),0,1)
                        err=(yy-pred)**2;be=(yy-base)**2
                        fm.append(dict(scope=scope,validation=mode,fold=f,target=target,model=model,n=len(te),brier=err.mean(),mae=np.abs(yy-pred).mean(),spearman=correlation(yy,pred),skill=1-err.mean()/be.mean()))
                        block=rows.iloc[te][['pair_id','cell_line','perturbation','time_hours','structure_available']].copy()
                        block=block.assign(scope=scope,validation=mode,fold=f,target=target,model=model,observed=yy,predicted=pred,training_mean=base)
                        predictions.append(block)
                print(scope,mode,f,'complete',flush=True)
    p=pd.concat(predictions,ignore_index=True);fm=pd.DataFrame(fm)
    summaries=[];comparisons=[];cors=[]
    contrast=[('structure_fp_time','structure_fp_mean_time'),('structure_fp_time','structure_fp_mean_spread_time'),('structure_fp_mean_time','structure_fp_mean_spread_time'),('categorical_identity_time','structure_fp_time'),('categorical_identity_mean_spread_time','structure_fp_mean_spread_time'),('physicochemical_time','physicochemical_mean_spread_time'),('categorical_identity_time','categorical_identity_mean_spread_time')]
    for (scope,mode,target),d in p.groupby(['scope','validation','target'],sort=True):
        wide=d.pivot(index='pair_id',columns='model',values='predicted')
        u=d.drop_duplicates('pair_id').set_index('pair_id').loc[wide.index]
        y=u.observed.to_numpy();g=u['cell_line' if mode=='leave_cell_line' else 'perturbation'].to_numpy()
        inv,w,den=bootstrap_weights(g,args.bootstrap_replicates,SEED)
        sums={m:np.bincount(inv,weights=(y-wide[m].to_numpy())**2) for m in models}
        baseerror=np.mean((y-u.training_mean.to_numpy())**2)
        for model in models:
            pp=wide[model].to_numpy();brier=np.mean((y-pp)**2)
            lo,hi=np.quantile((w@sums[model])/den,[.025,.975])
            foldm=fm[(fm.scope==scope)&(fm.validation==mode)&(fm.target==target)&(fm.model==model)]
            # Demeaning is diagnostic only and never used in model fitting.
            residuals=pd.DataFrame({'perturbation':u.perturbation.to_numpy(),'time':u.time_hours.to_numpy(),'y':y,'p':pp})
            for v in ['y','p']:residuals[v]-=residuals.groupby(['perturbation','time'],dropna=False)[v].transform('mean')
            summaries.append(dict(scope=scope,validation=mode,target=target,model=model,n=len(y),n_groups=len(np.unique(g)),brier=brier,ci_low=lo,ci_high=hi,mae=np.abs(y-pp).mean(),spearman=correlation(y,pp),pearson=correlation(y,pp,pearsonr),within_perturbation_time_spearman=correlation(residuals.y.to_numpy(),residuals.p.to_numpy()),skill=1-brier/baseerror,fold_mean_brier=foldm.brier.mean(),fold_mean_spearman=foldm.spearman.mean()))
        for base,added in contrast:
            dif=(y-wide[added].to_numpy())**2-(y-wide[base].to_numpy())**2
            bs=(w@(sums[added]-sums[base]))/den
            lo,hi=np.quantile(bs,[.025,.975])
            comparisons.append(dict(scope=scope,validation=mode,target=target,base_model=base,added_model=added,n=len(y),n_groups=len(np.unique(g)),delta_brier=dif.mean(),ci_low=lo,ci_high=hi,relative_brier_reduction=-dif.mean()/np.mean((y-wide[base])**2),fraction_pairwise_improved=np.mean(dif<0),bootstrap_replicates=args.bootstrap_replicates))
        for i,a in enumerate(models):
            for b in list(models)[i+1:]:cors.append(dict(scope=scope,validation=mode,target=target,model_a=a,model_b=b,n=len(y),spearman=correlation(wide[a].to_numpy(),wide[b].to_numpy()),pearson=correlation(wide[a].to_numpy(),wide[b].to_numpy(),pearsonr)))
    tables={'fold_metrics':fm,'summary':pd.DataFrame(summaries),'paired_comparisons':pd.DataFrame(comparisons),'model_correlations':pd.DataFrame(cors),'coverage':coverage,'fold_assignments':fold_table}
    for k,v in tables.items():write_new(v,output(k))
    write_new(p,output('predictions','parquet'))
    run={'seed':SEED,'modeling_table_path':str(args.modeling_table),'modeling_table_sha256':hashlib.sha256(args.modeling_table.read_bytes()).hexdigest(),'descriptor_table_sha256':hashlib.sha256(args.descriptor_table.read_bytes()).hexdigest(),'n_pairs':len(rows),'n_cell_lines':rows.cell_line.nunique(),'n_perturbations':rows.perturbation.nunique(),'n_structure_covered_pairs':int(covered.sum()),'n_structure_covered_perturbations':int(desc.structure_available.sum()),'time_covariate':'Original time_hours unchanged, including missing values. Training-fold mean imputation. Supplementary Table2 says experiment10 is24hours but this is not silently filled.','n_time_missing':int(rows.time_hours.isna().sum()),'outcomes':'Unchanged modeling_table columns. Brier is mean squared error of continuous threshold-defined program fractions; external resistance quantile is a ranked response, not an event probability.','outer_folds':'Original full1044-row GroupKFold5 by cell line; LeaveOneGroupOut by perturbation. Covered subset uses these same assignments, with training/test rows restricted to covered drugs.','matches_original_revision_predictions':reference_match,'ridge_alpha':10,'prediction_bounds':[0,1],'preprocessing':'SimpleImputer(mean) and StandardScaler fit separately on training rows for each target/fold/model. OneHotEncoder only in explicitly named categorical_identity controls. No identity, target label or outcome-based descriptor selection used in structure/physicochemical models.','missing_structure':'BRD3379 is184 pairs with unresolved structure. All881 bits and9 scalar properties missing, imputed from training-fold feature means. No missingness indicator is supplied. Covered860-pair sensitivity excludes BRD3379.','baseline_features':{'mean':means,'mean_spread':state},'models':{k:{'n_numeric':len(v[0]),'categorical':v[1]} for k,v in models.items()},'bootstrap':'2000 resamples of heldout grouping units, applied to fixed OOF errors (no refitting). Percentile95% CI; descriptive uncertainty,13drugs for leave-perturbation.','versions':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'sklearn':sklearn.__version__}}
    output('run','json').write_text(json.dumps(run,indent=2)+'\n')
    print(pd.DataFrame(summaries).to_string(index=False),flush=True)
if __name__=='__main__':main()
