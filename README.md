# Baseline Transcriptomes Encode Context-Dependent Perturbation Outcomes

Code and compact derived data for **Baseline transcriptomic states encode context-dependent information about lineage and perturbation outcomes across single-cell systems**.

The release contains the analysis code, processed outcome targets, out-of-fold predictions, comparator results, permutation results and numerical figure sources. It excludes manuscripts, peer-review correspondence, raw patient expression, raw single-cell expression, trained neural-network checkpoints and redundant images. Patient identifiers in released prediction rows are replaced by hashes; source datasets remain available from their original studies.

## Contents

- `code/revision_analysis.py`: additional matched comparators, independent control-half analysis, axis permutation null, early-progenitor and clone-depth controls, clinical models with training-fold ranks, and grouped bootstrap intervals.
- `source_workflow.zip`: original data staging, target construction, scoring, model and plotting scripts, together with source download manifests and gene-program definitions. These scripts document the original analyses; corrected clinical processing is in `code/revision_analysis.py`.
- `derived_data.zip`: processed targets and predictions, complete numerical controls, and figure source tables. It also contains the small MIX-seq modeling table for comparator replay.
- `code/replay_summary.py`: a lightweight check of the principal released estimates without downloading expression matrices.

## Quick numerical replay

```bash
python -m pip install -r requirements.txt
python code/replay_summary.py
```

This checks the released prediction/target results. It does not retrain encoders or substitute for an independent reanalysis.

## Figures from the compact numerical sources

Unpack `derived_data.zip` in this directory, then run `Rscript code/make_figures.R`. Required R packages are ggplot2, patchwork, dplyr, tidyr, scales, ggrepel, svglite and ragg. Outputs go to `output/`. This recreates the figures from the released numerical sources.

## Full reanalysis

Unpack `source_workflow.zip` and obtain source data using `config/download_manifest.tsv` and the data-staging scripts. Keep the original staged/processed/results layout documented by those scripts. Unpack `derived_data.zip` into the same project directory; its `tables` directory supplies the exact outcome definitions used here. Raw expression matrices are intentionally not duplicated in this release.

```bash
python code/revision_analysis.py --project-root /path/to/project --data-root /path/to/source_data --output-dir /path/to/output larry
```

Use `mixseq`, `axes`, `clinical` or `replogle` instead of `larry` for the other analyses. Full reanalysis additionally requires the source-specific packages in `requirements-core.txt` (including AnnData/scVI for source staging and the original encoder). The saved source encoder can be regenerated using `train_baseline_scvi.py`; the chemical sensitivity analyses remain conditional on that shared encoder.

## Interpretation and uncertainty

LARRY targets are experimentally lineage linked. MIX-seq and Replogle targets are population-level expression programs; Waddington-OT targets are expression-derived transport probabilities. Clone-condition rows, cells and patient observations have different inferential units. Bootstrap intervals hold fitted out-of-fold predictions fixed and do not include full model-retraining uncertainty. Axis nulls retain the same overlapping folds and permute outcomes within perturbation/time strata. The two reported clinical survival tests form one family.

## Source data

- LARRY: https://doi.org/10.1126/science.aaw3381
- MIX-seq: https://doi.org/10.1038/s41467-020-17440-w ; Figshare 10298696
- Waddington-OT: https://doi.org/10.1016/j.cell.2019.01.006
- Replogle Perturb-seq: https://doi.org/10.1016/j.cell.2022.05.013 ; Figshare 20029387
- Xenium: GEO GSE291246
- I-SPY2: GEO GSE194040
- BeatAML: cBioPortal study `aml_ohsu_2022`; https://doi.org/10.1016/j.ccell.2022.07.002

Source data retain their original terms. No raw patient or expression data are redistributed here.
