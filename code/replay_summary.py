from pathlib import Path
import zipfile, io
import pandas as pd
r=Path(__file__).resolve().parents[1]
with zipfile.ZipFile(r/'derived_data.zip') as z:
 def table(n):return pd.read_csv(io.BytesIO(z.read('results/'+n+'.tsv')),sep='\t')
 x=table('larry_bootstrap');print(x[(x['mode']=='original')&(x['model']=='ridge')][['dataset','n_clones','brier','ci_low','ci_high']].to_string(index=False))
 x=table('axis_permutation_summary');print(x[['target','observed','null_mean','qvalue']].to_string(index=False))
 x=table('clinical_bootstrap');print(x[x.metric.eq('auroc')][['cohort','model','estimate','ci_low','ci_high']].to_string(index=False))
 x=table('clinical_survival');print(x.to_string(index=False))
 assert set(x['n'])=={317,303}
