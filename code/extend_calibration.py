"""Complete fractional-target calibration and skill reporting across sources."""
from pathlib import Path
import numpy as np
import pandas as pd
import argparse
parser=argparse.ArgumentParser();parser.add_argument('--project-root',type=Path,required=True);parser.add_argument('--data-root',type=Path,required=True);parser.add_argument('--output-dir',type=Path,required=True);args=parser.parse_args()
P=args.project_root;E=args.data_root;D=args.output_dir;D.mkdir(parents=True,exist_ok=True)
cal=[];dist=[];skills=[]
def add(y,p,dataset,target,validation):
 f=pd.DataFrame({'observed':y,'predicted':p});f['bin']=pd.qcut(f.predicted,10,duplicates='drop')
 if f.bin.isna().all():f['bin']='all'
 for _,v in f.groupby('bin',observed=True):cal.append(dict(dataset=dataset,target=target,validation=validation,n=len(v),observed=v.observed.mean(),predicted=v.predicted.mean()))
 if validation in ['replicate','leave_perturbation']:
  for v in y:dist.append(dict(dataset=dataset,target=target,value=v))
w=pd.read_parquet(P/'tables/wot_future_fate_targets.parquet').query('day == 0').set_index('cell_id')
pred=pd.read_csv(E/'results/wot_baseline_benchmark/replicate_heldout_predictions.tsv',sep='\t').set_index('cell_id');w=w.loc[pred.index]
for c in [c for c in w if c.startswith('p_branch_')]:add(w[c].to_numpy(),pred[c.replace('p_branch_','pred_')].to_numpy(),'Waddington-OT',c.replace('p_branch_',''),'replicate')
out=w[[c for c in w if c.startswith('p_branch_')]].join(pred).reset_index();out.to_csv(D/'wot_day0_targets_predictions.tsv',sep='\t',index=False)
m=pd.read_csv(E/'results/wot_baseline_benchmark/replicate_heldout_metrics.tsv',sep='\t');mm=m.query("branch == 'macro' and metric == 'brier'").groupby('model').value.mean()
for model,value in mm.items():skills.append(dict(dataset='Waddington-OT',validation='replicate',target='macro branch',model=model,brier=value,skill=1-value/mm.training_mean))
for cell,folder in [('K562','replogle_k562_future_fate'),('RPE1','replogle_rpe1_generalization')]:
 p=pd.read_parquet(P/'results'/folder/'predictions.parquet').query("model == 'baseline_plus_perturbation'")
 for (v,t),f in p.groupby(['validation','target']):add(f.observed.to_numpy(),f.predicted.to_numpy(),cell,t.removeprefix('p_'),v)
 m=pd.read_csv(P/'results'/folder/'metrics.tsv',sep='\t').query("metric == 'brier'")
 mm=m.groupby(['validation','target','model']).value.mean().unstack('model')
 for (v,t),row in mm.iterrows():
  for model,value in row.items():skills.append(dict(dataset=cell,validation=v,target=t.removeprefix('p_'),model=model,brier=value,skill=1-value/row.training_mean))
for data,name in [(cal,'additional_calibration'),(dist,'additional_distributions'),(skills,'additional_brier_skill')]:pd.DataFrame(data).to_csv(D/(name+'.tsv'),sep='\t',index=False)
print('Added transport and genetic-perturbation calibration, distributions and skill.')
