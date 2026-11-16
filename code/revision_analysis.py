"""Reproduce reviewer-requested controls from preserved project data.

Run with Python 3.12, numpy, pandas, scipy, scikit-learn, h5py, pyarrow,
torch and statsmodels. Source data are read only; outputs stay beside this file.
"""
from pathlib import Path
import argparse, json, warnings, hashlib
import numpy as np
import pandas as pd
import h5py
from scipy import sparse
from scipy.stats import spearmanr
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.neighbors import KNeighborsRegressor
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut, StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.impute import SimpleImputer
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score, average_precision_score, balanced_accuracy_score

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
EXTERNAL = Path('source_data')
OUT = HERE / 'data'
OUT.mkdir(exist_ok=True)
SEED = 20261116
TARGETS = ['p_apoptosis', 'p_proliferation', 'p_stress', 'external_resistance_quantile']
KEYS = ['cell_line', 'experiment', 'time_hours']
warnings.filterwarnings('ignore', message='An input array is constant')

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
        g = h['obs']; key = g.attrs.get('_index', '_index')
        return pd.DataFrame({k:hread(g[k]) for k in g if k != key}, index=pd.Index(hread(g[key]), name='cell_id'))

def matrix(path, cells):
    with h5py.File(path) as h:
        g=h['obs']; ids=pd.Index(hread(g[g.attrs.get('_index', '_index')]))
        idx=ids.get_indexer(cells); assert (idx>=0).all()
        x=h['X']
        if isinstance(x,h5py.Dataset): return x[:][idx]
        full=sparse.csr_matrix((x['data'][:],x['indices'][:],x['indptr'][:]),shape=tuple(x.attrs['shape']))
        return full[idx].astype(np.float32)

def save(d, name):
    pd.DataFrame(d).to_csv(OUT/(name+'.tsv'),sep='\t',index=False)

def norm(p):
    p=np.clip(p,0,1); s=p.sum(1,keepdims=True)
    return np.divide(p,s,out=np.full_like(p,1/p.shape[1]),where=s>0)

def boot_metric(y,p,g,metric,reps=1000,seed=SEED):
    rng=np.random.default_rng(seed); ug=np.unique(g)
    indices=[np.flatnonzero(g==v) for v in ug]
    vals=[]
    for _ in range(reps):
        ii=np.concatenate([indices[k] for k in rng.integers(len(ug),size=len(ug))])
        try: vals.append(metric(y[ii],p[ii]))
        except ValueError: pass
    return np.quantile(vals,[.025,.975])

def larry():
    metrics=[]; details=[]; depth=[]; boots=[]
    for dataset in ['in_vitro','cytokine']:
        name='larry_'+dataset
        t=pd.read_csv(PROJECT/'tables'/f'{name}_clone_targets.tsv',sep='\t')
        a=pd.read_csv(PROJECT/'tables'/f'{name}_baseline_assignments.tsv',sep='\t')
        path=EXTERNAL/'staged'/name/f'{name}.h5ad'; meta=obs(path)
        cells=pd.Index(a.cell_id.unique()); x=matrix(path,cells)
        pos=pd.Series(np.arange(len(cells)),index=cells)
        branches=[v for v in t if v.startswith('p_branch_')]
        details.append(dict(dataset=dataset,all_profiled_clones=meta.loc[pd.to_numeric(meta.clone_id,errors='coerce').ge(0),'clone_id'].nunique(),eligible_clones=t.clone_id.nunique(),targets=len(t),baseline_days=','.join(map(str,sorted(a.baseline_time.unique())))))
        for target in ['clone_expansion','transition_entropy','p_rare_fate','p_differentiation']:
            depth.append(dict(dataset=dataset,target=target,rho=spearmanr(t.baseline_cells,t[target]).statistic))
        for mode in ['original','day2','day2_undifferentiated']:
            aa=a.copy()
            if mode!='original': aa=aa[aa.baseline_time.eq(2)]
            if mode.endswith('undifferentiated'):
                aa=aa[aa.cell_id.map(meta.cell_type_annotation.astype(str)).eq('Undifferentiated')]
            units=t.set_index('target_id').loc[pd.unique(aa.target_id)].reset_index()
            # Preserve original assignment order, including tied GroupKFold groups.
            rowpos=aa.cell_id.map(pos).to_numpy(); groups=units.clone_id.to_numpy()
            y=units[branches].to_numpy(float)
            conditions=pd.get_dummies(units.future_condition.astype(str),dtype=float).to_numpy()
            earliest=aa[aa.baseline_time.eq(aa.groupby('target_id').baseline_time.transform('min'))]
            annotation=pd.get_dummies(earliest.cell_id.map(meta.cell_type_annotation).astype(str),dtype=float)
            annotation['target_id']=earliest.target_id.to_numpy()
            ann=annotation.groupby('target_id').mean().loc[units.target_id].to_numpy()
            nbase=aa.groupby('target_id').cell_id.nunique().reindex(units.target_id).to_numpy()
            pred={k:np.zeros_like(y) for k in ['ridge','knn','earliest_annotation','training_mean']}
            scalar={}; fold_id=np.zeros(len(units),int)
            for fold,(tr,te) in enumerate(GroupKFold(5).split(units,groups=groups),1):
                train_ids=set(units.iloc[tr].target_id)
                train_cells=np.unique(rowpos[aa.target_id.isin(train_ids)])
                svd=TruncatedSVD(50,random_state=20260713+fold).fit(x[train_cells])
                z=svd.transform(x)
                zrows=pd.DataFrame(z[rowpos]);zrows['target_id']=aa.target_id.to_numpy()
                zz=zrows.groupby('target_id').mean().loc[units.target_id].to_numpy()
                f=np.column_stack([zz,conditions]); fa=np.column_stack([ann,conditions])
                pp=norm(Ridge(alpha=10).fit(f[tr],y[tr]).predict(f[te]))
                knn=make_pipeline(StandardScaler(),KNeighborsRegressor(n_neighbors=15,weights='distance'))
                kp=norm(knn.fit(f[tr],y[tr]).predict(f[te]))
                ap=norm(Ridge(alpha=10).fit(fa[tr],y[tr]).predict(fa[te]))
                mp=np.repeat(y[tr].mean(0,keepdims=True),len(te),0)
                for model,ppp in [('ridge',pp),('knn',kp),('earliest_annotation',ap),('training_mean',mp)]:
                    pred[model][te]=ppp
                    metrics.append(dict(dataset=dataset,mode=mode,fold=fold,model=model,metric='brier',value=np.mean((y[te]-ppp)**2),n=len(te)))
                    metrics.append(dict(dataset=dataset,mode=mode,fold=fold,model=model,metric='accuracy',value=np.mean(y[te].argmax(1)==ppp.argmax(1)),n=len(te)))
                fold_id[te]=fold
                if mode=='original':
                    for target in ['clone_expansion','transition_entropy','p_rare_fate','p_differentiation']:
                        yy=units[target].to_numpy(float); log=target=='clone_expansion'
                        yt=np.log1p(yy) if log else yy
                        nf=np.column_stack([np.log1p(nbase),conditions])
                        for model,ff in [('depth_only',nf),('depth_plus_state',np.column_stack([nf,zz]))]:
                            v=Ridge(alpha=10).fit(ff[tr],yt[tr]).predict(ff[te]);v=np.expm1(v) if log else v
                            if target.startswith('p_'):v=np.clip(v,0,1)
                            for metric,value in [('mae',np.mean(abs(yy[te]-v))),('spearman',spearmanr(yy[te],v).statistic)]:
                                metrics.append(dict(dataset=dataset,mode=mode,fold=fold,model=model,metric=metric,target=target,value=value,n=len(te)))
            for model,pp in pred.items():
                lo,hi=boot_metric(y,pp,groups,lambda yy,v:np.mean((yy-v)**2))
                boots.append(dict(dataset=dataset,mode=mode,model=model,n_targets=len(units),n_clones=len(np.unique(groups)),brier=np.mean((y-pp)**2),ci_low=lo,ci_high=hi,accuracy=np.mean(y.argmax(1)==pp.argmax(1)),skill=1-np.mean((y-pp)**2)/np.mean((y-pred['training_mean'])**2)))
                frame=pd.DataFrame(pp,columns=[c.replace('p_branch_','pred_') for c in branches]);frame['target_id']=units.target_id;frame['clone_id']=groups;frame['fold']=fold_id
                frame.to_parquet(OUT/f'{name}_{mode}_{model}_predictions.parquet',index=False)
            units.to_csv(OUT/f'{name}_{mode}_targets.tsv',sep='\t',index=False)
            print('LARRY',dataset,mode,len(units),flush=True)
        del x
    save(metrics,'larry_comparators');save(details,'larry_counts');save(depth,'larry_depth_correlations');save(boots,'larry_bootstrap')

def mix_design(rows,num,cat):
    return make_pipeline(ColumnTransformer([('num',make_pipeline(SimpleImputer(),StandardScaler()),num),('cat',OneHotEncoder(handle_unknown='ignore',sparse_output=False),cat)]),Ridge(alpha=10))

def mix_eval(rows,label,variants):
    metrics=[]; predictions=[]
    for mode in ['leave_cell_line','leave_perturbation']:
        g=rows.cell_line if mode=='leave_cell_line' else rows.perturbation
        splitter=GroupKFold(5) if mode=='leave_cell_line' else LeaveOneGroupOut()
        for fold,(train,test) in enumerate(splitter.split(rows,groups=g),1):
            for target in TARGETS:
                tr=train[rows.iloc[train][target].notna()];te=test[rows.iloc[test][target].notna()]
                if len(tr)<20 or not len(te):continue
                yy=rows.iloc[te][target].to_numpy(float);baseline=np.repeat(rows.iloc[tr][target].mean(),len(te))
                for model,(num,cat) in variants.items():
                    pp=baseline if model=='training_mean' else np.clip(mix_design(rows,num,cat).fit(rows.iloc[tr],rows.iloc[tr][target]).predict(rows.iloc[te]),0,1)
                    mse=np.mean((yy-pp)**2);bm=np.mean((yy-baseline)**2)
                    metrics.append(dict(analysis=label,validation=mode,fold=fold,target=target,model=model,n=len(te),brier=mse,mae=np.mean(abs(yy-pp)),spearman=spearmanr(yy,pp).statistic if np.std(pp)>1e-12 else np.nan,skill=1-mse/bm if bm else np.nan))
                    for k,i in enumerate(te): predictions.append(dict(analysis=label,validation=mode,fold=fold,target=target,model=model,pair_id=rows.iloc[i].pair_id,cell_line=rows.iloc[i].cell_line,perturbation=rows.iloc[i].perturbation,observed=yy[k],predicted=pp[k],training_mean=baseline[k]))
    return metrics,predictions

def mixseq():
    rows=pd.read_parquet(EXTERNAL/'results/mixseq_future_fate_benchmark_all/modeling_table.parquet')
    info=pd.read_csv(EXTERNAL/'processed/mixseq/all_CL_features.tsv',sep='\t')
    tissue=info[['CCLE_ID','Disease']].dropna(subset=['CCLE_ID']).drop_duplicates('CCLE_ID').set_index('CCLE_ID').Disease
    rows['tissue']=rows.cell_line.map(tissue).fillna('unknown')
    # A CCLE name suffix is an available tissue annotation when a Disease row is absent.
    missing=rows.tissue.eq('unknown'); rows.loc[missing,'tissue']=rows.loc[missing,'cell_line'].str.split('_',n=1).str[1].fillna('unknown')
    means=[c for c in rows if c.startswith('z') and c.endswith('_mean')]
    full=[c for c in rows if c.startswith('z')]
    variants={'training_mean':([],[]),'mean_only':(means+['time_hours'],['perturbation']),'mean_spread':(full+['time_hours'],['perturbation']),'perturbation_only':(['time_hours'],['perturbation']),'cell_line_only':(['time_hours'],['cell_line']),'tissue_plus_perturbation':(['time_hours'],['tissue','perturbation'])}
    m,p=mix_eval(rows,'original',variants)
    # Independent halves of the same selected control population, repeated five times.
    path=EXTERNAL/'processed/mixseq_baseline_all/baseline_scvi.h5ad'
    meta=obs(path)
    with h5py.File(path) as h:z=h['obsm']['X_baseline_scvi'][:]
    score=pd.read_parquet(EXTERNAL/'processed/mixseq/mixseq_fate_scores_all.parquet')
    match=meta.reset_index()[['cell_id','batch']].copy()
    match['cell_id']=match.cell_id.str.split(':').str[-1]
    cellscore=match.merge(score,on=['cell_id','batch'],how='left',validate='one_to_one')
    assert cellscore.score_apoptosis.notna().all()
    key=lambda v:tuple('NA' if pd.isna(x) else str(x) for x in v)
    controls={key(k):np.asarray(ii) for k,ii in meta.groupby(KEYS,observed=True,dropna=False).indices.items()}
    treated={key(k):d for k,d in score.groupby(KEYS+['perturbation'],observed=True,dropna=False)}
    for repeat in range(5):
        rng=np.random.default_rng(SEED+repeat);features={}
        for kk,ii in controls.items():
            jj=rng.permutation(ii);a,b=jj[:len(jj)//2],jj[len(jj)//2:]
            if min(len(a),len(b))<10:continue
            dd={}
            for j in range(z.shape[1]):
                dd.update({f'z{j}_mean':z[a,j].mean(),f'z{j}_sd':z[a,j].std(),f'z{j}_q10':np.quantile(z[a,j],.1),f'z{j}_q90':np.quantile(z[a,j],.9)})
            for axis in ['apoptosis','proliferation','stress']:
                v=cellscore.iloc[b]['score_'+axis].dropna().to_numpy();med=np.median(v)
                dd[axis+'_threshold']=med+1.4826*np.median(abs(v-med))
            features[kk]=dd
        split=[]
        for _,r in rows.iterrows():
            kk=key(r[KEYS]);ss=features.get(kk)
            if ss is None:continue
            v=r.to_dict();v.update(ss);td=treated[key(r[KEYS+['perturbation']])]
            for axis in ['apoptosis','proliferation','stress']:
                scores=td['score_'+axis].dropna().to_numpy();v['p_'+axis]=((scores>ss[axis+'_threshold']).sum()+1)/(len(scores)+2)
            split.append(v)
        sr=pd.DataFrame(split)
        # Same rows and folds for independent-half and all-control references.
        rr=rows[rows.pair_id.isin(sr.pair_id)].copy().reset_index(drop=True)
        sr=sr.set_index('pair_id').loc[rr.pair_id].reset_index()
        for label,df in [(f'split_{repeat}',sr),(f'matched_full_{repeat}',rr)]:
            mm,pp=mix_eval(df,label,{k:v for k,v in variants.items() if k in ['training_mean','mean_spread','mean_only']});m+=mm;p+=pp
        print('MIX control split',repeat,len(sr),flush=True)
    save(m,'mixseq_comparators');pd.DataFrame(p).to_parquet(OUT/'mixseq_predictions.parquet',index=False)
    save(rows.groupby('tissue').cell_line.nunique().reset_index(name='cell_lines'),'mixseq_tissues')

def axes():
    aa=pd.read_parquet(EXTERNAL/'processed/axes/mixseq_baseline_axes_all.parquet');features=[c for c in aa if '__' in c]
    base=aa.groupby(KEYS,observed=True,dropna=False)[features].mean().reset_index();base[features]=base[features].rank(pct=True)
    targets=pd.read_csv(PROJECT/'tables/mixseq_fate_targets_all_with_resistance.tsv',sep='\t');rows=targets.merge(base,on=KEYS,validate='many_to_one')
    rng=np.random.default_rng(SEED); output=[]
    # Permute outcomes within perturbation/time strata, retaining the same
    # features, overlapping training folds and regularization for every null fit.
    for target in TARGETS:
        f=rows[rows[target].notna()].copy().reset_index(drop=True)
        y=f[target]-f.groupby(['perturbation','time_hours'],dropna=False)[target].transform('mean')
        xx=f[features].to_numpy(float);folds=list(GroupKFold(5).split(f,groups=f.cell_line))
        fitted=[]
        for tr,te in folds:
            sx=StandardScaler().fit_transform(xx[tr]);
            # Precompute the linear ridge operator once for all null fits.
            centered=sx-sx.mean(0)
            op=np.linalg.solve(centered.T@centered+100*np.eye(sx.shape[1]),centered.T)
            fitted.append((tr,op))
        strata=list(f.groupby(['perturbation','time_hours'],dropna=False).indices.values())
        # Duplicate response contexts retain their row-level pairing within each stratum.
        ys=[y.to_numpy()]
        for rep in range(499):
            yy=y.to_numpy().copy()
            for ii in strata: yy[ii]=rng.permutation(yy[ii])
            ys.append(yy)
        for rep,yy in enumerate(ys):
            weights=np.stack([op@(yy[tr]-yy[tr].mean()) for tr,op in fitted])
            stable=((weights>0).all(0)|(weights<0).all(0)).sum()
            output.append(dict(target=target,replicate=rep-1,stable_axes=int(stable),total_axes=len(features)))
        print('AXES',target,output[-500]['stable_axes'],flush=True)
    null=pd.DataFrame(output);save(null,'axis_permutation')
    summary=[]
    for target,f in null.groupby('target'):
        observed=f[f.replicate.eq(-1)].stable_axes.iloc[0];nn=f[f.replicate.ge(0)].stable_axes
        summary.append(dict(target=target,observed=observed,null_mean=nn.mean(),null_low=nn.quantile(.025),null_high=nn.quantile(.975),pvalue=(1+(nn>=observed).sum())/500))
    from statsmodels.stats.multitest import multipletests
    ss=pd.DataFrame(summary);ss['qvalue']=multipletests(ss.pvalue,method='fdr_bh')[1];save(ss,'axis_permutation_summary')

def rank_fit(a,b):
    # Empirical CDF with midranks, fitted solely on training patients.
    aa=np.empty_like(a,dtype=float);bb=np.empty_like(b,dtype=float)
    for j in range(a.shape[1]):
        ref=np.sort(a[np.isfinite(a[:,j]),j]);n=len(ref)
        for v,dest in [(a[:,j],aa[:,j]),(b[:,j],bb[:,j])]:
            dest[:]=(np.searchsorted(ref,v,'left')+np.searchsorted(ref,v,'right')+1)/(2*max(n,1));dest[~np.isfinite(v)]=.5
    return aa,bb

def clinical():
    import torch
    from torch import nn
    from statsmodels.duration.hazard_regression import PHReg
    torch.set_num_threads(1)
    class Encoder(nn.Module):
        def __init__(self,n):
            super().__init__();self.network=nn.Sequential(nn.Linear(n,64),nn.ReLU(),nn.Dropout(.15),nn.Linear(64,8))
        def forward(self,x): return nn.functional.normalize(self.network(x),dim=1)
    def fitenc(x,y,seed):
        torch.manual_seed(seed);rng=np.random.default_rng(seed);net=Encoder(x.shape[1]);opt=torch.optim.AdamW(net.parameters(),lr=.002,weight_decay=.001)
        for _ in range(80):
            order=rng.permutation(len(x))
            for start in range(0,len(x),128):
                ii=order[start:start+128];batch=torch.tensor(x[ii],dtype=torch.float32);lab=torch.tensor(np.tile(y[ii],2));v=torch.cat([batch+.05*torch.randn_like(batch),batch+.05*torch.randn_like(batch)]);z=net(v);sim=z@z.T/.15;sim=sim-sim.max(1,keepdim=True).values.detach();eye=torch.eye(len(z),dtype=torch.bool);positive=lab[:,None].eq(lab[None,:])&~eye;lp=sim-torch.exp(sim).masked_fill(eye,0).sum(1).clamp_min(1e-12).log()[:,None];loss=-(lp*positive).sum(1).div(positive.sum(1).clamp_min(1)).mean();opt.zero_grad();loss.backward();opt.step()
        net.eval();return net
    weights=pd.read_csv(PROJECT/'results/intrinsic_fate_axis_mapper_stability/intrinsic_fate_axis_weights.tsv',sep='\t')
    selected=weights.assign(ab=weights.coefficient.abs()).sort_values('ab',ascending=False).groupby('target').head(20).feature.unique().tolist()
    allfeat=weights.feature.unique().tolist();out=[];preds=[];survival=[];annotations=[]
    for cohort in ['ispy2','beataml']:
        path=EXTERNAL/('staged/gse194040/ispy2_pretreatment_expression.h5ad' if cohort=='ispy2' else 'staged/aml_ohsu_2022/beataml_multiomics.h5ad')
        meta=obs(path);old=pd.read_parquet(EXTERNAL/f'results/patient_fate_embedding/{cohort}/patient_fate_embedding.parquet').set_index('cell_id')
        meta=meta.loc[old.index].copy()
        if cohort=='beataml':
            keep=meta.DISEASE_STAGE_AT_SPECIMEN_COLLECTION.eq('Initial Diagnosis') & meta.DX_AT_INCLUSION.eq('Acute Myeloid Leukaemia (AML) And Related Precursor Neoplasms')
            old=old.loc[keep];meta=meta.loc[keep]
        y=old.outcome.to_numpy(int);groups=old.patient.to_numpy()
        aa=pd.read_parquet(EXTERNAL/f'processed/axes/{cohort}_axes.parquet').set_index('cell_id').loc[old.index]
        common=[f for f in allfeat if f in aa];selected_here=[f for f in selected if f in common]
        raw=aa[common].to_numpy(float)
        with h5py.File(path) as h:
            var=h['var'];key=var.attrs.get('_index','_index');genes=pd.Index(hread(var[key]));gene_ids=genes.get_indexer(['BAX','BCL2','MKI67','PCNA','ATF3']);assert (gene_ids>=0).all()
        gx=matrix(path,old.index)[:,gene_ids]
        # Clinical metadata are transformed within folds using the same patients.
        cats=['hr','her2'] if cohort=='ispy2' else ['ELN2017']
        nums=[] if cohort=='ispy2' else ['AGE_AT_DIAGNOSIS']
        cdata=meta[cats+nums].copy()
        for c in cats:cdata[c]=cdata[c].astype(str).replace({'nan':'unknown','<NA>':'unknown'})
        for c in nums:cdata[c]=pd.to_numeric(cdata[c],errors='coerce')
        models=['PFE','clinical','clinical_plus_PFE','five_gene','clinical_plus_five_gene']
        pp={k:np.full(len(y),np.nan) for k in models};foldid=np.zeros(len(y),int)
        splitter=StratifiedGroupKFold(5,shuffle=True,random_state=20260713)
        for fold,(tr,te) in enumerate(splitter.split(raw,y,groups),1):
            trainrank,testrank=rank_fit(raw[tr],raw[te]);features=[]
            for rr in [trainrank,testrank]:
                potential=[]
                for target,ww in weights.groupby('target',sort=True):
                    ww=ww[ww.feature.isin(common)];ii=[common.index(v) for v in ww.feature];potential.append(((rr[:,ii]-ww.rank_reference_mean.to_numpy())/ww.rank_reference_scale.to_numpy())@ww.coefficient.to_numpy()+ww.intercept.iloc[0])
                features.append(np.column_stack([rr[:,[common.index(v) for v in selected_here]],*potential]))
            # Rank-derived potentials are ranked using training references as well.
            tx,vx=rank_fit(*features);sc=StandardScaler().fit(tx);tx,vx=sc.transform(tx),sc.transform(vx)
            encoder=fitenc(tx,y[tr],20260713+fold)
            with torch.no_grad():tz=encoder(torch.tensor(tx,dtype=torch.float32)).numpy();vz=encoder(torch.tensor(vx,dtype=torch.float32)).numpy()
            parts=[('cat',OneHotEncoder(handle_unknown='ignore',sparse_output=False),cats)]
            if nums:parts.append(('num',make_pipeline(SimpleImputer(add_indicator=True),StandardScaler()),nums))
            ct=ColumnTransformer(parts);tc,vc=ct.fit_transform(cdata.iloc[tr]),ct.transform(cdata.iloc[te])
            gt,gv=rank_fit(gx[tr],gx[te]);sg=StandardScaler().fit(gt);gt,gv=sg.transform(gt),sg.transform(gv)
            designs={'PFE':(tz,vz),'clinical':(tc,vc),'clinical_plus_PFE':(np.column_stack([tc,tz]),np.column_stack([vc,vz])),'five_gene':(gt,gv),'clinical_plus_five_gene':(np.column_stack([tc,gt]),np.column_stack([vc,gv]))}
            for model,(xt,xv) in designs.items():
                prob=LogisticRegression(C=.2,max_iter=2000).fit(xt,y[tr]).predict_proba(xv)[:,1];pp[model][te]=prob
                for metric,value in [('auroc',roc_auc_score(y[te],prob)),('average_precision',average_precision_score(y[te],prob)),('balanced_accuracy',balanced_accuracy_score(y[te],prob>=.5))]:out.append(dict(cohort=cohort,model=model,fold=fold,metric=metric,value=value))
            foldid[te]=fold
        for model,pr in pp.items():
            for i in range(len(y)):preds.append(dict(cohort=cohort,patient=groups[i],cell_id=old.index[i],outcome=y[i],fold=foldid[i],model=model,predicted=pr[i]))
            for metric,fn in [('auroc',roc_auc_score),('average_precision',average_precision_score)]:
                lo,hi=boot_metric(y,pr,groups,fn)
                annotations.append(dict(cohort=cohort,model=model,metric=metric,estimate=fn(y,pr),ci_low=lo,ci_high=hi,n=len(y)))
        for col in ['arm','hr','her2','platform_id','tissue','SPECIMEN_TYPE','TYPE_INDUCTION_TX','ELN2017','AGE_AT_DIAGNOSIS']:
            if col in meta:meta[[col]].assign(patient=groups).to_csv(OUT/f'{cohort}_{col}_metadata.tsv',sep='\t',index=False)
        # Inspect old fully supervised display against the labels used to fit it.
        from sklearn.cluster import KMeans
        cl=KMeans(2,n_init=20,random_state=SEED).fit_predict(old[[f'PFE_{i}' for i in range(1,9)]])
        save(pd.crosstab(cl,y).reset_index(),'clinical_display_clusters_'+cohort)
        if cohort=='beataml':
            time=pd.to_numeric(meta.OS_MONTHS,errors='coerce').to_numpy();event=meta.OS_STATUS.astype(str).str.startswith('1:').to_numpy(int);age=pd.to_numeric(meta.AGE_AT_DIAGNOSIS,errors='coerce').to_numpy();eln=meta.ELN2017.astype(str)
            pred=pp['PFE'];sd=(pred-pred.mean())/pred.std(ddof=1)
            for term,v in [('response_coordinate',sd)]:
                valid=np.isfinite(time)&(time>0);fit=PHReg(time[valid],v[valid,None],status=event[valid]).fit(disp=False)
                survival.append(dict(term=term,model='unadjusted',n=valid.sum(),estimate=np.exp(fit.params[0]),ci_low=np.exp(fit.params[0]-1.96*fit.bse[0]),ci_high=np.exp(fit.params[0]+1.96*fit.bse[0]),pvalue=fit.pvalues[0]))
            known=eln.isin(['Favorable','Intermediate','Adverse']).to_numpy();valid=np.isfinite(time)&(time>0)&np.isfinite(age)&known
            if valid.sum()>20:
                design=np.column_stack([sd[valid],(age[valid]-age[valid].mean())/age[valid].std(),pd.get_dummies(eln[valid],drop_first=True,dtype=float).to_numpy()]);fit=PHReg(time[valid],design,status=event[valid]).fit(disp=False)
                survival.append(dict(term='response_coordinate',model='age_ELN_adjusted',n=valid.sum(),estimate=np.exp(fit.params[0]),ci_low=np.exp(fit.params[0]-1.96*fit.bse[0]),ci_high=np.exp(fit.params[0]+1.96*fit.bse[0]),pvalue=fit.pvalues[0]))
        print('CLINICAL',cohort,len(y),flush=True)
    save(out,'clinical_fold_metrics');save(preds,'clinical_predictions');save(annotations,'clinical_bootstrap');save(survival,'clinical_survival')
    from statsmodels.stats.multitest import multipletests
    ss=pd.DataFrame(survival);ss['qvalue']=multipletests(ss.pvalue,method='fdr_bh')[1];save(ss,'clinical_survival')

def replogle():
    output=[]
    for cell,folder in [('K562','replogle_k562_future_fate'),('RPE1','replogle_rpe1_generalization')]:
        d=pd.read_parquet(PROJECT/'results'/folder/'predictions.parquet')
        table=pd.read_parquet(PROJECT/'results'/folder/'modeling_table.parquet')
        for (validation,target),f in d.groupby(['validation','target']):
            wide=f.pivot(index=['pair_id','observed'],columns='model',values='predicted').reset_index()
            split=wide.pair_id.str.split('|')
            groups=split.str[1].to_numpy() if validation=='leave_gem_group' else split.str[-1].str.split('->').str[-1].to_numpy()
            a=wide.baseline_plus_perturbation.to_numpy();b=wide.perturbation_descriptor.to_numpy();y=wide.observed.to_numpy()
            def delta(yy,pp):return spearmanr(yy,pp[:,0]).statistic-spearmanr(yy,pp[:,1]).statistic
            lo,hi=boot_metric(y,np.column_stack([a,b]),groups,delta)
            rng=np.random.default_rng(SEED);ug=np.unique(groups);ii=[np.flatnonzero(groups==v) for v in ug];values=[]
            for _ in range(5000):
                ix=np.concatenate([ii[k] for k in rng.integers(len(ug),size=len(ug))]);values.append(delta(y[ix],np.column_stack([a[ix],b[ix]])))
            family_lo,family_hi=np.quantile(values,[.05/24,1-.05/24])
            output.append(dict(cell_line=cell,validation=validation,target=target,n_groups=len(np.unique(groups)),n_pairs=len(y),delta_spearman=delta(y,np.column_stack([a,b])),ci_low=lo,ci_high=hi,family_ci_low=family_lo,family_ci_high=family_hi))
    save(output,'replogle_group_bootstrap')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--project-root',type=Path,default=PROJECT);parser.add_argument('--data-root',type=Path,default=EXTERNAL);parser.add_argument('--output-dir',type=Path,default=OUT);parser.add_argument('phase',choices=['larry','mixseq','axes','clinical','replogle']);args=parser.parse_args()
    PROJECT=args.project_root;EXTERNAL=args.data_root;OUT=args.output_dir;OUT.mkdir(parents=True,exist_ok=True)
    globals()[args.phase]()
