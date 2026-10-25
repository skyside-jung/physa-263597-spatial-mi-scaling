from pathlib import Path
import zipfile,io,csv,json,math,collections,shutil,hashlib,subprocess,sys
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
base=Path(__file__).resolve().parent; out=base;out.mkdir(exist_ok=True)
for folder in ['figures_recreated','checks','supplementary','source']:(out/folder).mkdir(exist_ok=True)
source=base/'source'/'analysis_only.zip';z=zipfile.ZipFile(source)
cfg=json.loads(z.read('protocol.json'))
rows=list(csv.DictReader(io.StringIO(z.read('merged/profiles.csv').decode('utf-8-sig'))))
original_scores=list(csv.DictReader(io.StringIO(z.read('analysis/collapse_scores.csv').decode('utf-8-sig'))))
by=collections.defaultdict(list)
for x in rows:by[(x['model'],int(x['L']),float(x['T']),int(x['seed']))].append(x)
rhos=np.linspace(.15,.45,13)
def averaged(model,L,T,field,grid=rhos):
 allseeds=[]
 for seed in cfg['seeds']:
  chunk=sorted(by[(model,L,float(T),seed)],key=lambda a:float(a['rho']))
  allseeds.append(np.interp(grid,[float(x['rho']) for x in chunk],[float(x[field]) for x in chunk]))
 return np.mean(allseeds,axis=0)
def collapse(model,sizes,T,field='mi_symmetry_pooled',exponent=None):
 if exponent is None:exponent=cfg['inference']['exponents'][model][field]
 a=np.array([averaged(model,L,T,field)*L**exponent for L in sizes]);return float((np.std(a,axis=0)/np.mean(a,axis=0)).mean())
calcs=[];reports=[]
for model in cfg['models']:
 for group,sizes in [('development',cfg['inference']['development_sizes']),('heldout',cfg['inference']['heldout_sizes'])]:
  for field in cfg['inference']['estimators']:
   for T in cfg['temperatures'][model]:
    score=collapse(model,sizes,T,field)
    old=next(a for a in original_scores if a['model']==model and a['group']==group and a['estimator']==field and abs(float(a['T'])-T)<1e-8)
    calcs.append(dict(model=model,group=group,estimator=field,T=T,score=score,reference=float(old['score']),abs_delta=abs(score-float(old['score']))))
  subset=[r for r in calcs if r['model']==model and r['group']==group and r['estimator']=='mi_symmetry_pooled']
  best=min(subset,key=lambda x:x['score']);reports.append(dict(model=model,group=group,T_hat=best['T'],score=best['score']))
print('max score mismatch',max(a['abs_delta'] for a in calcs))
assert len(calcs)==(6+7)*2*2 and max(a['abs_delta'] for a in calcs)<1e-10
with open(out/'checks'/'recomputed_collapse_scores.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=calcs[0]);w.writeheader();w.writerows(calcs)
plt.rcParams.update({'font.size':9,'figure.dpi':150})
def figure_model(model, groups, name):
 if model=='ising': fig,axes=plt.subplots(1,3,figsize=(12,3.6));sets=['development','heldout'];temps=cfg['temperatures'][model];colors=['C0','C1'];
 else: fig,axes=plt.subplots(1,2,figsize=(9,3.8));sets=[groups];temps=cfg['temperatures'][model];colors=['C0']
 a=axes[0]
 for group,color in zip(sets,colors):
  ds=[c for c in calcs if c['model']==model and c['group']==group and c['estimator']=='mi_symmetry_pooled']
  ds=sorted(ds,key=lambda x:x['T']);a.plot([x['T'] for x in ds],[x['score'] for x in ds],marker='o',label=group,color=color)
 exact=2/math.log(1+math.sqrt(2)) if model=='ising' else 1/math.log(1+math.sqrt(3));a.axvline(exact,color='gray',ls=':',label='Exact Tc');a.set(xlabel='Temperature T',ylabel='Cross-size collapse score',title='(a) Temperature selection');a.legend(fontsize=7)
 for j,group in enumerate(sets):
  ax=axes[j+1];sizes=cfg['inference']['development_sizes'] if group=='development' else cfg['inference']['heldout_sizes'];T=2.27 if model=='ising' else .995;expo=cfg['inference']['exponents'][model]['mi_symmetry_pooled']
  for L in sizes: ax.plot(rhos,averaged(model,L,T,'mi_symmetry_pooled')*L**expo,label=f'L={L}')
  ax.set(xlabel='Relative distance r/L',ylabel=f'Scaled MI, L^{expo:.4f} I',title=f'({chr(98+j)}) {group.title()} sizes');ax.legend(fontsize=7,ncol=2)
 fig.tight_layout();fig.savefig(out/'figures_recreated'/(name+'.png'),dpi=180);fig.savefig(out/'figures_recreated'/(name+'.pdf'));plt.close(fig)
figure_model('ising',None,'Figure2_Ising_recreated')
figure_model('potts3','development','Figure3_Potts_Development_recreated')
figure_model('potts3','heldout','Figure4_Potts_Validation_recreated')
# Figure 5: per-seed matched data from preserved profiles, point locations from the saved exact grid
m='potts3';t=.995
fig,ax=plt.subplots(1,2,figsize=(10,3.8))
fieldmap={'mi_per_snapshot_mean':'Configuration mean','mi_pooled':'Pooled MI','mi_symmetry_pooled':'Symmetry-pooled MI'}
for field,label in fieldmap.items():
 x=np.linspace(.15,.5,15);ax[0].plot(x,averaged(m,128,t,field,x),label=label)
ax[0].set(xlabel='Relative distance r/L',ylabel='Mutual information (bits)',title='(a) Potts, L=128, T=0.995');ax[0].legend(fontsize=8)
counts=[]
for field,label in fieldmap.items():
 count=0
 for L in cfg['sizes']:
  for seed in cfg['seeds']:
   values=sorted(by[(m,L,t,seed)],key=lambda a:float(a['rho']))
   xs=np.array([float(x['rho']) for x in values]);ys=np.array([float(x[field]) for x in values]);v=np.interp([.3,.5],xs,ys)
   if v[1]-v[0]>0:count+=1
 counts.append(count)
ax[1].bar(list(fieldmap.values()),counts);ax[1].set(ylabel='Count of positive changes',ylim=(0,100),title='(b) I(0.50L)-I(0.30L) > 0');ax[1].tick_params(axis='x',rotation=15)
for i,count in enumerate(counts):ax[1].text(i,count+2,f'{count}/90',ha='center',fontsize=9)
fig.tight_layout();fig.savefig(out/'figures_recreated'/'Figure5_Aggregation_Diagnostic_recreated.png',dpi=180);fig.savefig(out/'figures_recreated'/'Figure5_Aggregation_Diagnostic_recreated.pdf');plt.close(fig)
print('Figure5 delta sign counts',counts)
assert counts==[90,0,0]

result={'score_checks':len(calcs),'max_score_abs_difference':max(a['abs_delta'] for a in calcs),'figure5_counts':counts,'status':'PASS'}
(out/'checks'/'rerun_figure_check.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
