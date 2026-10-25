"""Recalculate Table 1 directly from preserved per-seed distance profiles.
Usage: python recompute_table1.py physa_main_v11_results.zip
Requires numpy only, avoids reading saved collapse scores as computational input.
"""
import csv,io,json,sys,zipfile,collections
import numpy as np

def main(path):
 with zipfile.ZipFile(path) as z:
  base='physa_main_v11/'
  p=json.loads(z.read(base+'protocol.json'))
  profiles=list(csv.DictReader(io.StringIO(z.read(base+'merged/profiles.csv').decode('utf-8-sig'))))
  expected=list(csv.DictReader(io.StringIO(z.read(base+'analysis/selected_temperatures.csv').decode('utf-8-sig'))))
 byseed=collections.defaultdict(list)
 for r in profiles:byseed[(r['model'],int(r['L']),float(r['T']),int(r['seed']))].append(r)
 settings=p['inference'];rho=np.linspace(settings['rho_min'],settings['rho_max'],settings['rho_points'])
 sets={'development':settings['development_sizes'],'heldout':settings['heldout_sizes']}
 out=[]
 for model in p['models']:
  for group,sizes in sets.items():
   for estimator in settings['estimators']:
    exponent=settings['exponents'][model][estimator]; scores=[]
    for T in p['temperatures'][model]:
     curves=[]
     for L in sizes:
      per_seed=[]
      for seed in p['seeds']:
       rows=sorted(byseed[(model,L,float(T),seed)],key=lambda a:float(a['rho']))
       x=[float(r['rho']) for r in rows]; y=[float(r[estimator]) for r in rows]
       per_seed.append(np.interp(rho,x,y))
      curves.append(np.mean(per_seed,axis=0)*L**exponent)
     arr=np.array(curves)
     cv=np.std(arr,axis=0,ddof=0)/np.mean(arr,axis=0)
     scores.append((float(np.mean(cv)),float(T)))
    score,T=min(scores)
    row=next(r for r in expected if (r['model'],r['group'],r['estimator'])==(model,group,estimator))
    delta=abs(score-float(row['min_score']))
    ok=abs(T-float(row['T_hat']))<1e-12 and delta<1e-10
    out.append(dict(model=model,group=group,estimator=estimator,T=T,score=score,score_delta=delta,pass_check=ok))
 status='PASS' if len(out)==8 and all(r['pass_check'] for r in out) else 'FAIL'
 print(json.dumps(dict(status=status,results=out,source='merged/profiles.csv; protocol.json; saved selections used only for comparison'),indent=2))
 return 0 if status=='PASS' else 1
if __name__=='__main__':sys.exit(main(sys.argv[1]))
