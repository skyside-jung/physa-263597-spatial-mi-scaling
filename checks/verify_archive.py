"""Check source ZIP integrity, SHA256 manifest, job completeness and Table 1.
Usage: python verify_archive.py /path/to/physa_main_v11_results.zip
No external Python dependencies.
"""
import csv,hashlib,io,json,sys,zipfile,collections
from pathlib import Path

def main(filename):
    with zipfile.ZipFile(filename) as z:
        root='physa_main_v11/'
        members={x.filename:x for x in z.infolist() if not x.is_dir()}
        manifest=json.loads(z.read(root+'manifest.json'))
        jobs=json.loads(z.read(root+'jobs.json'))
        protocol=json.loads(z.read(root+'protocol.json'))
        errors=[]; sha_checked=0
        for rel,expected in manifest.items():
            name=root+rel.replace('\\','/')
            if name not in members:
                errors.append('Missing manifest entry: '+name);continue
            digest=hashlib.sha256(z.read(name)).hexdigest()
            sha_checked+=1
            if digest.lower()!=expected.lower():errors.append('SHA256 mismatch: '+name)
        for j in jobs:
            folder=f"{j['model']}_L{j['L']}_T{j['T']:g}_seed{j['seed']}"
            for file in ('joint_counts.csv','observables.csv','DONE.json'):
                if root+'jobs/'+folder+'/'+file not in members:
                    errors.append('Missing job file: '+folder+'/'+file)
        def table(name):
            return list(csv.DictReader(io.StringIO(z.read(root+name).decode('utf-8-sig'))))
        scores=table('analysis/collapse_scores.csv')
        selections=table('analysis/selected_temperatures.csv')
        grid=collections.defaultdict(list)
        for s in scores:grid[(s['model'],s['group'],s['estimator'])].append(s)
        results=[]
        for row in selections:
            key=(row['model'],row['group'],row['estimator'])
            if not grid[key]:errors.append('Missing collapse score: '+str(key));continue
            best=min(grid[key],key=lambda a:(float(a['score']),float(a['T'])))
            success=(abs(float(best['T'])-float(row['T_hat']))<1e-12 and abs(float(best['score'])-float(row['min_score']))<1e-11)
            if not success:errors.append('Score selection mismatch: '+str(key))
            results.append({'model':key[0], 'group':key[1], 'estimator':key[2], 'temperature':best['T'],'min_score':best['score'],'matches_saved':success})
        result={'source':Path(filename).name,'zip_member_count':len(members),'manifest_entries':len(manifest), 'sha256_checked':sha_checked,'jobs_declared':len(jobs),'joint_count_files':sum(n.endswith('/joint_counts.csv') for n in members),'observation_files':sum(n.endswith('/observables.csv') for n in members),'completed_markers':sum(n.endswith('/DONE.json') for n in members),'selection_checks':results,'protocol_id':protocol.get('experiment_id'),'errors':errors,'status':'PASS' if not errors else 'FAIL'}
    print(json.dumps(result,ensure_ascii=False,indent=2))
    return 0 if not errors else 1
if __name__=='__main__':
    if len(sys.argv)!=2:raise SystemExit('Usage: python verify_archive.py physa_main_v11_results.zip')
    sys.exit(main(sys.argv[1]))
