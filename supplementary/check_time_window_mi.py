#!/usr/bin/env python3
"""PHYSA-263597: time-window sensitivity from original sample-wise joint counts.
No resimulation. No modification of production files.
Run in the directory containing physa_main_v11_results.zip.
"""
import argparse
import csv
import io
import math
import re
import zipfile
from collections import defaultdict
from pathlib import Path

DEV = (32, 48, 64, 96, 128)
VAL = (72, 104, 144, 208)
RHOS = tuple(0.15 + 0.025 * i for i in range(13))
EXPONENT = {'ising': 0.5, 'potts3': 8.0 / 15.0}
REFERENCE = {('ising', 'development'): (2.270, 0.02370056),
             ('ising', 'validation'): (2.270, 0.03210092),
             ('potts3', 'development'): (0.995, 0.00987459),
             ('potts3', 'validation'): (0.995, 0.00585080)}
WINDOWS = {'full': (0, 500), 'drop100': (100, 500),
           'drop200': (200, 500), 'first250': (0, 250),
           'last250': (250, 500)}
PATTERN = re.compile(r'^(ising|potts3|potts)_L(\d+)_T([0-9.]+)_seed(\d+)$')


def sym_mi(counts, q):
    total = sum(counts)
    if total <= 0:
        raise ValueError('Empty counts')
    p = sum(counts[a*q+a] for a in range(q)) / total
    value = 0.0
    if p > 0:
        value += p * math.log2(q*p)
    if p < 1:
        value += (1-p) * math.log2(q*(1-p)/(q-1))
    return value


def seed_window_profiles(zf, path, q, L):
    cols = [f'n{a}{b}' for a in range(q) for b in range(q)]
    ncol = q*q
    tables = {w: [[0]*ncol for _ in range(L//2)] for w in WINDOWS}
    sample_ids = set()
    seen_pairs = 0
    with zf.open(path) as raw:
        reader = csv.DictReader(io.TextIOWrapper(raw, encoding='utf-8-sig', newline=''))
        if not set(cols+['sample','r']).issubset(reader.fieldnames or []):
            raise ValueError(f'Unexpected CSV columns: {path}')
        for line in reader:
            k, r = int(line['sample']), int(line['r'])
            if not 0 <= k < 500 or not 1 <= r <= L//2:
                raise ValueError(f'Unexpected sample/r in {path}: {k}/{r}')
            sample_ids.add(k)
            seen_pairs += 1
            vals = [int(line[c]) for c in cols]
            if sum(vals) != 2*L*L:
                raise ValueError(f'Pair count mismatch: {path}: sample {k}, r {r}')
            for name, (start, stop) in WINDOWS.items():
                if start <= k < stop:
                    acc = tables[name][r-1]
                    for j, v in enumerate(vals):
                        acc[j] += v
    if len(sample_ids) != 500 or seen_pairs != 500*(L//2):
        raise ValueError(f'Incomplete sample coverage in {path}: {len(sample_ids)} samples, {seen_pairs} rows')
    result = {}
    for name, table in tables.items():
        original = [sym_mi(counts, q) for counts in table]
        values = []
        for rho in RHOS:
            x = rho*L
            lo = int(math.floor(x)); hi = int(math.ceil(x))
            values.append(original[lo-1] if lo == hi
                          else original[lo-1]*(hi-x)+original[hi-1]*(x-lo))
        result[name] = values
    return result


def score(averages, model, sizes, T):
    cv = []
    for j in range(len(RHOS)):
        xs = [L**EXPONENT[model]*averages[(model, L, T)][j] for L in sizes]
        mu = sum(xs)/len(xs)
        if mu <= 0:
            raise ValueError(f'Invalid CV mean: {model} {T} {j}')
        cv.append(math.sqrt(sum((v-mu)**2 for v in xs)/len(xs))/mu)
    return sum(cv)/len(cv)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--zip', default='physa_main_v11_results.zip')
    ap.add_argument('--out', default='time_window_mi_sensitivity.csv')
    args = ap.parse_args()
    zpath = Path(args.zip)
    if not zpath.is_file():
        raise SystemExit(f'ZIP not found: {zpath.resolve()}')
    profiles = {name: defaultdict(list) for name in WINDOWS}
    seeds = defaultdict(set)
    with zipfile.ZipFile(zpath) as zf:
        jobs = []
        for path in zf.namelist():
            if not path.endswith('/joint_counts.csv'):
                continue
            match = PATTERN.fullmatch(path.split('/')[-2])
            if not match:
                raise ValueError(f'Unknown job naming: {path}')
            model, L, T, seed = match.groups()
            model = 'potts3' if model == 'potts' else model
            jobs.append((path, model, int(L), float(T), int(seed)))
        print(f'Found {len(jobs)} joint-count jobs; expected 1170')
        if len(jobs) != 1170:
            raise ValueError('Expected 1170 complete jobs')
        for i, (path, model, L, T, seed) in enumerate(jobs, 1):
            key = (model, L, T)
            if seed in seeds[key]:
                raise ValueError(f'Duplicate seed {key}/{seed}')
            seeds[key].add(seed)
            this = seed_window_profiles(zf, path, 2 if model == 'ising' else 3, L)
            for w, values in this.items():
                profiles[w][key].append(values)
            if i % 100 == 0 or i == len(jobs):
                print(f'Processed {i}/{len(jobs)} jobs', flush=True)
    if len(seeds) != 117:
        raise ValueError(f'Expected 117 model-size-temp conditions, found {len(seeds)}')
    for key, vals in seeds.items():
        if len(vals) != 10:
            raise ValueError(f'Expected 10 independent seeds for {key}')

    averages = {}
    for w in WINDOWS:
        averages[w] = {key: [sum(v[j] for v in value)/len(value) for j in range(len(RHOS))]
                       for key, value in profiles[w].items()}

    rows = []
    all_pass = True
    for model in ('ising', 'potts3'):
        temps = sorted({T for (m, L, T) in seeds if m == model})
        for subset, sizes in (('development', DEV), ('validation', VAL)):
            ref_T, ref_score = REFERENCE[model, subset]
            for w in WINDOWS:
                curve = [(T, score(averages[w], model, sizes, T)) for T in temps]
                best_T, min_score = min(curve, key=lambda pair: pair[1])
                score_at_ref = next(s for t,s in curve if abs(t-ref_T) < 1e-10)
                margin = sorted(s for _,s in curve)[1] - min_score
                if w == 'full':
                    match = abs(best_T - ref_T) < 1e-10 and abs(min_score-ref_score) < 1e-5
                    all_pass &= match
                    print(f'REFERENCE {model}/{subset}: reproduced={match}; T={best_T:.3f}; score={min_score:.8f}')
                rows.append(dict(model=model, subset=subset, window=w,
                                 sample_start=WINDOWS[w][0], sample_stop=WINDOWS[w][1],
                                 n_samples_per_seed=WINDOWS[w][1]-WINDOWS[w][0],
                                 selected_T=f'{best_T:.3f}', min_score=f'{min_score:.9f}',
                                 score_at_full_selected_T=f'{score_at_ref:.9f}',
                                 gap_to_second_best=f'{margin:.9f}',
                                 unchanged_from_full=(abs(best_T-ref_T)<1e-10)))
    out = Path(args.out)
    with out.open('w', newline='', encoding='utf-8-sig') as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wr.writeheader(); wr.writerows(rows)
    print('\nRESULTS (window / Ising dev, Ising val, Potts dev, Potts val):')
    for w in WINDOWS:
        sel = [r for r in rows if r['window'] == w]
        print(f'{w:10s}: ' + ' | '.join(f"{r['selected_T']} (S={r['min_score']})" for r in sel))
    print(f'\nReference full-data reconstruction: {"PASS" if all_pass else "FAIL -- do not interpret sensitivity"}')
    print('Saved:', out.resolve())
    if not all_pass:
        raise SystemExit('ERROR: Full-data reference mismatch; validate estimator first')
    print('NOTE: post-hoc time-window robustness != proof of equilibrium or independent MI samples.')


if __name__ == '__main__':
    main()
