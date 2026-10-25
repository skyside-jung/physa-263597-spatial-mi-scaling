# PHYSA-263597 reproducibility and figure re-analysis package (draft)

Manuscript: *Inverse Identification of Critical Transitions from Spatial Information Scaling: Benchmark Tests in Ising and Three-State Potts Models*.

## What is present

- `source/analysis_only.zip`: protocol, saved seed-level profiles and diagnostics; NOT the 1,170 raw job directories.
- `source/mmi.zip`: original Python simulation/inference package supplied by the author.
- `recreate_figures.py`: self-contained secondary script reconstructing the numerical curves in Figures 2--5 from `merged/profiles.csv`. Does not call a previously prepared score CSV for calculations.
- `figures_recreated/`: **independently rendered verification figures**, not publisher-ready original plots.
- `checks/recomputed_collapse_scores.csv`: 52 independently recomputed temperature-by-observable scores compared with saved scores.
- `checks/figure_reproduction_results.json`: audit results, including 90/90, 0/90, 0/90 Figure 5 test.
- `checks/verify_archive.py`, `checks/recompute_table1.py`: audits requiring the separately preserved `physa_main_v11_results.zip`.
- `supplementary/`: additional post hoc sensitivity results and one source program, as supplied by the author. These sensitivity analyses have **not all been independently rerun** here.

## Re-create figures without the raw 484MB archive

```bash
python -m pip install -r requirements.txt
python recreate_figures.py
```

This produces Figure 2--5 verification figures and a JSON test report. Requires only `source/analysis_only.zip` and numpy/matplotlib. It should reproduce all 52 saved collapse scores to floating-point tolerance and the Figure 5 sign counts [90,0,0].

## Verify original raw-job archive (separate file)

```bash
python checks/verify_archive.py /path/to/physa_main_v11_results.zip
python checks/recompute_table1.py /path/to/physa_main_v11_results.zip
```

The earlier audit verified 3,531 archived SHA256-tracked files and 1,170 joint-count records and regenerated all eight Table 1 numerical results. The full original archive is **not included** in this light-weight package.

## Important limits

- This package does **not** reproduce all 1,170 Monte Carlo jobs from scratch.
- It does **not** verify the original Figure 2--5 plotting script, which was not among the supplied files; instead it reconstructs the numerical panels from preserved seed-level profiles.
- It does **not** recompute every supplementary bootstrap, temporal-window, correction, and nested-size analysis; supplied output files are preserved for audit.
- A public repository/DOI has **not yet been created**. Do not state that the data have been publicly deposited until publication is confirmed.
- Avoid exposing locally configured credentials or unrelated files if depositing these materials.
