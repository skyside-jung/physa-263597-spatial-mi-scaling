# PHYSA-263597 reproducibility and figure re-analysis package

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
- The lightweight code and analysis repository is public: https://github.com/skyside-jung/physa-263597-spatial-mi-scaling. **No Zenodo DOI has been issued and the full raw dataset has not yet been publicly deposited.** Do not claim publication of the full raw dataset until that occurs.
- Avoid exposing locally configured credentials or unrelated files if depositing these materials.

## How to cite and reuse

- Software citation metadata: [`CITATION.cff`](CITATION.cff). A Zenodo DOI will be added once a software release is archived.
- Original Python source code and verification scripts: [MIT License](LICENSE).
- Research-generated data tables, protocol metadata, and verification figure outputs: [CC BY 4.0](DATA_LICENSE.md). The license does not cover third-party content or materials for which the depositor lacks rights.
- **Full raw archive:** retained separately; its public dataset DOI and download URL will be added after an actual deposit.

This repository provides an auditable analysis-stage reproduction of the principal results; it does not claim complete re-execution of the original 1,170 Monte Carlo production jobs.
