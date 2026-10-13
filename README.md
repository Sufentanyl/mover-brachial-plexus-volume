# Local anesthetic volume and early rescue analgesia after single-shot brachial plexus block

Analysis code for the retrospective cohort study:

> Feng Q, Jing Y, Han C. Local anesthetic volume and early rescue analgesia after single-shot brachial plexus block: a retrospective cohort study. *Prepared for submission to Minerva Anestesiologica.*

The study used the Medical Informatics Operating Room Vitals and Events Repository (MOVER), version 2 (Epic component, University of California, Irvine, 2017–2023), to examine the association between cumulative perineural local anesthetic volume and time to in-hospital rescue analgesia in 571 adults who received a single-shot brachial plexus block.

## Data

**No patient data are included in this repository.** MOVER is available from the UCI Machine Learning Repository (doi:10.24432/C5VS5G) to researchers who sign the MOVER data use agreement, which does not permit redistribution of the data. Please cite the database paper when using MOVER:

> Samad M, Angel M, Rinehart J, Kanomata Y, Baldi P, Cannesson M. Medical Informatics Operating Room Vitals and Events Repository (MOVER): a public-access operating room database. JAMIA Open 2023;6:ooad084.

## Repository structure

```
data_prep/   00–07  upstream pipeline from the MOVER Epic extracts to the cohort,
                    exposure, pain scores, primary outcome and dexamethasone files
analysis/           analysis dataset (n = 571), primary model, sensitivity and
                    supplementary analyses; run_all.sh runs the listed steps
figures/            Figure 1 (flowchart) and Figures 2, 3 and Supplementary Figure 1
config.yaml         description of the path variables
```

Code comments are partly in Chinese. `data_prep/README.md` documents the inputs, outputs and key definitions of each upstream step.

## Environment

- Python 3.12.3; packages and versions in `requirements.txt` (lifelines 0.30.3, pandas 2.3.3, numpy 2.5.3, scikit-learn 1.9.1, statsmodels 0.15.0)
- R 4.3.3 with `survival`, `cmprsk` and `coxphf` (1.13.4). Multiple imputation is implemented in base R (fully conditional specification) in `analysis/01_mice_extended.R`.
- Random seed 20261010 for bootstrap and Monte Carlo procedures.

## Running the available analysis steps

`analysis/run_all.sh` requires the locally prepared input files described in
`data_prep/README.md`, including `cohort_info.json`, `pi_pain_meas.txt`, and the two
Module B extracts. It stops when required inputs are missing. It does not run the
subgroup and interaction scripts or generate the figures; those scripts are
listed separately below.

All paths are set by environment variables (see `config.yaml`):

```bash
export MOVER_DIR=/path/to/MOVER   # MOVER v2 Epic extracts
export WORK_DIR=/path/to/work     # intermediate files
export OUT_DIR=/path/to/output    # analysis dataset and results
bash analysis/run_all.sh
```

Steps:

1. `data_prep/00–07`: CPT screening (64415/64416), medication extraction, cohort construction, case-level (LOG_ID) linkage, cumulative volume within 60 min, pain scores, primary outcome, dexamethasone.
2. `analysis/build_571_dataset.py`: analysis dataset; excludes blocks given after the end of anesthesia (n = 573) and blocks given more than 24 h before the start of anesthesia (n = 571). `verify_trigger_condition.py` checks that requiring the triggering pain score to be recorded after the end of anesthesia leaves the events unchanged.
3. `analysis/run_v26.R`: primary cause-specific Cox model, Schoenfeld test, time-split Firth model, sensitivity analyses (Fine–Gray, Firth, overlap weighting and others).
4. Sensitivity analyses: `01_mice_extended.R` (extended model with multiple imputation), `02_no_shoulder.R`, `03_dedup.py` + `03_dedup_cox.R` (same-minute duplicate records), `interaction_v26.R` (volume × surgical category), `v26_subg2.R` / `v26_subg3.R` (subgroups), `sens_t0_early_excl.R` (check of the 24-h rule), `bench_v26.R` (benchmark).
5. Supplementary analyses: `c1_abc.py` (key secondary outcome, any opioid use: models A, B and C), `c2_riskdiff.py` (cumulative incidence differences), `d1_ipw_291.py` (inverse probability weighting for missing pain scores), `e3_bias_v2.py` (probabilistic bias analysis), `evalue.py` (E-value), `04_moduleB_revised.py` (conversion to general anesthesia and intraoperative opioid use), `supp_table1.py` (analgesics defining the rescue events, Supplementary Table I), `lower_extremity.py` (lower-extremity cohort, CPT 64447/64445, Supplementary Figure 1), `export_cif_inputs.py` (input files for the cumulative incidence figures).

Upstream extraction scripts: `data_prep/extract_pain_meas.py` (pain measurements from the Epic flowsheet table, input to `data_prep/05`), `data_prep/extract_proc_meds.py` (procedure events and medication records for Module B).

## Correspondence with the manuscript

| Manuscript | Script |
|---|---|
| Table I | Descriptive table-generation script not included; `run_v26.R` does not generate Table I |
| Table II | `analysis/run_v26.R` |
| Table III | `analysis/run_v26.R`, `01_mice_extended.R`, `02_no_shoulder.R`, `03_dedup*.{py,R}` |
| Interaction and subgroups | `analysis/interaction_v26.R`, `v26_subg2.R`, `v26_subg3.R` |
| Supplementary Table IV | `analysis/c1_abc.py`, `c2_riskdiff.py`, `04_moduleB_revised.py` |
| Supplementary Table V | `analysis/d1_ipw_291.py`, `e3_bias_v2.py` |
| Figures 1–3, Supplementary Figure 1 | `figures/` |

## Figures

```bash
python figures/plot_cif_figures.py primary.csv lower.csv out_dir
python figures/fig1_flowchart.py out_dir/Figure1.png out_dir/Figure1.tif
```

The input files contain individual-level data and are not distributed.

## License

MIT License (see `LICENSE`).
