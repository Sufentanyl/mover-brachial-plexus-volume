# Local anesthetic volume and early rescue analgesia after single-shot brachial plexus block

Analysis code for the retrospective cohort study:

> Local anesthetic volume and early rescue analgesia after single-shot brachial plexus block: a retrospective cohort study. *Submitted to Minerva Anestesiologica.*

The study used the Medical Informatics Operating Room Vitals and Events Repository (MOVER), version 2 (Epic component, University of California, Irvine, 2017–2023), to examine the association between cumulative perineural local anesthetic volume and time to in-hospital rescue analgesia in 571 adults who received a single-shot brachial plexus block.

## Data

**No patient data are included in this repository.** MOVER is available from the UCI Machine Learning Repository (doi:10.24432/C5VS5G) to researchers who sign the MOVER data use agreement, which does not permit redistribution of the data. Please cite the database paper when using MOVER:

> Samad M, Angel M, Rinehart J, Kanomata Y, Baldi P, Cannesson M. Medical Informatics Operating Room Vitals and Events Repository (MOVER): a public-access operating room database. JAMIA Open 2023;6:ooad084.

After obtaining access, set the path to the downloaded MOVER files as described under *Running the analysis*.

## Repository structure

```
analysis/   cohort construction, exposure and outcome definitions, primary,
            sensitivity and supplementary analyses
figures/    Figure 1 (flowchart) and Figures 2, 3 and Supplementary Figure 1
            (cumulative incidence) of the manuscript
```

## Environment

- Python 3.12.3 with lifelines 0.30.3 (see `requirements.txt` for all package versions)
- R 4.3.3 with coxphf 1.13.4 (Firth-penalized Cox regression)

## Running the analysis

<!-- To be completed with the script order and the data path setting once the analysis scripts are added. -->

## Figures

`figures/plot_cif_figures.py` produces Figures 2, 3 and Supplementary Figure 1 from two files exported by the analysis scripts (`primary.csv`: `time_h`, `status`, `vol_group`; `lower.csv`: `time_h`, `status`, `block`). These files contain individual-level data and are not distributed.

```
python figures/plot_cif_figures.py primary.csv lower.csv out_dir
python figures/fig1_flowchart.py out_dir/Figure1.png out_dir/Figure1.tif
```

## Correspondence between scripts and the manuscript

<!-- To be completed: Table I–III, Figures 1–3, Supplementary Tables I–V, Supplementary Figure 1. -->

## License

MIT License (see `LICENSE`).
