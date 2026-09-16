# Polygenic Risk Scores for Hypertensive Disorders of Pregnancy and Later-Life Cardiovascular Disease in Women

This repository contains code supporting the manuscript:

**Polygenic Risk Scores for Hypertensive Disorders of Pregnancy and Later-Life Cardiovascular Disease in Women**

## Study overview

This study investigated whether genetic predisposition to hypertensive disorders of pregnancy (HDP), measured using polygenic risk scores (PRS), was associated with later-life cardiovascular disease (CVD) outcomes in women from the Busselton Health Study (BHS).

The analytical cohort included 1,284 women from the 1994-95 Busselton Health Study with genotype data, pregnancy history, cardiovascular risk factors, and linked hospitalisation and death records.

Four HDP-related PRS were evaluated:

* PGS003586 - preeclampsia
* PGS003587 - gestational hypertension
* PGS004593 - preeclampsia-related PRS
* Broad HDP PRS derived from FinnGen-based HDP GWAS summary statistics

## Repository contents

This repository contains code for:

1. Calculation of HDP-related PRS using `pgsc_calc`.
2. Fine-Gray competing-risk analyses of first CVD outcomes.
3. Generation of manuscript result tables and forest plots.

Repository structure:

```text
HDP-CVD-PRS/
├── README.md
├── .gitignore
├── CITATION.cff
├── code/
│   ├── 01_calculate_PRS.sh
│   └── 02_FineGray_PRS_analysis.R
├── data/
│   └── README.md
└── results/
    └── README.md
```

## Polygenic risk score calculation

PRS were calculated using `pgsc_calc` v2.2.0.

Three scores were obtained from the Polygenic Score Catalog:

* PGS003586
* PGS003587
* PGS004593

A separate broad HDP PRS was calculated using a FinnGen-based HDP score.

PRS calculation used:

* `pgsc_calc` v2.2.0
* Nextflow v25.04.6
* Singularity v4.1.0
* GRCh38 as the target genome build
* HGDP + 1000 Genomes reference data for ancestry normalisation

The PRS calculation commands used in this study are provided in:

`code/01_calculate_PRS.sh`

Input sample sheets, genotype data, ancestry reference files and other restricted or locally stored resources are not distributed with this repository. Example paths in the script should be replaced with paths appropriate to the user's computing environment.

### Broad HDP PRS

The broad HDP PRS was developed using summary statistics from the FinnGen Consortium Data Freeze R8v4 HDP GWAS and the PRS-CS Bayesian polygenic prediction method.

The HDP-PRS model and associated resources are publicly available from:

* GitHub: https://github.com/dokyoonkimlab/hdp-prs-finngen-r8
* Zenodo: https://doi.org/10.5281/zenodo.14211027

The broad HDP score required genome-build liftover during the `pgsc_calc` workflow before scoring against the GRCh38 target genotype data.

Associated publication:

Jung SH, Kim H, Jung YM, Shivakumar M, Xiao B, Kim J, Jang B, Yun JS, Won HH, Park CW, Park JS, Jun JK, Kim D, Lee SM. Healthy lifestyle reduces cardiovascular risk in women with genetic predisposition to hypertensive disorders of pregnancy. *Nature Communications*. 2025;16:1463. doi:10.1038/s41467-025-56107-2.

## Fine-Gray competing-risk analysis

PRS were standardised to mean 0 and standard deviation 1 in the analytical cohort before analysis.

Fine-Gray competing-risk models were used to evaluate associations between each PRS and specific first CVD outcomes, with group-specific PRS effects estimated for women with and without a history of HDP.

For each specific first CVD outcome:

* the outcome of interest was treated as the event of interest;
* other first CVD outcomes were treated as competing events; and
* women without an incident CVD event were censored.

Models included the PRS, HDP history and a PRS-by-HDP-history term to obtain PRS effect estimates separately for women with and without a history of HDP. 

Models were adjusted for:

* age
* age squared
* systolic blood pressure
* waist circumference
* total cholesterol to HDL-cholesterol ratio
* diabetes
* lipid-lowering medication use
* ever-smoking status
* genetic principal components 1-10

Group-specific effect estimates are reported as subdistribution hazard ratios (sHR) with 95% confidence intervals and p-values per 1-standard deviation increase in PRS.

The Fine-Gray analysis and manuscript figure-generation code is provided in:

`code/02_FineGray_PRS_analysis.R`

## Analysis outputs

The Fine-Gray analysis script generates non-identifiable aggregate outputs including:

* group-specific sHR estimates
* 95% confidence intervals
* p-values
* aggregate event counts
* manuscript-ready result tables
* individual PRS forest plots
* combined manuscript Figures 2 and 3
* R session information

Participant-level analysis datasets and individual PRS values are not written to the public repository.

## Data availability

Individual-level Busselton Health Study genotype, phenotype, pregnancy-history and linked health data are not included in this repository because they are subject to data-access, governance and privacy restrictions.

The analysis scripts require authorised access to the corresponding BHS datasets.

No participant-level BHS data, participant identifiers, genotype files, individual PRS values, linked hospitalisation records, mortality records, pregnancy-history records or participant mapping files are distributed through this repository.

Researchers interested in accessing BHS data should follow the relevant BHS data-access and governance procedures.

## Software

The principal software used included:

* `pgsc_calc` v2.2.0
* Nextflow v25.04.6
* Singularity v4.1.0
* R v4.5.2
* `cmprsk`
* `dplyr`
* `tidyr`
* `readr`
* `ggplot2`
* `patchwork`
* `stringr`

## Reproducibility

This repository provides the code used for PRS calculation and the Fine-Gray competing-risk analyses reported in the manuscript.

The restricted individual-level BHS datasets required to reproduce the complete analyses are not publicly distributed. Users with authorised access to the relevant data can adapt the input paths in the scripts to their local computing environment.

## Citation

If using this code, please cite the associated manuscript:

**Krishnamurthy G, Temple JA, Ahmed MA, Hui J, Meikle PJ, Giles C, Phan HT, Brennecke SP, Moses EK, Melton PE. Polygenic Risk Scores for Hypertensive Disorders of Pregnancy and Later-Life Cardiovascular Disease in Women.**

Journal and DOI details will be added following publication.

For the broad HDP PRS, please cite:

Jung SH, Kim H, Jung YM, Shivakumar M, Xiao B, Kim J, Jang B, Yun JS, Won HH, Park CW, Park JS, Jun JK, Kim D, Lee SM. Healthy lifestyle reduces cardiovascular risk in women with genetic predisposition to hypertensive disorders of pregnancy. *Nature Communications*. 2025;16:1463. doi:10.1038/s41467-025-56107-2.

The corresponding HDP-PRS resource is available at:

https://doi.org/10.5281/zenodo.14211027

Users should also cite the original publications and Polygenic Score Catalog resources for PGS003586, PGS003587 and PGS004593 where appropriate.
