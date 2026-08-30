# COVID-19 Hospital Burden Modelling: Age, Deprivation, and Region in England (2020–21)

## Overview

This repository contains the R code and outputs for an MSc dissertation project at Imperial College London (MSc Health Data Analytics and Machine Learning). The project constructs a deterministic SEIRD compartmental model stratified by age, IMD deprivation decile, and ITL1 region to quantify the joint contribution of these factors to COVID-19 hospital burden in England during the 2020–21 autumn–winter wave.

## Key Findings

- Regional differences in COVID-19 hospitalisation burden were almost entirely explained by age × deprivation population composition rather than differential transmission dynamics
- Approximately 36% of modelled hospital admissions and in-hospital deaths were attributable to deprivation-related disparities in clinical severity (counterfactual analysis)
- School closure reduced R₀ below the epidemic threshold in less deprived deciles (5–10) but failed to do so in the four most deprived deciles, where non-school contacts alone sustained transmission

## Repository Structure
Dissertation/
├── Rscripts/ # Analysis scripts (numbered in order)
├── data/ # Input data (geography, population, contacts)
├── output/ # Model outputs and figures
└── docs/ # Documentation including pipeline log

## Scripts

| Script | Description |
|--------|-------------|
| `01_build_geography_lookup.R` | Build LSOA → LAD → ITL1 geography lookup |
| `02_add_imd_clean.R` | Link IMD 2019 deprivation scores to LADs |
| `03_aggregate_lad_imd.R` | Aggregate LAD-level deprivation deciles |
| `04_build_regression_data.R` | Prepare regression input data |
| `05_build_population_strata.R` | Build age × IMD × ITL1 population matrix |
| `06_seird_unstratified.R` | Unstratified SEIRD model prototype |
| `07_prepare_model_inputs.R` | Prepare contact matrices and model inputs |
| `08_refine_gamma_hd_hr.R` | Derive hospital discharge rates from Knock et al. |
| `09_age_imd_stratified_odin.R` | Age- and IMD-stratified SEIRD model (odin2/dust2) |
| `10_region_stratified_model.R` | Regional stratification |
| `11_verify_R0.R` | Verify basic reproduction number computation |
| `12_regression_stringency.R` | NPI stringency regression |
| `13_monty_fitting.R` | Bayesian MCMC fitting via monty package |
| `14_age_imd_gradient_plots.R` | Age and deprivation gradient visualisations |
| `15_region_age_imd_plots.R` | Regional × age × IMD burden analysis |
| `16_counterfactual.R` | Counterfactual analysis (CF1 and CF2 scenarios) |
| `17_school_closure.R` | School closure scenario and R₀ analysis |

## Data Sources

- **NHS England COVID-19 Hospital Admissions**: Weekly trust-level admissions via `epiforecasts/covid19.nhs.data` R package (March 2022 release)
- **IMD 2019**: English Indices of Deprivation, Ministry of Housing, Communities and Local Government
- **ONS Geography**: LSOA → LAD → ITL1 lookup (April 2021 boundaries)
- **ONS Population**: Mid-2020 population estimates
- **POLYMOD**: Social contact matrices via `socialmixr` R package (Mossong et al. 2008)

## Model

The model is implemented using:
- [`odin2`](https://github.com/mrc-ide/odin2) — ODE model specification
- [`dust2`](https://github.com/mrc-ide/dust2) — Stochastic simulation
- [`monty`](https://github.com/mrc-ide/monty) — Bayesian MCMC fitting

## Key Parameters

- **β**: Transmission probability per contact — fitted via MCMC, fixed at decile 1 posterior median (0.031) for all analyses
- **π_a**: Age- and IMD-specific clinical fraction — from Goodfellow et al. (2024) via LOESS regression of 2021 Census health prevalence
- **h_a**: Age-specific hospitalisation probability — derived from Knock et al. (2021) IHR estimates

## Supervisors

Prof Marc Baguelin and Dr Christian Morgenstern, Imperial College London

## Reference

Huang S. *How Age, Deprivation, and Region Shaped COVID-19 Hospital Burden in England: A Compartmental Modelling Study, 2020–2021*. MSc Dissertation, Imperial College London, 2025.
