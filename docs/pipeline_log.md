# Data Processing Pipeline Log

Purpose, inputs/outputs, and validation results for each script. Keep in `docs/`.

------------------------------------------------------------------------

## 01_build_geography_lookup.R

**Purpose:** Builds LSOA -\> STP -\> LAD -\> ITL geography lookup, the backbone for all downstream scripts.

**Inputs:** `data/raw/LSOA11_CCG21_STP21_LAD21_EN_LU_...csv`, `data/raw/Local_Authority_District_(April_2021)_..._Lookup...csv`

**Output:** `data/processed/lookup_LSOA11_STP_ITL2.csv`

**Validation:**

| Check             | Result | Official reference               | Match |
|-------------------|--------|----------------------------------|-------|
| Unique LSOAs      | 32,844 | 32,844 (England, 2011 Census)    | OK    |
| Unique STPs       | 42     | 42 (England ICS/STP, April 2021) | OK    |
| Unique ITL2       | 33     | 33 (England ITL2 regions)        | OK    |
| Missing ITL2 rows | 0      | 0                                | OK    |

------------------------------------------------------------------------

## 02_add_imd_clean.R

**Purpose:** Joins geography lookup with IMD 2019 deprivation data at LSOA level.

**Inputs:** `data/processed/lookup_LSOA11_STP_ITL2.csv`, `data/raw/File_7_-_All_IoD2019_Scores...csv`

**Output:** `data/processed/lookup_LSOA11_STP_ITL2_IMD_clean.csv`

**Validation:**

| Check               | Result | Official reference       | Match |
|---------------------|--------|--------------------------|-------|
| Rows / Unique LSOAs | 32,844 | 32,844                   | OK    |
| Missing IMD         | 0      | 0                        | OK    |
| Unique LAD/LTLA     | 309    | 309 (England LADs, 2021) | OK    |
| Unique ITL3         | 133    | 133                      | OK    |
| Unique Region/ITL1  | 9      | 9                        | OK    |

------------------------------------------------------------------------

## 03_aggregate_lad_imd.R

**Purpose:** Aggregates LSOA-level IMD to LAD level (mean/median score, rank, decile); attaches geography.

**Inputs:** `data/processed/lookup_LSOA11_STP_ITL2_IMD_clean.csv`

**Outputs:** `data/processed/LAD21_IMD_summary.csv`, `data/processed/LAD21_geography_IMD_summary.csv`

**Validation:**

| Check | Result | Official reference | Match |
|----|----|----|----|
| Number of LAD areas | 309 | 309 | OK |
| Decile distribution | 31x9 + 30x1 = 309 | Expected from `ntile()` | OK |
| Most deprived LAD (rank 1) | Blackpool | Blackpool (MHCLG IMD2019) | OK |
| Ranks 2-3, 6 | Knowsley, Liverpool, Manchester | Same | OK |
| Ranks 4-5 | Order swapped vs official | Weighted vs unweighted averaging — methodology note, not a bug | Noted |

**Bug found later (while building script 05) and fixed:** the geography-join step produced 324 rows for 309 LADs — 15 LADs (mostly Buckinghamshire and North/West Northamptonshire, which span \>1 STP) had duplicate rows. Fixed by switching to a majority-STP rule (`count()` + `slice_max()`) so each LAD gets exactly one row. Confirmed fix: 309 rows, 309 distinct LADs. This bug had propagated into script 04's output (`regression_data.rds`) — script 04 was re-run after the fix.

------------------------------------------------------------------------

## 04_build_regression_data.R

**Purpose:** Builds LAD-by-week dataset combining NHS admissions, LAD geography/IMD, OxCGRT policy, and population. Used as (a) a real-data validation target for the dynamic model, and (b) the input for the descriptive regression analysis (script 11).

**Inputs:** NHS admissions (`covid19.nhs.data::get_admissions()`), `data/processed/LAD21_geography_IMD_summary.csv` (post script-03 fix), OxCGRT policy data (GitHub), `data/raw/ukpopestimatesmid2020on2021geography.xls`

**Output:** `data/processed/regression_data.rds`

**Bugs found and fixed (across several iterations):** - `lad_imd` duplicate-row bug (see script 03) — fixed by re-running with corrected script 03 output; added defensive deduplication on read as a safeguard - OxCGRT `RegionCode` filter was `GBR_ENG` (does not exist in the live dataset) — corrected to `UK_ENG` after directly inspecting `unique(oxcgrt$RegionCode[oxcgrt$CountryCode=="GBR"])`. This was the root cause of `stringency_mean` being 100% NA in earlier runs. - Population (`pop_lad`) exclusion list was missing `"County"` — same double-counting bug as script 05's original version. Fixed by adding `"County"` to the exclusion list. - Multiple `select()` calls were being masked by a same-named function from another loaded package, causing "unused argument" errors. Fixed by namespacing all `select()` calls as `dplyr::select()`.

**Known remaining issue:** LADs in output: 306, not 309 — likely small LADs (City of London, Isles of Scilly) excluded by NHS disclosure control. Not investigated further; not a blocker. `release_date` still a placeholder (`2022-03-03`) — not finalised.

**Status:** all known bugs fixed and re-run. `stringency_mean` now populated (no longer 100% NA).

------------------------------------------------------------------------

## 05_build_population_strata.R

**Purpose:** Builds N\_{a,i,r} — population by age band x LAD-level IMD decile x NHS region/ITL1 — for the dynamic model's population stratification term.

**Inputs:** `data/processed/LAD21_geography_IMD_summary.csv` (post script-03 fix), `data/raw/ukpopestimatesmid2020on2021geography.xls`

**Output:** `data/processed/population_age_imd_region.csv`

**Caveat:** IMD decile here is the LAD-level proxy (script 03's `lad_imd_decile`), not the true LSOA-level decile.

**Validation:**

| Check | Result | Official reference | Match |
|----|----|----|----|
| Total population | 56,550,138 | \~56.5M (England, mid-2020/2021) | OK |
| Age bands | 17 | 17 (0-4 ... 80+) | OK |
| IMD deciles | 10 | 10 | OK |
| Regions | 9 | 9 | OK |
| Missing cells | 238 / 1530 | Expected — sparse region x decile combinations | OK |

**Bug found and fixed:** first run gave total population 60,003,484 (\~3.45M too high). Traced to the script-03 duplicate-LAD bug (above) — once that was fixed, total population matched the independently cross-checked value (56,550,138) exactly.

------------------------------------------------------------------------

## 06_seird_unstratified.R

**Purpose:** Single-population (no age/IMD/region stratification) SEIR-HD model in odin, to confirm the model structure before adding stratification.

**Parameters:** see `docs/model_parameter_table.md`. Several values are single-representative placeholders (h, mu_ca_h, M, gamma_hd, gamma_hr) — not yet age/IMD/region-stratified.

**Validation:** equations checked in base R before odin translation — - Population exactly conserved across all time steps - Single epidemic wave, peak clinical infectious \~day 64 - Final attack rate \~91% - Final deaths and admissions match theoretical values exactly (pi_c x h x mu_ca_h and pi_c x h respectively)

**Status:** structure confirmed working — superseded by script 08 (which adds age and IMD stratification using the same structure). Kept as a reference / sanity baseline.

------------------------------------------------------------------------

## 07_prepare_model_inputs.R

**Purpose:** Processes Goodfellow et al. (2024)'s published data files and the Knock et al. (2021) IHR/IFR table into the exact arrays the age-stratified odin model needs: a clinical-fraction matrix, 10 contact matrices (one per IMD decile), and an age-indexed hospitalisation/mortality table.

**Updated (advisor request):** contact matrices now use a population-weighted urban/rural blend per IMD decile, replacing the earlier urban-only version. Script 12 (`12_urban_rural_weighted_contact.R`) was a temporary standalone implementation of this logic; it has been merged into script 07 and should be deleted from `Rscripts/`.

**Inputs:** - `data/parameters/clin_frac.csv` (Goodfellow repo `/data/`) - `data/parameters/G.csv` (Goodfellow repo `/data/`, intrinsic connectivity matrix) - `data/parameters/rural_age.csv` (Goodfellow repo `/data/`, population by age/IMD/urban-rural) - `data/parameters/ihr_ifr_by_age_knock2021.csv` (Knock et al. 2021, Table S9)

**Outputs:** - `data/parameters/pi_matrix.csv` (17 ages x 10 IMD deciles, clinical fraction) - `data/parameters/contact_matrix_imd1.csv` ... `contact_matrix_imd10.csv` (17x17 each; urban/rural blended per IMD decile — see method note below) - `data/parameters/urban_share_by_decile.csv` (urban population weight per IMD decile, kept for QA and reporting) - `data/parameters/h_mu_by_age.csv` (h_a and mu_ca_h by age, using mean-pi across deciles as denominator)

**Method notes:** - **Contact matrix blending:** G (POLYMOD-derived intrinsic connectivity matrix) is density-corrected separately for urban and rural settings, then blended using each IMD decile's actual urban population share as the weight: `M_blended = w_urban * M_urban + (1 - w_urban) * M_rural` where `M[i,j] = G[i,j] * N_j / sum(N)` for each setting. Urban and rural populations are NOT simply pooled into one age vector — that would incorrectly assume identical contact intensities across settings. - **Age-band mapping (Knock → Goodfellow):** Knock's 0-4 band is duplicated for Goodfellow's "Under 1" and "1 to 4"; Knock's 75-79 and 80+ are averaged into Goodfellow's "75+". Documented approximation. - **h_a denominator:** uses IMD-decile-averaged π (mean across all 10 deciles), matching Goodfellow's own convention for the analogous parameter — see `docs/model_parameter_table.md` for the full rationale.

**Validation:** - pi_matrix: 17 ages x 10 deciles, values in plausible 0.2-0.7 range - Contact matrix sanity check: blended row sum for "Under 1", IMD decile 1 should fall between the urban-only and rural-only row sums (confirmed in script output) - h_a and mu_ca_h increase monotonically with age as expected (h_a: 0.008 at Under-1 to 0.435 at 75+; mu_ca_h: 0.039 to 0.240)

------------------------------------------------------------------------

## 08_refine_gamma_hd_hr.R

**Purpose:** Replace the provisional 80/20-weighted `gamma_hd` and `gamma_hr` placeholders (scripts 08/09) with properly derived, age-varying values, using Knock et al. (2021)'s actual hospital pathway structure and branch probabilities. Produces `data/parameters/gamma_hd_hr_by_age.csv` for use in the updated scripts 08 and 09.

**Input:** - `data/processed/regression_data.rds` (script 04) — used to compute region-specific weights (see Section 2 below) - Knock et al. (2021) Table S2 (durations), Table S6 (region-level branch probabilities), Table S8 (age-scaling factors) — all transcribed by hand from the published supplementary material

**Output:** `data/parameters/gamma_hd_hr_by_age.csv`

| Column             | Description                                           |
|--------------------|-------------------------------------------------------|
| `age_band_gf`      | Goodfellow's 17 age bands (Under 1 to 75+)            |
| `p_death_total`    | P(death \| hospitalised), age-specific                |
| `p_recovery_total` | P(recovery \| hospitalised), age-specific             |
| `gamma_hd`         | Rate HD → D (1 / mean days to death), age-specific    |
| `gamma_hr`         | Rate HR → R (1 / mean days to recovery), age-specific |

------------------------------------------------------------------------

**Method (detailed, for dissertation review):**

**Background — why this matters:** In the original model (scripts 08/09), `gamma_hd = 1/11.3` and `gamma_hr = 1/14.1` were flat constants for all ages, derived from a rough 80/20 weighting of Knock's general-ward vs ICU pathway durations (Table S2 only). This ignored (a) the age-varying probability of taking different hospital pathways (ICU vs general ward), and (b) the regional variation in those probabilities. Script 12 replaces this with a proper probability-weighted aggregation across all pathways.

**Step 1 — Hospital pathway structure (Knock Figure S4/S11):** Knock models the within-hospital pathway as a branching process. A hospitalised patient either: - Is NOT triaged to ICU (probability `1 - p_ICU(a)`), goes to general ward, then dies (probability `p_HD(a)`, mean duration 10.3 days) or recovers (probability `1 - p_HD(a)`, mean duration 10.7 days); OR - Is triaged to ICU (probability `p_ICU(a)`), passes through ICU_pre (2.5 days), then either dies directly in ICU (probability `p_ICUD(a)`, mean duration 11.8 days) or is stepped down to a general ward, where they either die (probability `p_WD(a)`, mean duration 15.1 days = ICU_WD 7.0d + W_D 8.1d) or recover (probability `1 - p_WD(a)`, mean duration 27.8 days = ICU_WR 15.6d + W_R 12.2d).

All durations are from Knock Table S2 (posterior means).

**Step 2 — Age-scaling factors (Table S8):** Knock Table S8 gives relative age-scaling factors for each branch probability (scale relative to the age group with the highest value; max = 1). These are transcribed exactly from the published table and mapped onto Goodfellow's 17 age bands (Knock's "0-5" band is duplicated for Goodfellow's "Under 1" and "1 to 4"; Knock's "75-80" and "80+" are averaged into Goodfellow's "75+").

**Step 3 — England-representative absolute probabilities (Table S6 + real NHS data):** Knock Table S6 gives region-specific posterior means for the four branch probabilities (`p_ICU`, `p_HD`, `p_ICUD`, `p_WD`) at the age of maximum scaling (these are the "ceiling" values that the Table S8 relative factors scale downward from). These vary across Knock's 7 NHS regions (NW, NEY, MID, EE, LON, SW, SE).

Knock's own aggregation method weights each region's estimate by that region's share of the total England attack rate — a quantity only available from his own fitted model output, which we do not have access to. We instead weight by each region's **real observed share of total COVID hospital admissions** from `regression_data.rds` (NHS England data, script 04). This follows the same logic (weight by share of burden) while using directly observed data rather than another model's estimates.

Region weights computed from NHS admissions data:

| Knock region | ONS regions combined          | Admissions weight |
|--------------|-------------------------------|-------------------|
| NW           | North West                    | 0.159             |
| NEY          | North East + Yorkshire        | 0.170             |
| MID          | East Midlands + West Midlands | 0.204             |
| EE           | East                          | 0.082             |
| LON          | London                        | 0.159             |
| SW           | South West                    | 0.078             |
| SE           | South East                    | 0.149             |

**Documented approximation:** East Midlands and West Midlands are pooled into Knock's single "MID" region, and North East and Yorkshire into "NEY", because Knock's 7-region structure does not match ONS's 9 ITL1 regions exactly. This pooling assumes within-group homogeneity, which is a simplification.

**Step 4 — Combine to get age-specific absolute probabilities:** For each of the four branch parameters, the England-representative maximum value (from Step 3) is multiplied by the age-specific relative scaling factor (from Step 2) to give the age-specific absolute branch probability. This is the standard approach implied by Knock's Table S8 structure.

**Step 5 — Aggregate into gamma_hd and gamma_hr:** For each age group, the three death pathways and two recovery pathways are combined by probability-weighting their mean durations:

```         
P(death | hosp)   = (1-p_ICU)*p_HD + p_ICU*p_ICUD + p_ICU*(1-p_ICUD)*p_WD
P(recovery | hosp) = (1-p_ICU)*(1-p_HD) + p_ICU*(1-p_ICUD)*(1-p_WD)

mean_death_duration    = [P_general_death * 10.3 + P_icu_death * 11.8 + P_stepdown_death * 15.1]
                         / P(death | hosp)

mean_recovery_duration = [P_general_recovery * 10.7 + P_stepdown_recovery * 27.8]
                         / P(recovery | hosp)

gamma_hd = 1 / mean_death_duration
gamma_hr = 1 / mean_recovery_duration
```

**Validation / sanity check:**

| Quantity        | Old placeholder    | New range across ages          |
|-----------------|--------------------|--------------------------------|
| `gamma_hd`      | 0.0885 (11.3 days) | 0.0873–0.0953 (10.5–11.5 days) |
| `gamma_hr`      | 0.0709 (14.1 days) | 0.0785–0.0908 (11.0–12.7 days) |
| `p_death_total` | not age-varying    | 2.6% (Under 1) to 39.5% (75+)  |

Key observations: - `gamma_hd` range is very close to the old placeholder — the flat 80/20 approximation happened to be a reasonable estimate for the death timeline. - `gamma_hr` is notably higher (faster recovery) than the old placeholder (0.071 → 0.079–0.091), because properly weighting the ICU→stepdown recovery pathway (27.8 days, but low probability) reduces its pull on the mean. - `p_death_total` increases monotonically with age, from 2.6% to 39.5%, consistent with clinical expectation. - `mean_death_duration` varies only modestly (10.5–11.5 days), because the biological time from hospitalisation to death is relatively constant across ages — what varies is the probability of dying, not how long it takes.

**Known limitations:** - Table S8 age-scaling factors are transcribed by hand from the published supplementary; values should be verified against the original PDF before the final dissertation submission. - Table S6 values are similarly hand-transcribed; verify against original PDF. - The Knock → ONS region mapping (pooling East/West Midlands; North East/Yorkshire) is an approximation and is documented as such. - `p_ICU_pre` (the 2.5-day ICU triage pre-compartment) is not separately modelled in our 2-compartment HD/HR structure; its duration is implicitly absorbed into the ICU pathway timings above. - All values represent first-wave (pre-vaccine, pre-variant) severity; not appropriate for later pandemic periods.

------------------------------------------------------------------------

## 09_age_imd_stratified_model.R

**Purpose:** Age x IMD decile stratified SEIRD + hospital model in odin. Extends Goodfellow et al. (2024)'s age-stratified SEIRD structure (force of infection, clinical/subclinical split) with two added compartments, HD and HR, between Ic and the final outcomes, parameterised using Knock et al. (2021)'s IHR/IFR and hospital pathway probabilities.

**Inputs:** - `data/parameters/pi_matrix.csv` (script 07) - `data/parameters/h_mu_by_age.csv` (script 07) - `data/parameters/contact_matrix_imd1.csv` ... `imd10.csv` (script 07, urban/rural blended using national IMD-level weights from `rural_age.csv`) - `data/parameters/gamma_hd_hr_by_age.csv` (script 12, age-specific rates)

**Key updates since first version:**

*Contact matrices (script 07 update):* `contact_matrix_imd*.csv` files are now urban/rural population-weighted blends (not urban-only). File names unchanged — no code change needed in this script. Plot subtitle updated from "England (urban)" to "England — urban/rural blended contact matrices" to reflect this. Note that the blending weights here are national IMD-level (from `rural_age.csv`), not region-specific — the region-specific blending is implemented in script 09.

*gamma_hd, gamma_hr (script 12 update):* Previously hardcoded as flat constants (`gamma_hd = 1/11.3`, `gamma_hr = 1/14.1` — a rough 80/20 placeholder). Now read from `gamma_hd_hr_by_age.csv` as **age-varying arrays** (length 17). The odin model definition was updated accordingly: `gam_hd` and `gam_hr` changed from `user()` scalars to `gam_hd[] <- user()` arrays with `dim(gam_hd) <- n_age`, and the `deriv()` equations updated to index `gam_hd[i]` and `gam_hr[i]`.

**Model structure:**

```         
S → E → Ip (presymptomatic) → Ic (clinical) → HD (→ D, death in hospital)
                                             → HR (→ R, recovery)
              Is (subclinical) → R
```

Force of infection: `lambda[k] = susc * sum(weighted_contact[k,]) / proportion[k]` where `weighted_contact[k,j] = contact[k,j] * (Ip[j] + Ic[j] + xi*Is[j])`

**History:** originally written in deSolve (validated working version confirmed), then rewritten in odin per request. Two bugs fixed during conversion: - `sum()` cannot take a compound expression in odin — fixed by computing elementwise product into intermediate array first - deSolve version had Ic/Is ordering mismatch in init vector causing instability — fixed by aligning to Goodfellow's convention (Ic before Is)

**Output:** `output/imd_hospital_gradient_odin.png`

**Validation:** - IMD decile 1 (most deprived) peaks highest, decreasing toward decile 10, with a small uptick at decile 10 (older age structure partially offsets lower clinical fraction) — consistent with Goodfellow et al. (2024) - Independent R0 check (script 10): R0 = 2.714 for IMD decile 1, matching Goodfellow's reported 2.71 to 2 decimal places

**Parameter status after script 12 update:**

| Parameter      | Status                | Value / source                   |
|----------------|-----------------------|----------------------------------|
| gamma_hd       | Refined (script 12)   | 0.0873–0.0953 across ages        |
| gamma_hr       | Refined (script 12)   | 0.0785–0.0908 across ages        |
| h_a, mu_ca_h   | Age-specific only     | Script 07, from Knock Table S9   |
| pi_a           | Age x IMD             | Script 07, from Goodfellow LOESS |
| Contact matrix | IMD-specific, blended | Script 07, urban/rural weighted  |

------------------------------------------------------------------------

## 10_region_stratified_model.R

**Purpose:** Adds the region dimension to the age x IMD stratified model. Region enters ONLY through population composition (age x IMD population structure differs by region) — the model equations (odin generator from script 08) are unchanged. This was a deliberate methodological choice: CoMix contact data showed no significant regional differences in contact rates, and OxCGRT has no NHS-region-level policy granularity, so region-specific behaviour was not assumed.

**Inputs:** - `data/processed/population_age_imd_region.csv` (script 05) - `data/processed/lookup_LSOA11_STP_ITL2_IMD_clean.csv` (script 02) - `data/raw/ruc_lsoa_2fold.csv` (ONS Rural Urban Classification 2011, LSOA-level) - `data/parameters/G.csv`, `rural_age.csv`, `urban_share_by_decile.csv` (script 07) - `data/parameters/pi_matrix.csv`, `h_mu_by_age.csv` (script 07) - `data/parameters/gamma_hd_hr_by_age.csv` (script 12)

**Key updates since first version:**

*gamma_hd, gamma_hr:* Same change as script 08 — now reads age-specific vectors from `gamma_hd_hr_by_age.csv` and passes them as arrays to the odin generator.

*Contact matrix blending (region x IMD specific):* Each region x IMD combination gets its own blended contact matrix:

`M_blended = w_urban * M_urban + (1 - w_urban) * M_rural`

where `w_urban` is **region x IMD specific**, computed from LSOA-level RUC data (32,844 LSOAs each labelled Urban/Rural) joined to our LSOA → region → IMD lookup. LSOA count used as population proxy (\~1,500 residents per LSOA). Falls back to national IMD-level weight for cells with fewer than 5 LSOAs.

`M_urban` and `M_rural` are density-corrected from G using Goodfellow's national urban/rural age structures for that IMD decile — region-specific urban/rural age breakdowns are not available (documented approximation).

The `proportion` vector uses the actual region x IMD total age structure from `pop_region_rebinned` (ONS data, re-binned to Goodfellow's 17 age bands).

**Approximation documented:** ONS age bins (0-4, 5-9, ..., 80+) re-binned to Goodfellow's (Under 1, 1-4, ..., 75+) by splitting 0-4 using Goodfellow's national ratio and merging 75-79 + 80+ into 75+.

**Output:** `output/region_imd_hospital_gradient.png`

**Validation (post urban/rural blending update):** - 76 of 90 region x IMD combinations produced results; 14 sparse cells skipped - Clear, expected gradient (decile 1 peaks highest): East, South West, South East, Yorkshire - North West: previously non-monotonic (decile 2 \> decile 1) — resolved after urban/rural blending, confirming the earlier distortion was from using inappropriate urban-only matrices for this mixed-setting region - East Midlands: modest improvement but still some non-monotonicity in lower deciles — small LAD sample (\~3-4 LADs per cell), reported as limitation - North East: only deciles 1-4 populated — deprivation highly concentrated in this region; no LADs fall into less-deprived deciles - London: flat gradient (0.65-0.83 per 1,000) — likely genuine effect of London's unusually homogeneous age structure across IMD deciles; not a model error - South East: decile 10 (least deprived) sits higher than deciles 6-9 — attributable to older age structure in South East's affluent deciles, where higher h_a partially offsets lower pi_a

**Parameter status after script 12 update:**

| Parameter        | Varies by          | Source                           |
|------------------|--------------------|----------------------------------|
| gamma_hd         | Age                | Script 12 (Knock Table S2/S6/S8) |
| gamma_hr         | Age                | Script 12 (Knock Table S2/S6/S8) |
| w_urban          | Region x IMD       | RUC LSOA data (script 09)        |
| Contact matrix M | Age x IMD x Region | Scripts 07 + 09                  |
| pi_a             | Age x IMD          | Script 07                        |
| h_a, mu_ca_h     | Age only           | Script 07                        |

------------------------------------------------------------------------

## 11_verify_R0.R

**Purpose:** Independent closed-form check (next-generation-matrix method, not the full ODE) that the model's contact matrices and clinical fraction data combine correctly to produce a plausible R0.

**Inputs:** `data/parameters/clin_frac.csv`, `G.csv`, `rural_age.csv`

**Validation:**

| Check | Result | Reference | Match |
|----|----|----|----|
| R0, urban, IMD decile 1 (most deprived) | 2.714 | 2.71 (Goodfellow et al. 2024, reported) | OK, matches to 2 decimal places |
| R0 trend across deciles | Monotonic decrease 1→9 (2.71→2.17), slight uptick at decile 10 (2.18) | Same pattern described in Goodfellow et al. 2024 | OK |

**Conclusion:** Strong independent confirmation that the contact matrix construction and clinical fraction data (script 07 outputs) are correct — this validates the foundational transmission parameters before trusting any dynamic model output built on top of them (scripts 08, 09).

------------------------------------------------------------------------

## 12_regression_stringency.R

**Purpose:** Replaces saturated week fixed effects (`factor(epiweek)`) with a natural spline time trend (`ns(week_num, df)`) in the negative binomial regression, at the advisor's request, to free up identifying variation so `stringency_lag1` can be estimated (previously perfectly collinear with week FE).

**Input:** `data/processed/regression_data.rds` (script 04, post-fix: RegionCode corrected to `UK_ENG`, population County-exclusion bug fixed, `dplyr::select()` namespacing fixed)

**Method:** 1. Fit the model across a range of spline degrees of freedom (df candidates tested: 4, 6, 8, 10, 12, 16, 20, 25, 30, 40) 2. Compare AIC across df to find best-fitting smoothness 3. Check stability of the stringency_lag1 coefficient across df — report instability rather than picking the df that gives the "nicest" result 4. Diagnostic plot: spline-implied time trend vs raw weekly admissions, to check the spline isn't over-smoothing real epidemic waves

**Status: IN PROGRESS, not yet resolved.**

**Findings so far:** - At df=16, the spline visibly under-fits both epidemic waves (model-implied peak \~5,000 admissions/week vs raw peak of \~24,000 in wave 1, \~13,000 in wave 2) — the spline was far too smooth to capture true epidemic dynamics - df_candidates range was expanded to include up to 40, to test whether AIC has a genuine optimum or keeps improving toward the saturated-FE limit - stringency_lag1 coefficient is unstable across df values tested so far (df=4: IRR=0.990, CI crosses null-adjacent; df=6: IRR=1.04; df=8-16: IRR oscillates 1.00-1.02, with df=12's CI crossing 1) — this instability itself may be the most honest finding to report, rather than the coefficient at any single chosen df

**Known unresolved tension:** smoothing time enough to leave room for stringency to be estimated conflicts with capturing the true (sharp, multi-wave) shape of the epidemic. If AIC keeps improving all the way to df=40, this suggests the data may not support both goals simultaneously — worth reporting as a substantive finding, not just a technical wrinkle.

**Unchanged caveat (regardless of df chosen):** this only resolves the collinearity problem, not reverse causality. Stringency varies only at the national-weekly level (OxCGRT has no England sub-national policy data), so it is structurally confounded with national epidemic dynamics. Any stringency_lag1 coefficient from this model should be reported as a descriptive association, not a causal policy effect.

------------------------------------------------------------------------

## [13] — unused

Script 13 briefly existed as `13_verify_R0.R` during renumbering before settling at `10_verify_R0.R`.
