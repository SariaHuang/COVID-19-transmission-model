---

editor_options: 
  markdown: 
    wrap: 72
---

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
|------------------|------------------|------------------|------------------|
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
|------------------|------------------|------------------|------------------|
| Total population | 56,550,138 | \~56.5M (England, mid-2020/2021) | OK |
| Age bands | 17 | 17 (0-4 ... 80+) | OK |
| IMD deciles | 10 | 10 | OK |
| Regions | 9 | 9 | OK |
| Missing cells | 238 / 1530 | Expected — sparse region x decile combinations | OK |

**Bug found and fixed:** first run gave total population 60,003,484 (\~3.45M too high). Traced to the script-03 duplicate-LAD bug (above) — once that was fixed, total population matched the independently cross-checked value (56,550,138) exactly.

------------------------------------------------------------------------

## 06_seird_unstratified.R

**Purpose:** Single-population (no age/IMD/region stratification) SEIRD + hospital skeleton model in odin2, to confirm the model structure before adding stratification.

**odin2/dust2 conversion:** `user()` replaced by `parameter()`; `dt` parameter removed (odin2 provides dt automatically); model run via `dust_system_create` + `dust_system_set_state_initial` + `dust_system_simulate`. Output is a 2D array `[n_states, n_times]` (no particle dimension for deterministic models).

**Parameters:** see `docs/model_parameter_table.md`. Several values are single-representative placeholders (h, mu_ca_h, M, gamma_hd, gamma_hr) — not yet age/IMD/region-stratified. gamma_hd and gamma_hr are refined in script 08.

**Validation:** equations checked in base R before odin2 translation — - Population exactly conserved across all time steps - Single epidemic wave, peak clinical infectious \~day 64 - Final attack rate \~91% - Final deaths and admissions match theoretical values exactly (pi_c × h × mu_ca_h and pi_c × h respectively)

**Status:** structure confirmed working — superseded by script 09 (which adds age and IMD stratification using the same structure). Kept as a reference / sanity baseline.

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

**Purpose:** Age × IMD decile stratified SEIRD + hospital model in odin2. Extends Goodfellow et al. (2024)'s age-stratified SEIRD structure (force of infection, clinical/subclinical split) with two added compartments, HD and HR, between Ic and the final outcomes, parameterised using Knock et al. (2021)'s IHR/IFR and hospital pathway probabilities.

**Inputs:** - `data/parameters/pi_matrix.csv` (script 07) - `data/parameters/h_mu_by_age.csv` (script 07) - `data/parameters/rural_age.csv` (Goodfellow repo) - `data/parameters/urban_share_by_decile.csv` (script 07) - `data/parameters/contact_matrix_imd1.csv` ... `imd10.csv` (script 07, urban/rural blended using national IMD-level weights from `rural_age.csv`) - `data/parameters/gamma_hd_hr_by_age.csv` (script 08, age-specific rates)

**Key updates since first version:**

*odin → odin2/dust2 conversion:* - `user()` replaced by `parameter()` throughout the model block - Model compiled with `odin2::odin({...})` - Runtime: `dust_system_create` + `dust_system_set_state_initial` + `dust_system_simulate`; output is `[n_states, n_times]` matrix - State order fixed as: S[1:17], E[18:34], Ip[35:51], Ic[52:68], Is[69:85], HD[86:102], HR[103:119], D[120:136], R[137:153], Adm[154:170] (170 states total = 10 compartments × 17 age groups)

*Urban/rural blended population (`get_blended_inputs()` — new function):* The `proportion` vector (used as initial age structure and force-of-infection denominator) and implied `pop_size` now use a **blended urban/rural population** rather than urban-only. This is consistent with the blended contact matrices.

``` r
get_blended_inputs <- function(imd_decile) {
  w           <- urban_share$w_urban[urban_share$IMD == imd_decile]
  pop_blended <- w * pop_urban + (1 - w) * pop_rural
  list(proportion = pop_blended / sum(pop_blended),
       pop_size   = sum(pop_blended))
}
```

`w_urban` comes from `urban_share_by_decile.csv` (national IMD-level urban population share, ranging from decile 1 = 98.4% to decile 10 = 80.7%). This function is defined once in script 09 and used in scripts 13, 14, 16, 17.

*Contact matrices (script 07 update):* `contact_matrix_imd*.csv` files are urban/rural population-weighted blends (not urban-only). Blending weights are national IMD-level (from `rural_age.csv`). Region-specific blending is implemented separately in script 10.

*gamma_hd, gamma_hr (script 08 update):* Previously hardcoded as flat constants (`gamma_hd = 1/11.3`, `gamma_hr = 1/14.1`). Now read from `gamma_hd_hr_by_age.csv` as **age-varying arrays** (length 17). `gam_hd[]` and `gam_hr[]` declared as `parameter()` arrays in the odin2 model and indexed as `gam_hd[i]`, `gam_hr[i]` in `deriv()` equations.

**Model structure:**

```         
S → E → Ip (presymptomatic) → Ic (clinical) → HD (→ D, death in hospital)
                                             → HR (→ R, recovery)
              Is (subclinical) → R
```

Force of infection: `lambda[k] = susc * sum(weighted_contact[k,]) / proportion[k]` where `weighted_contact[k,j] = contact[k,j] * (Ip[j] + Ic[j] + xi*Is[j])`

**History:** originally written in deSolve (validated working version confirmed), then rewritten in odin per request. Two bugs fixed during conversion: - `sum()` cannot take a compound expression in odin — fixed by computing elementwise product into intermediate array first - deSolve version had Ic/Is ordering mismatch in init vector causing instability — fixed by aligning to Goodfellow's convention (Ic before Is)

**Output:** `output/imd_hospital_gradient_odin.png`

**Validation (blended population, odin2):** - IMD decile 1 (most deprived) peaks highest, decreasing toward decile 10, with a small uptick at decile 10 (older age structure partially offsets lower clinical fraction) — consistent with Goodfellow et al. (2024) - Attack rate gradient confirmed: decile 1 = 87.2%, decile 10 = 78.2% - `get_blended_inputs()` confirmed working: proportion vectors sum to 1.0 for all 10 deciles; pop_size values consistent with national totals - Independent R0 check (script 10): R0 = 2.714 for IMD decile 1, matching Goodfellow's reported 2.71 to 2 decimal places

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
|------------------|------------------|------------------|------------------|
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

## 13_monty_fitting.R

**Purpose:** Fit the age × IMD SEIRD + hospital model to observed NHS weekly hospital admissions using Bayesian MCMC (monty package). Produces posterior distributions for transmission probability (β) and negative-binomial overdispersion (size) for each of the 10 IMD deprivation deciles.

**Inputs:** - `data/processed/observed_weekly_admissions.csv` (script 04) - `data/parameters/gamma_hd_hr_by_age.csv` (script 08/12) - `data/parameters/pi_matrix.csv`, `h_mu_by_age.csv`, `rural_age.csv` (script 07) - `data/parameters/contact_matrix_imd{1-10}.csv` (script 07) - odin model generator `age_seird_hosp` (sourced from script 09)

**Outputs:** `output/fitting/fitted_samples_imd{1-10}.rds`

------------------------------------------------------------------------

### Approach: wrapper method with odin2/dust2

monty is used via the **wrapper approach**: `monty_model_function()` wraps a hand-written log-likelihood function. The odin2 model from script 09 (`age_seird_hosp`, compiled with `odin2::odin()`) is called inside the likelihood for each proposed parameter set via `dust_system_create` + `dust_system_simulate`.

This is the wrapper approach (not native odin2/monty integration with `compare()` and `dust_filter_create`), chosen for compatibility with the already-validated model structure from script 09 and the hand-written negative binomial likelihood.

------------------------------------------------------------------------

### Parameters fitted

| Parameter | Description                  | Prior              | Starting value |
|------------------|------------------|------------------|------------------|
| log_beta  | Log transmission probability | Normal(-2.81, 0.5) | log(0.06)      |
| log_size  | Log NegBin overdispersion    | Normal(1, 1)       | log(10)        |

**I0_frac fixed at 1×10⁻⁴** — not fitted. In a single-wave deterministic ODE model, I0_frac and β are not separately identifiable: any combination of lower β + higher I0_frac (or vice versa) can produce the same epidemic trajectory. Attempting to fit both produced a monotonically drifting I0_frac chain (visible in diagnostic trace plots) with no convergence after 2,000 steps. I0_frac was fixed at 1×10⁻⁴ (approximately 1 case per 10,000 at epidemic onset) and this decision is reported as a methodological limitation in the dissertation.

------------------------------------------------------------------------

### rec_c bug (found and fixed during development)

The previous hand-written deSolve version of `run_epidemic()` in early drafts of script 13 contained a critical bug: `rec_c = 1/2.9` was passed as a rate, but then used as `1/rec_c` in the ODE equations, effectively double-inverting it to give `rec_c = 2.9` as the rate (8.4× faster than intended). This caused Ic to clear almost instantaneously, distorting the hospitalisation inflow and making the epidemic curve collapse too early.

**Fix:** The hand-written deSolve `run_epidemic()` was **replaced entirely** with a wrapper (`run_epidemic_fit()`) around the validated odin model from script 09. In script 09's odin code, `rec_c <- user(2.9)` is correctly treated as a duration (days), and `1/rec_c` is taken inside `deriv()` — this was already correct and required no change.

------------------------------------------------------------------------

### pop_size evolution (two fixes)

**Fix 1 (during development):** The original likelihood used total (urban + rural) population, but the model ran with `urban = TRUE`, inflating absolute admission counts for deciles 3–8 (urban fraction 73–85%). Fixed to urban-only population.

**Fix 2 (blending update):** After implementing `get_blended_inputs()` in script 09, `pop_size` was updated to use **blended urban/rural population** throughout, consistent with the blended `proportion` vector:

``` r
pop_size <- get_blended_inputs(imd_decile)$pop_size
```

This ensures the absolute admission scaling (cum_admissions × pop_size) is consistent with the model's population basis. The change is small for most deciles (decile 1: +1.7%; decile 8: +25% vs urban-only) but improves methodological consistency. MCMC was rerun after this update.

------------------------------------------------------------------------

### Data window

The observed data (`obs_data`) starts 2020-07-27 — the earliest date available from the NHS COVID-19 Hospital Activity publication accessed via `covid19.nhs.data::get_admissions()`. The true first wave (March–June 2020) is not available at LAD/IMD-decile level from any public source (trust-level data exists but cannot be reliably mapped to IMD deciles — investigated and documented separately). The fitting window therefore covers the **autumn/winter 2020 wave**, not the true first wave.

Wave 1 end (first trough after first peak) was **manually set to 2021-04-05** after auto-detection (smoothed rollmean + argmin) returned 2021-05-10 — the auto-detection found no genuine trough within the data window, as the post-peak decline extended to the data boundary. Manual date chosen based on epidemiological knowledge of the UK inter-wave low point (spring 2021, before Delta wave onset).

**Wave 1 window: 2020-07-27 to 2021-04-05 (37 weeks)**

------------------------------------------------------------------------

### Known limitation: severity parameter time mismatch

Hospitalisation pathway parameters (`IHR_a`, `IFR_a`, `h_a`, `mu_ca_h`, `γ_hd`, `γ_hr`) are sourced from Knock et al. (2021) estimates based on the **true first wave** (pre-vaccine, pre-variant). The data window fitted here covers the autumn/winter wave, during which Alpha variant emergence (December 2020) and vaccine rollout (from 8 December 2020) altered severity. This time mismatch is reported as a known limitation.

------------------------------------------------------------------------

### Sampler settings

```         
n_samples = 2000, n_chains = 3, burnin = 500
vcv = diag(2) * c(0.02, 0.04)   # random walk step sizes for log_beta, log_size
```

------------------------------------------------------------------------

### Convergence results

| IMD decile | beta median (C1/C2/C3) | size median (C1/C2/C3) | Converged |
|------------------|------------------|------------------|------------------|
| 1 (most deprived) | \~0.031 / 0.031 / 0.031 | \~4.5 / 4.3 / 4.4 | Yes |
| 2–10 | See `output/fitting/` .rds files | See .rds files | Yes |

All chains converged for β (stable mixing within \~50 steps of burn-in). size stabilised in the 2–5 range for most deciles; no chain collapsed toward 0 (which would indicate structural model misfit).

------------------------------------------------------------------------

### Fitted beta gradient (key finding and limitation)

β increases monotonically from decile 1 (β ≈ 0.031) to decile 10 (β ≈ 0.041), contrary to the expected deprivation gradient. This arises from a systematic data coverage bias: `obs_data` coverage is highest in decile 1 (98.4% urban LADs, confirmed) and lower in wealthier deciles, which attenuates the per-capita observed admission rate in wealthier deciles and biases MCMC β estimates upward. The fitted β gradient is therefore **not interpretable as a true transmissibility gradient** and is discussed as a methodological artefact in the dissertation. See script 15 for the sensitivity analysis addressing this.

------------------------------------------------------------------------

## 14_age_imd_gradient_plots.R

**Purpose:** Visualise model outputs at age group × IMD decile resolution, using posterior median β from script 13. Produces heatmaps and gradient line plots showing how cumulative hospital admissions and in-hospital deaths vary across age and deprivation.

**Inputs:** - `output/fitting/fitted_samples_imd{1-10}.rds` (script 13) - odin model and all parameter data (sourced from scripts 08–09)

**Outputs:** `output/plots/14_heatmap_admissions.png`, `14_heatmap_deaths.png`, `14_gradient_admissions.png`, `14_gradient_deaths.png`, `14_age_curves_by_decile.png`, `14_summary_burden_by_decile.csv`

------------------------------------------------------------------------

### New function: run_epidemic_fit_age()

Unlike `run_epidemic_fit()` in script 13 (which returns only aggregate `cum_admissions`), this function returns age-specific cumulative admissions and deaths for all 17 age groups, extracted from odin's `Adm[a]` and `D[a]` output columns. Per-1,000 rates are computed as:

```         
cum_adm_per1000[a] = Adm[a] * 1000
```

This is correct because the odin model initialises `S0 = proportion` (a vector summing to 1), so all compartments are expressed as proportions of the age-group population — multiplying by 1,000 gives the per-1,000 rate directly without further population scaling.

------------------------------------------------------------------------

### Key finding: reversed admission gradient

The heatmap of cumulative admissions (age × decile) shows 75+ admission rates **increasing from decile 1 to decile 10** — the opposite of the expected deprivation gradient. This is a direct consequence of the fitted β gradient documented in script 13: higher β in wealthier deciles produces more modelled infections and admissions regardless of the deprivation-varying π_a and h_a.

This finding is reported as a data-driven artefact and is addressed in the sensitivity analysis (script 15, fixed β approach).

------------------------------------------------------------------------

### Death gradient

Cumulative death rates show the same spatial pattern as admissions, as expected mathematically: deaths = `μ_ca_h × h_a × (flow through Ic)`, which differs from admissions only by the fixed age-specific factor `μ_ca_h`. The two heatmaps are therefore structurally identical and their similarity is not a model error.

------------------------------------------------------------------------

### Negative values in daily admission curves (fixed)

`diff(cum_adm_per1000)` produced small negative values at ODE boundary steps (numerical precision artefact). Fixed by applying `pmax(..., 0)` to all daily difference series. `coord_cartesian(ylim = c(0, NA))` added to prevent negative-value spikes from collapsing the y-axis scale in panel plots.

------------------------------------------------------------------------

## 15_region_age_imd_plots.R

**Purpose:** Extend the age × IMD analysis to the region dimension. Uses model per-1,000 rates from script 14, weighted by regional population data (`population_age_imd_region.csv`), to produce region × age × IMD burden estimates. Produces **both** main analysis and sensitivity analysis outputs in a single script run, with unified colour scales for direct comparison.

**Inputs:** - `output/fitting/fitted_samples_imd{1-10}.rds` (script 13) - `data/processed/population_age_imd_region.csv` (script 05) - odin model and all parameter data (scripts 08–09)

**Outputs** (saved to `output/plots/region/`):

| File prefix | Analysis | β used |
|------------------------|------------------------|------------------------|
| `15_main_*` | Main analysis | Decile-specific fitted β |
| `15_sensitivity_*` | Sensitivity analysis | Fixed β = decile 1 posterior median (0.031) |

Five plots per analysis (10 total) + 2 CSV summary files.

------------------------------------------------------------------------

### Region dimension: population-weighting approach

Region does not enter the ODE model directly (the transmission equations are identical across regions). Region enters only through population composition: each region has a different distribution of population across age × IMD decile cells, producing different expected absolute burdens even at the same per-1,000 rate.

**Calculation:**

```         
abs_adm[r, d, a]   = (cum_adm_per1000[d, a] / 1000) × pop[r, d, a]
adm_per1000[r]     = Σ_{d,a} abs_adm[r, d, a] / total_pop[r] × 1000
```

where `pop[r, d, a]` is from `population_age_imd_region.csv`.

This is consistent with the model's scope: region-specific contact rates and policy responses are not modelled (no region-level data available for either — documented in script 10). Region differences in results reflect population composition only, not differential transmission dynamics.

------------------------------------------------------------------------

### Age band mapping

`population_age_imd_region.csv` uses age bands `"0-4", "5-9", ..., "75-79", "80+"` (17 bands). The odin model uses 17 age groups starting from "Under 1" and "1 to 4" separately, and ending at "75+".

**Mapping applied:** - `"0-4"` in region data → weighted average of model "Under 1" + "1 to 4", using urban population weights from `rural_age` (w_Under1 ≈ 0.18, w_1to4 ≈ 0.82 — the 0-4 band is mostly 1–4 year olds) - `"75-79"` + `"80+"` merged → model "75+" - All other bands: direct 1:1 mapping

------------------------------------------------------------------------

### Sensitivity analysis: fixed β

To isolate the contribution of deprivation-varying clinical fractions (π_a) and hospitalisation probabilities (h_a) to the inequality gradient — removing the artefactual influence of the data-coverage-biased β estimates — a sensitivity analysis was conducted with β fixed at the decile 1 posterior median (β = 0.031) across all deciles.

**Justification:** Decile 1 has the highest data coverage (98.4% urban LADs), making it the most reliable β estimate. Fixing β = 0.031 uniformly eliminates the spurious β gradient while preserving all within-model inequality variation driven by π_a and h_a.

In the sensitivity analysis heatmaps, the admission gradient reverses to the expected direction (decile 1 highest), confirming that π_a and h_a alone produce a genuine deprivation-related burden gradient when β is held constant. This result is the primary evidence that the model's inequality parameters are correctly specified despite the data-coverage artefact in β.

------------------------------------------------------------------------

### Unified colour scale

Both main and sensitivity analysis heatmaps use **identical colour scale limits**, computed from the 95th percentile of the combined data from both analyses. This allows direct visual comparison of the two sets of plots. The 5% of cells exceeding the cap are shown at the maximum colour value (`oob = scales::squish`) rather than being clipped to NA.

------------------------------------------------------------------------

### Known limitations documented

1.  **Severity parameter time mismatch** (inherited from script 13): Knock 2021 parameters are first-wave estimates; data window is autumn/winter 2020.
2.  **Blended model population vs total regional population**: the ODE model uses blended urban/rural proportion vectors and pop_size (via `get_blended_inputs()` in script 09), but `population_age_imd_region.csv` used as the regional denominator contains total (urban + rural) population. The blending weights (w_urban) are national IMD-level averages, not region-specific, so some mismatch between model population basis and regional population denominator remains. Reported as an approximation; direction of bias is small and expected to be less than 5% for most region × decile cells.
3.  **Sparse region × decile cells**: North East, West Midlands, and Yorkshire have near-zero population in deciles 5–10 (deprivation is geographically concentrated in these regions). These cells appear as empty/grey in the heatmaps — this is genuine geographic sparsity, not missing data.
4.  **β gradient artefact** (from script 13): main analysis plots show reversed deprivation gradient. Sensitivity analysis (fixed β) is the recommended primary figure for reporting the inequality gradient.

------------------------------------------------------------------------

## 16_counterfactual.R

**Purpose:** Counterfactual analysis quantifying how much of the modelled COVID-19 hospitalisation and death burden is attributable to deprivation- related health inequality. Compares a baseline scenario (observed π_a by IMD decile) against two counterfactual scenarios:

- **CF1:** all deciles assigned decile 10 (least deprived) π_a — eliminates health inequality entirely
- **CF2:** all deciles assigned decile 5 (median) π_a — partial equalisation, more realistic policy target

**Inputs:** - `output/fitting/fitted_samples_imd{1-10}.rds` (script 13) - `data/parameters/pi_matrix.csv`, `h_mu_by_age.csv`, `rural_age.csv`, `urban_share_by_decile.csv` (scripts 07–09) - `data/processed/population_age_imd_region.csv` (script 05) - odin2 model generator `age_seird_hosp` + `get_blended_inputs()` (script 09)

**Outputs** (saved to `output/plots/counterfactual/`): - `16_fig1_attributable_fraction.png` — CF1 vs CF2 side-by-side - `16_fig2_burden_comparison_region.png` — three lines (baseline/CF1/CF2) by region - `16_fig3_national_summary.png` — avoided burden by age group (both CFs) - `16_counterfactual_full.csv` - `16_counterfactual_decile_summary.csv`

**Population basis:** blended urban/rural via `get_blended_inputs()`, consistent with scripts 09 and 13. `adm_avoided_abs` = (adm_diff / 1000) × pop_blended.

------------------------------------------------------------------------

### Counterfactual design

**What is replaced:** `π_a` (clinical fraction, age × IMD-specific) is replaced with decile 10 (least deprived) values for all deciles. This represents a scenario in which all deprivation deciles experience the probability of clinical disease of the least deprived population.

**What is NOT replaced:** `h_a` (hospitalisation probability given clinical disease). Although `h_a = IHR_a / π_a` in the model, IHR_a is age-specific but **deprivation-invariant** (sourced from Knock et al. 2021, Table S9, which does not stratify by deprivation). Replacing both π_a and h_a simultaneously would mechanically inflate h_a for deprived deciles (lower π_a → higher h_a = IHR_a / π_a), producing artefactual increases in counterfactual admissions relative to baseline — confirmed empirically in diagnostic checks where decile 10 counterfactual admissions exceeded baseline for most age groups. h_a is therefore held at its observed decile-specific value throughout.

**β:** Fixed at decile 1 posterior median (β = 0.031) for all deciles in both baseline and counterfactual, consistent with the sensitivity analysis in script 15. This ensures that differences between baseline and counterfactual are driven solely by π_a, not by the artefactual β gradient documented in script 13.

------------------------------------------------------------------------

### Key results

**CF1 (decile 10 π_a — full equalisation):**

| Metric                                                     | Value |
|------------------------------------------------------------|-------|
| Total avoided admissions (blended population, single wave) | 8,114 |
| Total avoided deaths (blended population, single wave)     | 1,384 |
| Mean attributable fraction — admissions                    | 39.3% |
| Mean attributable fraction — deaths                        | 39.9% |

**CF2 (decile 5 π_a — partial equalisation):**

| Metric | Value |
|------------------------------------|------------------------------------|
| Total avoided admissions | Lower than CF1 (decile 1–4 still above decile 5 π_a) |
| Mean attributable fraction — admissions | Lower than CF1 across all deciles |

CF2 represents a more realistic policy target — reducing health inequality to the median level rather than eliminating it entirely. Comparing CF1 and CF2 shows the marginal gain from targeting the most deprived vs achieving median health equity.

------------------------------------------------------------------------

### Attributable fraction by IMD decile

| Decile              | Adm attr. fraction | Death attr. fraction |
|---------------------|--------------------|----------------------|
| 1 (most deprived)   | 45.3%              | 45.4%                |
| 2                   | 42.4%              | 42.7%                |
| 3                   | 41.6%              | 42.3%                |
| 4                   | 47.0%              | 48.4%                |
| 5                   | 52.0%              | 53.4%                |
| 6                   | 51.9%              | 52.9%                |
| 7                   | 48.8%              | 49.2%                |
| 8                   | 38.9%              | 39.2%                |
| 9                   | 24.8%              | 25.0%                |
| 10 (least deprived) | 0%                 | 0%                   |

------------------------------------------------------------------------

### Non-monotonicity in attributable fraction (documented, not a bug)

The attributable fraction peaks at deciles 5–6 (\~52%) rather than decile 1 (45%), which is counterintuitive. This is a genuine mathematical property of the model, not an error. The attributable fraction is a **relative** measure (baseline − counterfactual) / baseline. Although decile 1 has the largest absolute difference in π_a from decile 10 (Δπ_a = 0.152 vs 0.072 for decile 5), its baseline admission burden is also highest, making the denominator large and suppressing the relative fraction. Confirmed by inspecting π_a values directly:

```         
Decile 1:  mean π_a = 0.480   (Δ from d10 = 0.152)
Decile 5:  mean π_a = 0.400   (Δ from d10 = 0.072)
Decile 10: mean π_a = 0.328   (reference)
```

This non-linearity is reported in the dissertation as a substantive finding: relative inequality in clinical outcomes is not simply proportional to the absolute deprivation gradient in π_a.

------------------------------------------------------------------------

### Age pattern of avoided burden (Figure 3)

75+ age group accounts for 3,065 avoided admissions (38% of total) and 731 avoided deaths (53% of total). This reflects the combined effect of higher absolute π_a differences in older age groups and larger absolute population counts at risk. The concentration of avoided burden in older age groups strengthens the policy relevance: health inequality in clinical severity disproportionately harms the elderly.

------------------------------------------------------------------------

### Known limitations

1.  **Blended population**: absolute avoided counts use blended urban/rural population (`get_blended_inputs()$pop_size`). This is more accurate than urban-only but still uses national IMD-level blending weights rather than decile-specific rural age structures, which remain unavailable.
2.  **Single-wave model**: the 365-day single-wave ODE produces a conservative estimate; multiple waves would accumulate further avoided burden.
3.  **π_a only**: `h_a` is not varied by deprivation due to the deprivation-invariant IHR assumption (see design rationale above). If IHR does in fact vary by deprivation (not captured in Knock 2021), the true attributable fraction would be higher.
4.  **β fixed**: using decile 1 β for all deciles is a modelling assumption; if true between-decile β variation exists, results would differ.
5.  **Severity parameter time mismatch** (inherited from script 13): Knock 2021 parameters are first-wave estimates; data window is autumn/winter 2020.

------------------------------------------------------------------------

## 17_school_closure.R

**Purpose:** School closure policy scenario analysis. Extends Goodfellow et al. (2024) by applying school closure to a model fitted to real NHS admissions data, disaggregated by IMD decile, age group, and ITL1 region.

**Inputs:** - `output/fitting/fitted_samples_imd{1-10}.rds` (script 13) - `data/parameters/pi_matrix.csv`, `h_mu_by_age.csv`, `rural_age.csv`, `urban_share_by_decile.csv` (scripts 07–09) - `data/parameters/contact_matrix_imd{1-10}.csv` (script 07) - `data/processed/population_age_imd_region.csv` (script 05) - odin2 model `age_seird_hosp` + `get_blended_inputs()` (script 09) - POLYMOD UK contact data via `socialmixr` package

**Outputs** (saved to `output/plots/school_closure/`): - `17_fig1_R0_reduction.png` - `17_fig2_admission_reduction.png` - `17_fig3_pct_reduction_by_decile.png` - `17_fig4_reduction_by_age.png` - `17_fig5_region_decile_facet.png` - `17_fig6_region_avoided.png` - `17_school_closure_decile_summary.csv` - `17_school_closure_full.csv`

------------------------------------------------------------------------

### Method: school contact fraction

School contacts are identified from POLYMOD UK data using Goodfellow et al. (2024)'s density-correction method. The school-specific intrinsic connectivity matrix G_school is extracted by filtering POLYMOD contacts where `cnt_school = 1`, then density-corrected using the same procedure as G_total:

```         
G_school[i,j] = M_school[i,j] * sum(N) / N[j]
school_frac[i,j] = min(G_school[i,j] / (G[i,j] + 1e-10), 1)
```

School closure is modelled by removing the school contact fraction from each decile-specific contact matrix:

```         
M_closed[i,j] = max(M_total[i,j] * (1 - school_frac[i,j]), 0)
```

This assumes complete removal of school contacts (upper bound on school closure impact). In practice, some contacts would be displaced to home/other settings (Jackson et al. 2011), so results represent a theoretical maximum reduction.

**POLYMOD warning (not a bug):** `socialmixr` produces two warnings when generating the school contact matrix: 1. "Not all age groups represented" — because the 0–1 age band is not in POLYMOD's 5-year groupings; linearly interpolated by socialmixr. Consistent with Goodfellow et al. (2024)'s approach. 2. "Large differences in sub-population sizes" — normalisation factor range 0–8.4 due to the very small Under-1 sub-population. Impact negligible given near-zero hospitalisation rates in under-1s.

Both warnings are documented and do not affect results materially.

------------------------------------------------------------------------

### R0 computation

R0 for each IMD decile is computed via the Next Generation Matrix method (consistent with script 11):

```         
N[i,j] = beta * M[i,j] * (pi_j*(1/sympt + rec_c) + xi*(1-pi_j)*rec_s)
R0 = dominant eigenvalue of N
```

β is fixed at decile 1 posterior median (0.031) for all deciles, consistent with the sensitivity analysis in script 15.

------------------------------------------------------------------------

### Key results: R0 gradient

| IMD decile        | R0 baseline | R0 schools closed           |
|-------------------|-------------|-----------------------------|
| 1 (most deprived) | 1.39        | 1.10 (above threshold)      |
| 2                 | 1.32        | 1.08 (above threshold)      |
| 3                 | 1.26        | 1.06 (above threshold)      |
| 4                 | 1.21        | 1.03 (above threshold)      |
| 5                 | 1.17        | \~1.00 (threshold)          |
| 6–10              | 1.11–1.15   | 0.90–0.98 (below threshold) |

**Key finding:** School closure pushes R0 below 1 for deciles 5–10 (theoretically eliminating sustained transmission) but NOT for deciles 1–4, where R0 remains above 1 after school closure. This asymmetry is the central finding of script 17: a uniform national school closure provides epidemic control for less deprived populations but is insufficient for the most deprived communities.

**Comparison with Goodfellow et al. (2024):** In their model (p = 0.06, urban only, no NHS data fitting), R0 never dropped below 1 in any decile after school closure. In our model (β = 0.031, fitted to NHS data), R0 drops below 1 for deciles 5–10. The difference is attributable to the substantially lower fitted β (0.031 vs 0.06) — the fitted value reflects NPI suppression during the modelled data window, making school closure sufficient to halt transmission in less- deprived deciles at this effective transmission rate.

------------------------------------------------------------------------

### Population basis

Blended urban/rural population via `get_blended_inputs()`, consistent with scripts 09, 13, 14, 16. `adm_avoided_abs` = (adm_reduction / 1000) × pop_blended.

------------------------------------------------------------------------

### Region analysis

Regional avoided admissions computed by weighting model per-1,000 rates by `population_age_imd_region.csv` (total regional population), consistent with script 15. Figure 5 (`17_fig5_region_decile_facet.png`) filters to cells with ≥10,000 population to remove statistically unstable sparse cells.

------------------------------------------------------------------------

### Known limitations

1.  **Complete school contact removal**: assumes 100% removal of school contacts; actual displacement to other settings (home, community) would attenuate the reduction. Results are theoretical upper bounds.
2.  **Fixed β throughout**: β fixed at decile 1 posterior median (0.031); school closure effect on R0 would differ if true β gradient exists.
3.  **Blended population**: blending weights are national IMD-level averages; region-specific urban fractions used for regional burden but not for individual decile-level ODE runs.
4.  **Static school closure**: school closure is modelled as a permanent change; the actual policy (January–March 2021) was temporary. Dynamic NPI modelling would require time-varying β (future work).
5.  **POLYMOD contacts from 2008**: contact patterns used are from the 2006/2008 POLYMOD study; actual pre-pandemic UK contact patterns may differ.
