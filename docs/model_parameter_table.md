# Model Parameters — Extended IMD-COVID Model

Note: f, γ, r_c, r_s are durations in days, used directly (not as reciprocals).

## Transmission / progression (Goodfellow 2024, Table 1)

| Symbol | Meaning | Value | Source |
|---|---|---|---|
| p | Transmission probability per contact | 0.06 | Goodfellow 2024 / Davies et al. 2020 |
| M_ak | Daily contacts, age a–k | Varies by age & IMD | Goodfellow 2024, POLYMOD via socialmixr |
| ξ | Relative infectiousness, subclinical | 0.5 | Goodfellow 2024, assumption |
| π_a | Clinical fraction | Varies by age & IMD | Goodfellow 2024, LOESS on Census 2021 + Davies et al. |
| f | Latent period duration | 3 days | Goodfellow 2024 / Davies et al. |
| γ | Preclinical infectious period duration | 2.1 days | Goodfellow 2024 / Davies et al. |
| r_c | Clinical infectious period duration | 2.9 days | Goodfellow 2024 / Davies et al. |
| r_s | Subclinical infectious period duration | 5 days | Goodfellow 2024 / Davies et al. |
| μ_sa | Subclinical mortality probability | 0 | Goodfellow 2024, assumption |

## Severity / hospitalisation (extension, combines both papers)

| Symbol | Meaning | Value | Source |
|---|---|---|---|
| IHR_a | Infection hospitalisation ratio by age | data/parameters/ihr_ifr_by_age_knock2021.csv | Knock 2021, Table S9 |
| IFR_a | Infection fatality ratio by age | data/parameters/ihr_ifr_by_age_knock2021.csv | Knock 2021, Table S9 |
| h_a | P(hospitalised \| clinical) | IHR_a / π_a | Derived |
| μ_ca_h | P(death \| hospitalised) | IFR_a / IHR_a | Derived, Knock 2021 |
| μ_ca_g | P(death \| clinical, not hospitalised) | 0 | Assumption |
| γ_hd | Rate, HD → D | 0.0885 (placeholder) | See calculation below |
| γ_hr | Rate, HR → R | 0.0709 (placeholder) | See calculation below |

## γ_hd / γ_hr calculation (placeholder method)

Weighted average of Knock (2021) Table S2 pathway durations, using a rough
80% general-ward / 20% ICU-pathway split (not derived from actual branch
probabilities — needs replacing, see TODO).

γ_hd: 0.8 × 10.3 (general ward death) + 0.2 × 13.45 (avg of ICU direct death
11.8, and ICU→stepdown death 11.1) ≈ 11.3 days → γ_hd = 1/11.3 ≈ 0.0885

γ_hr: 0.8 × 10.7 (general ward recovery) + 0.2 × 27.8 (ICU→stepdown recovery
15.6+12.2) ≈ 14.1 days → γ_hr = 1/14.1 ≈ 0.0709

## TODO — parameters to change for the age × IMD × region stratified model

- **π_a**: replace single value (0.55) with full age × IMD-decile array from
  Goodfellow's LOESS output (his repo `/data/` folder)
- **IHR_a, IFR_a / h_a, μ_ca_h**: replace single placeholder (h=0.036,
  μ_ca_h=0.1) with full age-stratified table (already have:
  `ihr_ifr_by_age_knock2021.csv`); h_a and μ_ca_h become IMD-decile-dependent
  too once π_a varies by decile
- **M_ak**: replace single average contact value (M=11) with full age-specific
  contact matrix from Goodfellow's repo (POLYMOD via socialmixr, density-
  corrected by IMD/age population structure)
- **γ_hd, γ_hr**: replace placeholder 80/20 weighting with proper probability-
  weighted aggregation using actual age-varying branch probabilities (Knock
  2021 Table S6/S8)
- **N / initial conditions**: replace single 100,000 toy population with real
  region × IMD × age population (`population_age_imd_region.csv`, script 05)
- **Region dimension**: enters via population composition only (age × IMD
  distribution per region), not via contact matrix or policy stringency —
  decision made earlier based on lack of region-level data for either
