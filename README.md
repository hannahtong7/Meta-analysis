# Meta-analysis: Thermal Sensitivity of Ocean Phytoplankton – Meta-Analysis

This repository contains all code and data for a meta-analysis investigating how phytoplankton growth responds to ocean warming, with a focus on variation across major taxonomic groups, latitudes, and thermal origins. The analysis quantifies growth responses (as log response ratios) under high-emission scenarios (+3 to +5 °C) and examines the influence of moderators such as **Taxa**, **Latitude**, and **Control Temperature**.

**Files:**

-   `Meta-analysis-final-code.R` – main analysis script

-   `Meta-analysis.csv` – original raw dataset

-   `Meta-analysis_filtered_3to5.csv` – dataset filtered for +3 to +5 °C treatments

-   `Meta-analysis_filtered_3to5_with_lnRR.csv` – final dataset including effect sizes (lnRR)

## Objectives

-   Quantify the **overall effect** of warming (3–5 °C) on phytoplankton growth.

-   Examine **strain-level**, **taxonomic**, and **latitudinal variation** in responses.

-   Assess the role of **thermal origin (control temperature)** in modulating thermal sensitivity.

-   Visualize results via **forest plots**, **moderator regressions**, and **taxon-specific summaries**.

-   Evaluate **publication bias** using funnel plots and Egger's test.
