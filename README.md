# MS252: **"Effects of ambient temperature on road traffic mortality in Latin America: Individual- and city-level variations across 272 SALURBAL cities"**

This repository contains the full modeling pipeline used for our multi-city time-stratified case-crossover analysis of ambient temperature and road-traffic deaths across 272 cities in Latin America. The workflow covers data processing, DLNM modeling, stratified and interaction analyses, pooled effect estimation using Rubin’s Rule, sensitivity checks, and calculation of Excess Death Fractions (EDF).
---

## 📁 Script Overview

| Script | Description |
|--------|-------------|
| `00_data wrangling_1.R` | Loads and reshapes city-level 100-imputatoins road-traffic mortality data. |
| `01_data_wrangling_2.R` | Aggregates previously processed city-imputation CSV files into imputation-level datasets. |
| `02_median_imputation_main_effects.R` | Runs DLNM models using median-imputed data across all cities; estimates main effects. |
| `03_rubin_imputation_main_effect.R` | Applies Rubin’s Rules to combine results from 100 imputed datasets; produces pooled estimates for main effects. |
| `04_sensitivity_knots.R` | Conducts sensitivity analyses by testing alternative spline knot specifications in DLNM models. |
| `05_subgroup_sex.R` | Stratifies analysis by sex (male/female). |
| `06_subgroup_age.R` | Stratifies analysis by age group (<=9, 10-19, 20–39, 40–64, 65+). |
| `07_subgroup_mode.R` | Stratifies analysis by transport mode (e.g., pedestrian, motorcyclist, bicyclists, motor vehicle). |
| `08_interaction.R` | Explores effect modification by city-level characteristics (e.g., average temperature, peak-hour travel times). |
| `09_temp_cluster.R` | Explores effect modification by temperature clusters. |
| `10_EDF.R` | Generates excess death fractions (EDF) plots for main and subgroup analyses. |
| `11_non-imputation.R` | Runs non-imputed models as a sensitivity check and comparison against imputed results. |
| `11a_sensitivity_outcome.R` | Compares effect estimates between imputed vs. non-imputed outcomes. |
| `12_demo.R` | Simulates a lightweight simulated dataset (10 cities × 30 days × 50 imputations) to demonstrate core pipeline functions. |
---

## 💻 System Requirements

### **Software**
- **R version** ≥ 4.2.0 (tested with 4.3.1)

### **Operating Systems Tested**
- macOS Sequoia 15.2

## ⚙️ Installation Guide
- Clone the repository.
- Install required R packages, e.g., install.packages(c("data.table", "dlnm", "gnm")).
- Update file paths in scripts: from "/Volumes/TOSHIBA Kai/MS252/..." to "your_project_folder"
- Package installation: ~3–5 minutes.

## 🚀 Demo Instructions
- Run 12_demo.R: simulates 10 cities × 30 days × 50 imputations.
- It (1) runs DLNM main effect; (2) applies Rubin’s rule; (3) generates subgroup and interaction plots.
- Expected Output: 3-panel plot including main effect, subgroup (e.g., male deaths), interaction (e.g., city-level mean temperature) results.
- Run Time: ~1 minutes on a standard desktop (8GB RAM).

## 📘 Instructions for Use (Full Pipeline)
### Step 1: Data Preparation
- Run 00_data wrangling_1.R
- Run 01_data_wrangling_2.R

### Step 2: Main Effects 
- 02_median_imputation_main_effects.R
- 03_rubin_imputation_main_effect.R

### Step 3: Sensitivity Analyses
- 04_sensitivity_knots.R
- 11_non-imputation.R
- 11a_sensitivity_outcome.R

### Step 4: Stratification and interaction
- 05_subgroup_sex.R
- 06_subgroup_age.R
- 07_subgroup_mode.R
- 08_interaction.R
- 09_temp_cluster.R

### Step 5: Attributable Risk (EDF)
- 10_EDF.R 


