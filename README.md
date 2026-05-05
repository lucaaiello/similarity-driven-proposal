# Reproducibility Repository

This repository contains the code required to reproduce the results presented in the manuscript.

## 📂 Repository Structure

The repository is organized into two main components:

- **Real Data Analysis**
- **Simulation Studies**

---

## 📊 Real Data Analysis

The folder `Real_Data_Analysis/` contains all materials needed to reproduce the empirical results.

### Contents

- `data/`  
  Input data used for the real data analysis.

- `rda.R`  
  Main script to reproduce the real data analysis results.

  Running this script produces:
  - **Main manuscript (Section 5):**
    - Figure 4  
    - Table 1  
  - **Supplementary material (Section S3):**
    - Figures S6–S8  
    - Table S1  

- `DM_RJMCMC_VS_sampler.cpp`  
  C++ implementation of the sampler used for variable selection in the **Dirichlet–Multinomial regression model**.

---

## 🧪 Simulation Studies

The folder `Simulations/` contains all code required to reproduce the simulation experiments presented in Section 3 of the manuscript.

### Samplers

- `sampler_F_stat.cpp`  
  C++ implementation of the sampler for variable selection based on **F-statistics**.

- `sampler_LR.cpp`  
  C++ implementation of the sampler for variable selection in **linear regression settings**.

---

### R Scripts

- `simulations_F_stat.R`  
  Generates results for the F-statistic simulation study:
  - **Main manuscript:**
    - Figures 1–3  
  - **Supplementary material:**
    - Figure S1  

- `simulations_LR.R`  
  Generates results for the linear regression simulation study:
  - **Main manuscript:**
    - Figure 3  
  - **Supplementary material:**
    - Figures S2–S4  

---

## ⚙️ Requirements

- **R** (≥ 4.x recommended)
- Standard R packages (depending on your scripts)
- C++ compiler (e.g., via `Rcpp` or system toolchain)

---

## ▶️ How to Reproduce Results

### Real Data Analysis

```r
source("Real_Data_Analysis/rda.R")
