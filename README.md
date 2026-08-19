# Semi-Markov Surface Residence and Memory Transmission in Confined Bacterial Chemotaxis

This repository contains the MATLAB codes and numerical data used in the manuscript:

> **Semi-Markov Surface Residence and Memory Transmission in Confined Bacterial Chemotaxis**

## Repository Contents

### 1. Calculation Scripts

Scripts beginning with `BS_` and `Parallel_` perform the numerical calculations and data processing used in the study.

- `BS_*.m` — Calculation scripts for the corresponding analyses.
- `Parallel_*.m` — Parallelized calculation scripts for computationally intensive analyses.

### 2. Figure Scripts

Scripts beginning with `fig` generate the corresponding figures in the manuscript.

- `fig2.m` — Code for generating **Figure 2**.
- `fig3.m` — Code for generating **Figure 3**.
- `fig4.m` — Code for generating **Figure 4**.
- `fig5.m` — Code for generating **Figure 5**.

Thus, `figXX.m` corresponds to **Figure XX** in the manuscript.

### 3. Auxiliary Functions

Scripts beginning with `func_` contain auxiliary MATLAB functions used by the calculation and figure-generation scripts.

- `func_*.m` — Auxiliary MATLAB functions.

### 4. Data Files

The `.mat` files contain numerical data generated from simulations and calculations. These data are used as inputs for subsequent analysis and figure generation.

## Requirements

- MATLAB
- Required MATLAB toolboxes depend on the specific scripts being used.

## Usage

Set the repository directory as the MATLAB working directory or add it to the MATLAB path.

To generate the corresponding manuscript figures, run:

```matlab
fig2
fig3
fig4
fig5
