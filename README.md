
# ClonalHaematopoiesis

This repository contains the R code required to reproduce the analyses presented in the article:

**"Accounting for sequencing depth when inferring clonal haematopoiesis dynamics from variant allele recordings under a hierarchical mixed-effects logistic model".**

## Data Preprocessing

The raw data are preprocessed using the R script `CH-make-dataset.R`.

The script generates a CSV file containing the processed dataset, which is subsequently used as input for the data analysis.

## Data Analysis

The data analysis is performed in two steps.

### Step 1. Model-specific analysis

The data are analysed separately using each of the 16 candidate models (Models 1–16) through the R script `CH-data-main.R`.

The script takes the parameter `nModel` as input, specifying the candidate model to be fitted.

For each model, the script automatically creates dedicated output folders in the current working directory and exports the corresponding model-specific results. These results are subsequently used for model averaging (Step 2).

### Step 2. Model averaging

The model-specific results obtained in Step 1 are aggregated using the R script `CH-data-mAVG-main.R`.

The script automatically creates the required output folders in the current working directory and exports the aggregated results.

## Simulation Studies

The simulation studies are performed in two steps.

### Step 1. Simulation-specific analysis

The R script `CH-sim-main.R` runs an individual simulation study and exports the corresponding results to automatically generated folders in the current working directory.

### Step 2. Aggregation of simulation results

The simulation-specific results are subsequently combined using the R script `CH-sim-combine.R`.

This script aggregates the results across simulations to facilitate their visual inspection.

## Supporting Functions

The R script `CH-fun.R` contains the supporting functions required for both the simulation studies and the data analyses.

These functions are automatically loaded by the corresponding main scripts.

## Requirements

The analyses were performed using **R version 4.2.1 (2022-06-23)**.

The required R packages are specified in the respective R scripts.
