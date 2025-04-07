Code for running the transient analysis of seagrass-soil model
==============

***Antoine Le Vilain

<ins>Contact</ins>: antoine.levilain18@gmail.com

This repository contains the code used to generate the figures and analyses supporting: Temporal Lags in the Dynamics of Carbon Storage Following the Collapse of Seagrass Ecosystems Depend on Type of Stressor. By running the included R and Python scripts, you can replicate the transient dynamics simulations of a seagrass–soil ecosystem under gradual changes in environmental stressors. The files are organised so that each figure (or set of figures) has associated scripts for both data generation (solving and storing model outputs) and plotting the final results.

Below is an overview of the key scripts and their roles in creating each figure:

figure1_equilibrium.R & figure1_plot.R
Identify equilibrium solutions of the seagrass–soil model along different stressor gradients and produce the bifurcation diagrams and transient plots (Figure 1).

figure2_plot.R
Simulate and visualise how seagrass biomass (S) and stable carbon (CB) transients respond over time to fast vs. slow rates of environmental change (Figure 2).

figure3_lag.R & figure3_plot.R
Investigate the time lag between seagrass collapse and subsequent carbon loss across a range of stressor-increase rates (Figure 3).

figure4_pattern_lag_2D.R & figure4_plot.py
Explore pairwise stressor scenarios and record lag and patterns of S and CB changes (Figure 4).

figure5_lag_3D.R & figure5_plot.py
Extend the lag analysis to a 3D setting where all stressors change together (Figure 5).

How to Reproduce the Results
Clone this repository to your local machine or download the zip file.

Install dependencies:

The R scripts rely on packages such as deSolve, rootSolve, dplyr, nleqslv, ggplot2, doMPI, and others (see the package loading sections near the top of each .R file).

The Python scripts typically require matplotlib and pandas (and possibly others). Check the imports at the top of each .py file.

Set your working directory in each R or Python script as needed.

Run the scripts in the order corresponding to the figure you want to reproduce. For instance, to replicate Figure 1 results, first run figure1_equilibrium.R (to generate or confirm equilibrium data) and then run figure1_plot.R to create the final plot.

We performed our analysis using R version 4.1.1 and Python 3.7.3, but more recent versions should generally be compatible. If you encounter any issues with package installations or script errors, please reach out via email.
