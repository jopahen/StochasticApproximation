# Stochastic Approximation in Mathematical Finance
This Git houses the resources of my end of year project for the course Stochastic Simulations (EPFL, autumn 2021).

# Overview of Files
The /Documentation folder has project description and a copy of my report as PDF.

The /legacy folder is old code, not used in the final project. It is not well documented/commented and can be ignored.

The /Plots folder contains the .eps files used for the plots in the report.

All R-scripts in the main folder are related to the implementations & experiments carried out with respect to the project description and as recorded in the project report (see /Documentation). In short:

BS_functions.r defines auxiliary functions for all other scripts

BS_IV_rootfind.r uses Brent algorithm to determine implied volatility

IS_pricing.r contains implementation and analysis of MC-pricers using importance sampling with modified drift

optimal_drift.r contains a routine to compute optimal importance sampling drifts for different volatilities

RM_IV_Put_xxx.r are implementations of RM-algorithms with different modifications according to the experiments carried out during the project

RM_xxx_Analysis.r contain the respective performance analyses conducted in the project
