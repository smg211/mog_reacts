# mog_reacts_2024
This repository contains MATLAB code for detecting and decoding reactivations as in Griffin et al., 2024 (Nature). A complete description of the method is described in the paper.

System Requirements:
  - MATLAB (code was tested on R2020a)
  - ??

There are three main functions included in the main directory of the repository:
- PCAReact. This function implements PCA-based detection of reactivation events, and calculates the reactivation rate for each PC during each break. Optionally, one can also determine the reactivation rate of shuffled data to determine the null distribution of reactivation rates.
- DecodeReach. This function implements Bayesian decoding of reach direction during task execution and performs k-fold cross-validation at various lags to determine the optimal lag for the decoder.
- DecodeReact. This function utilizes the tuning profiles at the specified (i.e., optimal) lag derived in DecodeReach to attempt to decode reach directions from reactivation events. It further calculates the reactivation rate for each direction during each break. Optionally, one can also determine the reach reactivation rate of shuffled data to determine the null distribution of reach reactivation rates.

These functions are combined in the example.m script, which further includes a SimData function for generating simulated kinematics and spiking data to test the pipeline, as well as some additional plots to demonstrate the efficacy of the method.
