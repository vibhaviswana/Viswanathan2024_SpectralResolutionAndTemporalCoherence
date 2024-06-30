# Viswanathan2024_SpectralResolutionAndTemporalCoherence

This repository contains the code used in 
[Viswanathan et al. (2024)](https://doi.org/10.1101/2024.03.11.584489).

Viswanathan, V., Heinz, M. G., & Shinn-Cunningham, B. G. (2024). Impact of Reduced 
Spectral Resolution on Temporal-Coherence-Based Source Segregation. 
bioRxiv 2024.03.11.584489. doi: https://doi.org/10.1101/2024.03.11.584489.	

## Basic usage

Before running the code in this repository, download the Viswanathan et al. (2022)
temporal-coherence-based source segregation model available from: 
https://github.com/vibhaviswana/modeling-consonant-confusions.

The code in this repository should be run using the following steps: 
- Add the Viswanathan et al. (2022) model code to the MATLAB path.
- Add the Bruce et al. (2018) auditory-nerve model code to the MATLAB path.
- Run ```predictCMR.m```. This runs the Bruce et al. (2018) auditory-nerve model,  
followed by the Viswanathan et al. (2022) temporal-coherence-based source segregation model 
on the CMR stimuli in ```CMRtestStimuli.mat```.

---
**Note**

The code calls the Bruce et al. (2018) model, which is available as part of the [Auditory Modeling Toolbox; AMT](https://amtoolbox.org).

---
