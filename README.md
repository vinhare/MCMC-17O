# MCMC-17O

# Bayesian Triple Oxygen Isotope Parameter Estimation

## Overview
This repository contains the MATLAB script **TripleOxygenBayes.m**, which implements a Bayesian Markov Chain Monte Carlo (MCMC) approach for estimating parameters related to the triple oxygen isotope composition of biominerals, specifically ratite eggshells. The methodology is based on Bayesian inference to model the isotopic fractionation processes affecting body water composition, providing insights into past CO₂ levels and global primary productivity (GPP).

The script is associated with the paper:

**Hare, V.J., Yarian, D.A., Faith, J.T., Harris, C., Lee-Thorp, J.A., Passey, B.H., Sokolowski, K.G., & Ségalen, L. (2025). Triple oxygen isotopes in eggshell carbonate as a proxy of late Cenozoic CO₂ and primary productivity. *under review*.**

## Features
- Implements Bayesian parameter estimation using the MCMC toolbox for MATLAB.
- Incorporates external functions for body water isotope modeling.
- Supports analysis of fossil and modern eggshell carbonate isotopes.
- Allows estimation of key physiological and environmental parameters.

## Dependencies
This script requires the following MATLAB toolboxes and external functions:

### Required MATLAB Toolboxes:
- **MCMC toolbox** by Marko Laine: [https://mjlaine.github.io/mcmcstat/](https://mjlaine.github.io/mcmcstat/)
- **DERIVESTsuite** for numerical differentiation.

### External Functions (included in the repository or must be provided by the user):
- `BWM_Emu_st.m`
- `TripleOxygenss.m`
- `Body_Water_Model_rh.m`
- `SS_over_rh.m`

## Input Data
The script requires a dataset in MATLAB `.mat` format with the following structure:
- `data.ydata`: Measured **Δ'¹⁷O_bw** values (per meg)
- `data.xdata`: Measured **δ'¹⁸O_bw** values (per mil)
- `data.yerr`: Uncertainty in **Δ'¹⁷O_bw** (1σ standard deviation)

Example dataset: `TripleO_Central_Modern.mat`

## Parameter Definitions
The model parameters estimated using MCMC are:
- **θ₁ (DeltaOatm)**: Atmospheric Δ'¹⁷O anomaly (per meg)
- **θ₂ (deltaOmw)**: δ¹⁸O of meteoric water (per mil)
- **θ₃ (rh)**: Relative humidity (fraction, 0-1)
- **θ₄ (Leafratio)**: Ratio of leaf water uptake
- **θ₅ (WEIndex)**: Water Economy Index (mL/kJ)
- **θ₆ (d¹⁸O_atmO₂)**: δ¹⁸O of atmospheric O₂ (per mil)

## Running the Script
1. Load MATLAB and navigate to the script directory.
2. Ensure that the required `.mat` dataset is available in the working directory.
3. Run the script by executing:
   ```matlab
   TripleOxygenBayes
4. The script will perform MCMC sampling and generate posterior distributions of the model parameters.

## Output
The script produces:
- **MCMC chains**: Stored in `chain` variable.
- **Plots**:
  - Trace plots of MCMC chains.
  - Posterior density estimates.
  - Parameter correlation scatter plots.
- **Final estimated atmospheric Δ'¹⁷O and uncertainty**, computed using bootstrapping.

For the Central Namib (Modern) dataset, and example of the Parameter correlation scatter plots is shown below:
  
<img src="https://github.com/user-attachments/assets/8263a360-64fb-4ed7-946b-94ffa32c20fc" width="600">

With the posterior density estimates:

<img src="https://github.com/user-attachments/assets/c846a0d1-d52c-4991-8d1f-f146f255fddf" width="600">

And finally, the "best fit" model, using the posterior parameters.

<img src="https://github.com/user-attachments/assets/b712231d-01b4-41a7-8e9f-8b7aacab52c5" width="600">

## Interpretation
The results can be used to:
- Infer past **CO₂** levels and **GPP** using Bayesian inverse modeling.
- Understand physiological and environmental factors influencing isotope ratios.
- Compare fossil eggshell isotope signatures with modern analogues.

## IMPORTANT NOTE!
The code here is optimised for ratites. If you are going to apply this to other animals, you will need to spend some time properly setting up the priors for that particular species, including (but not limited to) different body mass, and potentially also different priors for δ18Omw, WEI, rl/s, etc. To help with convergence, you may also find it useful to modify the range of RH values to a tigher range, informed by your study area. For animals that are taking in more meteoric water, you may also find it more useful to solve for f_H20-in-food, rather than rl/s, by modifying the scripts appropriately. In my experience, one can solve for either, but not both together (as they are strongly covariant). 
 
## Citation
If you use this script in your research, please cite:

Hare, V.J., Yarian, D.A., Faith, J.T., Harris, C., Lee-Thorp, J.A., Passey, B.H., Sokolowski, K.G., & Ségalen, L. (2025). Triple oxygen isotopes in eggshell carbonate as a proxy of late Cenozoic CO₂ and primary productivity. *Under review*.

## License
This code is provided under the MIT License. See `LICENSE` for details.

---

For any questions, please contact: **vincent.john.hare@gmail.com**

