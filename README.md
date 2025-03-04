# DREAM-Suite
DiffeRential Evolution Adaptive Metropolis algorithm: MATLAB and Python Toolbox

## Description

Bayesian inference has found widespread application and use in science and engineering to reconcile Earth system models with data, including prediction in space (interpolation), prediction in time (forecasting), assimilation of observations and deterministic/stochastic model output, and inference of the model parameters. Bayes theorem states that the posterior probability, $p(H\vert\tilde{\textbf{y}})$ of a hypothesis, $H$ is proportional to the product of the prior probability, $p(H)$ of this hypothesis and the likelihood, $L(H \vert \tilde{\textbf{y}})$ of the same hypothesis given the new observations, $\tilde{\textbf{y}}$, or $p(H\vert\tilde{\textbf{y}}) \propto p(H) L(H\vert\tilde{\textbf{y}})$. In science and engineering, $H$ often constitutes a mathematical model, $y \leftarrow F(\textbf{x},\cdot)$ which summarizes, in algebraic and differential equations, state variables and fluxes, all knowledge of the system of interest, and the unknown parameter values, $\textbf{x}$, are subject to inference using the training data $\tilde{\textbf{y}}$. Unfortunately, for complex system models the posterior distribution is often high dimensional and analytically intractable, and sampling methods are required to approximate the target. DREAM-Suite is a MATLAB and Python toolbox of the DiffeRential Evolution Adaptive Metropolis (DREAM) algorithm of Vrugt et al. 2008a, 2009a and implements in one package the DREAM$_{(LOA)}$ (Vrugt and Beven, 2018), DREAM$_{(ABC)}$ (Sadegh and Vrugt, 2014), DREAM$_{(BMA)}$ (Vrugt et al., 2008b), DREAM$_{(D)}$ (Vrugt and ter Braak, 2011), DREAM$_{(ZS)}$, MTDREAM$_{(ZS)}$ (Laloy and Vrugt, 2012) and DREAM$_{(KZS)}$ (Zhang et al., 2020) algorithms. These multi-chain algorithms handle high-dimensionality and multi-modality and evolve the chains to a limiting distribution via single- or multi-try sampling from an archive of current or past states using parallel direction, snooker and/or Kalman candidate points. DREAM-Suite provides scientists and engineers with an arsenal of options and capabilities for solving posterior estimation problems involving bounded and unbounded parameter spaces, continuous and discrete variables, formal & distribution-adaptive likelihood functions, informal & pseudo-likelihood functions, informative/noninformative univariate and multivariate prior distributions and summary statistics in the context of memory-free and dynamic simulation models, diagnostic model evaluation, data assimilation, Bayesian model averaging and diagnostic Bayes. DREAM-Suite implements multi-core (thread) evaluation of the chains and/or candidate points. The postprocessor automatically generates a suite of figures including traceplots of the univariate and multivariate scale reduction factors, autocorrelation and traceplots of the sampled parameters, prior and likelihood, histograms of the marginal parameter distributions, matrix plot of the posterior samples, confidence and prediction intervals and Tables with scoring rules and performance metrics of the Bayes predictive distribution. 36 built-in case studies illustrate the main capabilities and functionalities of DREAM-Suite. 

## Getting Started

### Installing: MATLAB

* Download and unzip the zip file 'MATLAB_code_DREAM_Suite_V2.0.zip' in a directory 'DREAM-Suite'
* Add the toolbox to your MATLAB search path by running the script 'install_DREAM_Suite.m' available in the root directory
* You are ready to run the examples

### Executing program

* After intalling, you can simply direct to each example folder and execute the local 'example_X.m' file
* Please make sure you read carefully the instructions (i.e., green comments) in 'install_DREAM_Suite.m' and the manual !!!  

### Installing: Python

* Download and unzip the zip file 'Python_code_DREAM_Suite_V2.0.zip' to a directory called 'DREAM-Suite'

### Executing program

* Go to Command Prompt and directory of example_X in the root of DREAM-Suite
* Now you can execute this example by typing 'python example_X.py'.
* Instructions can be found in the file 'DREAM_Suite.py' and in the manual !!!  

## Authors

* Vrugt, Jasper A. (jasper@uci.edu)

## Literature
1. Vrugt, J.A., R. de Punder, and P. Grünwald, A sandwich with water: Bayesian/Frequentist uncertainty quantification under model misspecification, Submitted to Water Resources Research, May 2024, https://essopenarchive.org/users/597576/articles/937008-a-sandwich-with-water-bayesian-frequentist-uncertainty-quantification-under-model-misspecification
2. Vrugt, J.A. (2024), Distribution-Based Model Evaluation and Diagnostics: Elicitability, Propriety, and Scoring Rules for Hydrograph Functionals, _Water Resources Research_, 60, e2023WR036710, https://doi.org/10.1029/2023WR036710
3. Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of distribution-adaptive likelihood functions: Generalized and universal likelihood functions, scoring rules and multi-criteria ranking, _Journal of Hydrology_, 615, Part B, 2022, https://doi.org/10.1016/j.jhydrol.2022.128542
4. Zhang, J., J.A. Vrugt, X. Shi, G. Lin, L. Wu, and L. Zeng (2020), Improving simulation efficiency of MCMC for inverse modeling of hydrologic systems with a Kalman-inspired proposal distribution, _Water Resources Research_, 56, e2019WR025474. https://doi.org/10.1029/2019WR025474
5. Vrugt, J.A., and K.J. Beven (2018), Embracing equifinality with efficiency: Limits of Acceptability sampling using the DREAM(LOA) algorithm, _Journal of Hydrology_, 559, pp. 954-971, https://doi.org/10.1016/j.jhydrol.2018.02.026.
6. Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software package: Theory, concepts, and MATLAB implementation, _Environmental Modeling and Software_, 75, pp. 273-316, https://doi.org/10.1016/j.envsoft.2015.08.013
7. Sadegh, M., and J.A. Vrugt (2014), Approximate Bayesian computation using Markov chain Monte Carlo simulation: DREAM_(ABC), _Water Resources Research_, https://doi.org/10.1002/2014WR015386
8. Vrugt, J.A., and M. Sadegh (2013), Toward diagnostic model calibration and evaluation: Approximate Bayesian computation, _Water Resources Research_, 49, pp. 4335–4345, https://doi.org/10.1002/wrcr.20354
9. Laloy, E., and J.A. Vrugt (2012), High-dimensional posterior exploration of hydrologic models using multiple-try DREAM_(ZS) and high-performance computing, _Water Resources Research_, 48, W01526, https://doi.org/10.1029/2011WR010608
10. Vrugt, J.A., and C.J.F. ter Braak (2011), DREAM_(D): An adaptive Markov chain Monte Carlo simulation algorithm to solve discrete, noncontinuous, and combinatorial posterior parameter estimation problems, _Hydrology and Earth System Sciences_, 15, pp. 3701-3713, https://doi.org/10.5194/hess-15-3701-2011
11. Vrugt, J.A., C.J.F. ter Braak, H.V. Gupta, and B.A. Robinson (2009), Equifinality of formal (DREAM) and informal (GLUE) Bayesian approaches in hydrologic modeling? _Stochastic Environmental Research and Risk Assessment_, 23(7), pp. 1011-1026, https://doi.org/10.1007/s00477-008-0274-y
12. Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and J.M. Hyman (2009a), Accelerating Markov chain Monte Carlo simulation by differential evolution with self-adaptive randomized subspace sampling, _International Journal of Nonlinear Sciences and Numerical Simulation_, 10(3), pp. 271-288
13. Vrugt, J.A., C.G.H. Diks and M.P. Clark (2008b), Ensemble Bayesian model averaging using Markov Chain Monte Carlo sampling, _Environmental Fluid Mechanics_, 8 (5), pp. 579-595, https://doi.org/10.1007/s10652-008-9106-3
14. Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and B.A. Robinson (2008a), Treatment of input uncertainty in hydrologic modeling: Doing hydrology backward with Markov chain Monte Carlo simulation, _Water Resources Research_, 44, W00B09, https://doi.org/10.1029/2007WR006720
15. Ter Braak, C.J.F., and J.A. Vrugt (2008), Differential Evolution Markov Chain with snooker updater and fewer chains, _Statistics and Computing_, https://doi.org/10.1007/s11222-008-9104-9
16. Ter Braak, C.J.F. (2006), A Markov Chain Monte Carlo version of the genetic algorithm differential evolution: easy Bayesian computing for real parameter spaces, _Statistics and Computing_, 16, pp. 239-249, doi:10.1007/s11222-006-8769-1
17. Vrugt, J.A., C.G.H. Diks, H.V. Gupta, W. Bouten, and J.M. Verstraten (2005), Improved treatment of uncertainty in hydrologic modeling: Combining the strengths of global optimization and data assimilation, _Water Resources Research_, 41, W01017, https://doi.org/10.1029/2004WR003059.
18. Vrugt, J.A., S.C. Dekker, and W. Bouten (2003), Identification of rainfall interception model parameters from measurements of throughfall and forest canopy storage, Water Resources Research, 39(9), 1251, https://doi.org/10.1029/2003WR002013.

## Version History

* 1.0
    * Initial Release
* 2.0
    * Implemented many new capabilities including distribution-adaptive likelihood functions and scoring rules and performance metrics of Bayesian predictive distribution
    * Python implementation
    * Source code in MATLAB and Python

## Built-in case studies

1. example 1: $d$-dimensional banana shaped Gaussian distribution
2. example 2: $d$-dimensional Gaussian distribution
3. example 3: $d$-dimensional multimodal normal mixture distribution
4. example 4: real-world example rainfall-runoff (hymod in C++/MATLAB/Python)
5. example 5: rainfall-runoff (hymod as external executable)
6. example 6: hmodel with distribution-adaptive likelihood functions (Vrugt et al., 2022)
7. example 7: HYDRUS-1D soil hydraulic model: multiplicative prior
8. example 8: Approximate Bayesian Computation: Benchmark function
9. example 9: Spectral likelihood function in watershed modeling
10. example 10: Gaussian mixture distibution: multivariate prior
11. example 11: $d$-variate t-distribution: df degrees of freedom & correlation matrix R
12. example 12: pedometrics problem involving variogram fitting
13. example 13: Nash-Cascade hydrograph
14. example 14: Approximate Bayesian Computation watershed signatures
15. example 15: Approximate Bayesian Computation bivariate normal benchmark test
16. example 16: Hydrogeophysical inversion
17. example 17: Watershed model, normal, AR(1) and heteroscedastic likelihood
18. example 18: Lotka-Volterra model: informal likelihood (GLUE)
19. example 19: Bayesian Model Averaging [I recommend MODELAVG toolbox]
20. example 20: Limits of acceptability: Soil temperature modeling
21. example 21: Limits of acceptability: Soil moisture model HYDRUS-1D
22. example 22: Limits of acceptability: Nash-Cascade hydrograph
23. example 23: Limits of acceptability: SAC-SMA (old C-code Euler int.)
24. example 24: Flow duration curve fitting
25. example 25: Bedrock depth from high-res topo data & geomorph model
26. example 26: Data assimilation Lorenz-1963 model (SODA: Vrugt et al. 2005)
27. example 27: Data assimilation interception model (Vrugt et al. 2003)
28. example 28: Rainfall and hmodel parameter estimation from streamflow
29. example 29: Gaussian mixture distribution & multiplicative prior
30. example 30: Predator prey interactions
31. example 31: AR(2)-parameter estimation: Test of distribution-adaptive likelihood functions
32. example 32: Distribution-adaptive likelihood functions
33. example 33: 2-dimensional rectangular target distribution
34. example 34: Haverkamp infiltration equation using HYDRUS-1D data
35. example 35: Haverkamp infiltration equation using SWIG database
36. example 36: Sediment transport modeling

## Acknowledgments
