# DREAM-Suite
DiffeRential Evolution Adaptive Metropolis algorithm: MATLAB and Python Toolbox

## Description

Bayesian inference has found widespread application and use in science and engineering to reconcile Earth system models with data, including prediction in space (interpolation), prediction in time (forecasting), assimilation of observations and deterministic/stochastic model output, and inference of the model parameters. Bayes theorem states that the posterior probability, $p(H|\tilde{\textbf{Y}})$ of a hypothesis, $H$ is proportional to the product of the prior probability, $p(H)$ of this hypothesis and the likelihood, $L(H|\tilde{\textbf{Y}})$ of the same hypothesis given the new observations, $\tilde{\textbf{Y}}$, or $p(H|\tilde{\textbf{Y}}) \propto p(H) L(H|\tilde{\textbf{Y}})$. In science and engineering, $H$ often constitutes some numerical model, $\mathcal{F}(\textbf{x},\cdot)$ which summarizes, in algebraic and differential equations, state variables and fluxes, all knowledge of the system of interest, and the unknown parameter values, \textbf{x} are subject to inference using the data $\tilde{\textbf{Y}}$. Unfortunately, for complex system models the posterior distribution is often high dimensional and analytically intractable, and sampling methods are required to approximate the target. DREAM-Suite implements many different implementations of the DiffeRential Evolution Adaptive Metropolis (DREAM) algorithm of \cite{vrugt2008a,vrugt2009a} including the DREAM, DREAM_{(LOA)}, DREAM_{(ABC)}, DREAM_{(BMA)}, DREAM_{(D)}, DREAM_{(ZS)}, DREAM_{(ZS)}, MTDREAM_{(ZS)} and DREAM_{(KZS)} algorithms. The toolbox in MATLAB and Python accomodates continuous and/or discrete variables and uses single- and multi-try sampling from an archive of current or past states using parallel direction, snooker and/or Kalman candidate points. DREAM Suite provides scientists and engineers with an arsenal of options and utilities to solve posterior sampling problems involving (among others) bimodality, high-dimensionality, summary statistics, bounded parameter spaces, dynamic simulation models, formal/informal likelihood functions (GLUE), distribution-adaptive likelihood functions, diagnostic model evaluation, data assimilation, Bayesian model averaging, distributed computation, and informative/noninformative prior distributions. DREAM-Suite supports parallel computing and includes tools for convergence analysis of the sampled chain trajectories and post-processing of the results. 36 built-in case studies illustrate the main capabilities and functionalities of DREAM-Suite. 

## Getting Started

### Installing: MATLAB

* Download and unzip the zip file 'MATLAB_code_DREAM_Suite_V2.0.zip' in a directory 'DREAM-Suite'
* Add the toolbox to your MATLAB search path by running the script 'install_DREAM_Suite.m' available in the root directory
* You are ready to run the examples

### Executing program

* After intalling, you can simply direct to each example folder and execute the local 'example_X.m' file
* Please Make sure you read carefully the instructions (i.e., green comments) in 'install_DREAM_Suite.m' and the manual !!!  

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
2. Vrugt, J.A. (2024), Distribution-Based Model Evaluation and Diagnostics: Elicitability, Propriety, and Scoring Rules for Hydrograph Functionals, Water Resources Research, 60, e2023WR036710, https://doi.org/10.1029/2023WR036710
3. Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of distribution-adaptive likelihood functions: Generalized and universal likelihood functions, scoring rules and multi-criteria ranking, Journal of Hydrology, 615, Part B, 2022, https://doi.org/10.1016/j.jhydrol.2022.128542
4. Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software package: Theory, concepts, and MATLAB implementation, Environmental Modeling and Software, 75, pp. 273-316, https://doi.org/10.1016/j.envsoft.2015.08.013
5. Sadegh, M., and J.A. Vrugt (2014), Approximate Bayesian computation using Markov chain Monte Carlo simulation: DREAM_(ABC), Water Resources Research, https://doi.org/10.1002/2014WR015386
6. Vrugt, J.A., and M. Sadegh (2013), Toward diagnostic model calibration and evaluation: Approximate Bayesian computation, Water Resources Research, 49, pp. 4335–4345, https://doi.org/10.1002/wrcr.20354
7. Laloy, E., and J.A. Vrugt (2012), High-dimensional posterior exploration of hydrologic models using multiple-try DREAM_(ZS) and high-performance computing, Water Resources Research, 48, W01526, https://doi.org/10.1029/2011WR010608
8. Vrugt, J.A., and C.J.F. ter Braak (2011), DREAM_(D): An adaptive Markov chain Monte Carlo simulation algorithm to solve discrete, noncontinuous, and combinatorial posterior parameter estimation problems, Hydrology and Earth System Sciences, 15, pp. 3701-3713, https://doi.org/10.5194/hess-15-3701-2011
9. Vrugt, J.A., C.J.F. ter Braak, H.V. Gupta, and B.A. Robinson (2009), Equifinality of formal (DREAM) and informal (GLUE) Bayesian approaches in hydrologic modeling? Stochastic Environmental Research and Risk Assessment, 23(7), pp. 1011-1026, https://doi.org/10.1007/s00477-008-0274-y
10. Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and J.M. Hyman (2009), Accelerating Markov chain Monte Carlo simulation by differential evolution with self-adaptive randomized subspace sampling, International Journal of Nonlinear Sciences and Numerical Simulation, 10(3), pp. 271-288
11. Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and B.A. Robinson (2008), Treatment of input uncertainty in hydrologic modeling: Doing hydrology backward with Markov chain Monte Carlo simulation, Water Resources Research, 44, W00B09, https://doi.org/10.1029/2007WR006720
12. Ter Braak, C.J.F., and J.A. Vrugt (2008), Differential Evolution Markov Chain with snooker updater and fewer chains, Statistics and Computing, https://doi.org/10.1007/s11222-008-9104-9
13. Ter Braak, C.J.F. (2006), A Markov Chain Monte Carlo version of the genetic algorithm differential evolution: easy Bayesian computing for real parameter spaces, Statistics and Computing, 16, pp. 239-249, doi:10.1007/s11222-006-8769-1

## Version History

* 1.0
    * Initial Release
* 2.0
    * Implemented many new capabilities including distribution-adaptive likelihood functions and scoring rules and performance metrics of Bayesian predictive distribution
    * Python implementation
    * Source code in MATLAB and Python

## Acknowledgments
