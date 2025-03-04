%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                                                                                    %%
%% DDDDD   RRRR    EEEEE    AA    MM   MM      SSSSSS  UU   UU   II   TTTTTTTT  EEEEE %%
%% DDDDDD  RRRR    EEEEE   AAAA   MM   MM      SSSSS   UU   UU   II   TTTTTTTT  EEEEE %%
%% DD  DD  RR RR   EE     AA  AA  MMM MMM      SS      UU   UU   II      TT     EE    %%
%% DD  DD  RR RR   EEE    AA  AA  MMMMMMM ---- SS      UU   UU   II      TT     EEE   %%
%% DD  DD  RRRRR   EEE    AAAAAA  MMM MMM ---- SSSSSS  UU   UU   II      TT     EEE   %%
%% DD  DD  RR RR   EE     AAAAAA  MM   MM          SS  UU   UU   II      TT     EE    %%
%% DDDDDD  RR  RR  EEEEE  AA  AA  MM   MM       SSSSS  UUUUUUU   II      TT     EEEEE %%
%% DDDDD   RR  RR  EEEEE  AA  AA  MM   MM      SSSSSS  UUUUUUU   II      TT     EEEEE %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                                                                                    %%
%% Example 5:  Gaussian likelihood applied to conceptual watershed model coded in     %%
%%             Fortran                                                                %%
%%                                                                                    %%
%% Check the following papers                                                         %%
%%   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software %%
%%       package: Theory, concepts, and MATLAB implementation, Environmental Modeling %%
%%       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             %%
%%   Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and        %%
%%       J.M. Hyman (2009), Accelerating Markov chain Monte Carlo simulation by       %%
%%       differential evolution with self-adaptive randomized subspace sampling,      %%
%%       International Journal of Nonlinear Sciences and Numerical Simulation, 10(3), %%
%%       271-288                                                                      %%
%%   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and B.A. Robinson (2008), %%
%%       Treatment of input uncertainty in hydrologic modeling: Doing hydrology       %%
%%       backward with Markov chain Monte Carlo simulation, Water Resources Research, %%
%%       44, W00B09, doi:10.1029/2007WR006720                                         %%
%%   Vrugt, J.A., H.V. Gupta, W. Bouten and S. Sorooshian (2003), A Shuffled Complex  %%
%%       Evolution Metropolis algorithm for optimization and uncertainty assessment   %%
%%       of hydrologic model parameters, Water Resour. Res., 39 (8), 1201,            %%
%%       doi:10.1029/2002WR001642.                                                    %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% Problem settings defined by user
DREAMPar.d = 5;                          % Dimension of the problem
DREAMPar.lik = 11;                       % Model output is simulation: Gaussian likelihood function

%% Provide information parameter space and initial sampling
Par_info.initial = 'latin';                 % Latin hypercube sampling
Par_info.boundhandling = 'reflect';         % Explicit boundary handling
Par_info.min = [1.0 0.10 0.10 0.00 0.10];   % If 'latin', min values
Par_info.max = [500 2.00 0.99 0.10 0.99];   % If 'latin', max values

%% Define name of function (.m file) for posterior exploration
Func_name = 'hymodFORTRAN';

%% Define the measured streamflow data
load bound.txt; Meas_info.Y = bound(65:795,4);

%% Optional settings
options.parallel = 'yes';               % Run each chain on a different core
options.IO = 'yes';                     % Input-output writing of model files (only for parallel!)
options.modout = 'yes';                 % Return model (function) simulations of samples (yes/no)?

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_d';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;                            % # Markov chains
        DREAMPar.T = 1500;                          % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 3;                             % # Markov chains
        DREAMPar.T = 5000;                          % # generations
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.steps = 1000*ones(1,DREAMPar.d);       % # discrete steps
end

%% Run DREAM-Suite
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar,...
    Par_info,Meas_info,options);
