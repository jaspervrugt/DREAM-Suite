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
%% Example 32: Distribution-adaptive likelihood functions applied to conceptual       %%
%%             watershed models using measured discharge data                         %%
%%                                                                                    %%
%% Check the following papers                                                         %%
%%   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of %%
%%       distribution-adaptive likelihood functions: Generalized and universal        %%
%%       likelihood functions, scoring rules and multi-criteria ranking, Journal of   %%
%%       Hydrology, 615, Part B, 2022, doi:10.1016/j.jhydrol.2022.128542.             %%
%%       https://www.sciencedirect.com/science/article/pii/S002216942201112X          %%
%%   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software %%
%%       package: Theory, concepts, and MATLAB implementation, Environmental Modeling %%
%%       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             %%
%%   Vrugt, J.A., C.G.H. Diks, H.V. Gupta, W. Bouten, and J.M. Verstraten (2005),     %%
%%       Improved treatment of uncertainty in hydrologic modeling: Combining the      %%
%%       strengths of global optimization and data assimilation, Water Resources      %%
%%       Research, 41, W01017, doi:10.1029/2004WR003059                               %%
%%   Schoups, G., J.A. Vrugt, F. Fenicia, and N.C. van de Giesen (2010), Corruption   %%
%%       of accuracy and efficiency of Markov Chain Monte Carlo simulation by         %%
%%       inaccurate numerical implementation of conceptual hydrologic models, Water   %%
%%       Resources Research, 46, W10530, doi:10.1029/2009WR008648                     %%
%%   Schoups, G., and J.A. Vrugt (2010), A formal likelihood function for parameter   %%
%%       and predictive inference of hydrologic models with correlated,               %%
%%       heteroscedastic and non-Gaussian errors, Water Resources Research, 46,       %%
%%       W10531, doi:10.1029/2009WR008933                                             %%
%%   Vrugt, J.A., C.J.F. ter Braak, H.V. Gupta, and B.A. Robinson (2009),             %%
%%       Equifinality of formal (DREAM) and informal (GLUE) Bayesian approaches in    %%
%%       hydrologic modeling?, Stochastic Environmental Research and Risk Assessment, %%
%%       23(7), 1011-1026, doi:10.1007/s00477-008-0274-y                              %%
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
DREAMPar.lik = 17;      % 13: Normal likelihood function with AR(2)-model of residuals
                        % 14: Generalized likelihood function (Obsolete)
                        % 16: Laplacian likelihood function with AR(2)-model of residuals
                        % 17: Skewed Student t likelihood function
                        % 44: Generalized likelihood function PLUS
                        % 45: Universal likelihood function
%% Please refer to Vrugt et al. (2022) for theory on likelihood functions

%% Provide information parameter space and initial sampling
Par_info.initial = 'latin';                 % Latin hypercube sampling
Par_info.boundhandling = 'reflect';         % Explicit boundary handling
Par_info.norm = 1;                          % Work in normalized space

%% Define name of function (.m file) for posterior exploration
Func_name = 'hymod';    %1 C++ implementation

%% Load the data
[Meas_info,plugin] = load_data;

%% Define names of parameters
par_names = {'S_{\rm u,max}','\beta','\alpha','K_{\rm s}','K_{\rm f}'};
% index:         1     2     3      4     5
% parname:     Sumax  beta  alfa    Ks    Kf
fpar_mod   =  [ 200    1     0.5   0.01  0.50 ];
parmin_mod =  [  50    0.10   0    1e-5  1e-1 ];
parmax_mod =  [ 1000   10.0   1    1e-1   5   ];

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_zs';
switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;            % # Markov chains
        DREAMPar.T = 3000;          % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 4;             % # Markov chains
        DREAMPar.T = 30000;         % # generations
end

%% Optional DREAM settings
options.parallel = 'yes';           % This example runs in parallel
options.save = 'yes';               % Save DREAM workspace during run
options.modout = 'no';              % Store model simulations

%% Define likelihood properties: [user needs to select idx_nuis]
[DREAMPar,Par_info,id_vpar,fpar,fname] = define_lik(DREAMPar,Par_info,...
    par_names,fpar_mod,parmin_mod,parmax_mod);
DREAMPar.d = size(Par_info.min,2);  % # parameters

%% Define global likelihood variables
global LV; LV.id_vpar = id_vpar; LV.fpar = fpar; LV.filename = fname;

%% Run DREAM-Suite
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar,...
    Par_info,Meas_info,options,[],plugin);
