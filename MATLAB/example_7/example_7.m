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
%% Example 7:  Gaussian likelihood applied to soil water flow model coded in Fortran  %%
%%                                                                                    %%
%% Check the following papers                                                         %%
%%   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software %%
%%       package: Theory, concepts, and MATLAB implementation, Environmental Modeling %%
%%       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             %%
%%   Scharnagl, B., J.A. Vrugt, H. Vereecken, and M. Herbst (2011), Bayesian inverse  %%
%%       modeling of soil water dynamics at the field scale: using prior information  %%
%%       on soil hydraulic properties, Hydrology and Earth System Sciences, 15,       %%
%%       3043–3059, doi:10.5194/hess-15-3043-2011.                                    %%
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
DREAMPar.d = 7;                         % Dimension of the problem
DREAMPar.lik = 11;                      % Model output is simulation (soil moisture values)

%% Provide information parameter space and initial sampling (use of informative multiplicative prior)
Par_info.initial = 'prior';             % Initial sample from prior distribution
Par_info.prior = { ...
    'normpdf(0.0670,0.0060)',...        % Marginal prior of theta_r
    'normpdf(0.4450,0.0090)',...        % Marginal prior of theta_s
    'normpdf(-2.310,0.0600)',...        % Marginal prior of log10(alpha)
    'normpdf(0.2230,0.0110)',...        % Marginal prior of log10(n)
    'normpdf(-1.160,0.2700)',...        % Marginal prior of log10(Ks)
    'normpdf(0.3900,1.4700)',...        % Marginal prior of L
    'unifpdf(-250,-50)' };              % Marginal prior of h_bot
Par_info.boundhandling = 'reflect';     % Enforce boundaries of parameter space

%% Define feasible parameter space (minimum and maximum values)
% parname		[ thetar thetas	log10(alpha) log10(n) log10(Ks)	  L    h_bot
Par_info.min =	[ 0.0430 0.4090	 -2.5528	 0.1790   -2.2366	-5.49	-250 ];  % For boundary handling
Par_info.max =	[ 0.0910 0.4810	 -2.0706	 0.2670   -0.0800	 6.27	-50  ];  % For boundary handling

%% Define name of function (.m file) for posterior exploration
Func_name = 'HYDRUS';

%% Define measurement data ( observed soil moisture contents )
data_hydrus = Load_data; Meas_info.Y = data_hydrus.water;

%% Optional settings
options.parallel = 'yes';              % Run chains in parallel
options.IO = 'yes';                    % Input-output writing of model files (for parallel setting only!)
options.save = 'yes';                  % Save memory of DREAM during trial

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs'}
method = 'dream';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;                            % Number of Markov chains
        DREAMPar.T = 2000;                          % Number of generations
    case {'dream_zs','dream_dzs'}
        DREAMPar.N = 3;                             % Number of Markov chains        
        DREAMPar.T = 5000;                          % Number of generations
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.steps = 200*ones(1,DREAMPar.d);        % Number of discrete steps
end

%% Run DREAM-Suite 
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar,...
    Par_info,Meas_info,options,[],data_hydrus);