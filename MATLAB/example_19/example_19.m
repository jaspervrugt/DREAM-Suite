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
%% Example 19: BMA mixture model estimation. I recommend using the MODELAVG toolbox!  %%
%%             This provides a much more comprehensive treatment of the BMA model,    %%
%%             including confidence/prediction intervals, scoring rules and other     %%
%%             metrics of the BMA distribution forecast                               %%
%%                                                                                    %%
%% Check the following papers                                                         %%
%%   Vrugt, J.A. (2024), Distribution-Based Model Evaluation and Diagnostics:         %%
%%       Elicitability, Propriety, and Scoring Rules for Hydrograph Functionals,      %%
%%       Water Resources Research, 60, e2023WR036710,                                 %%
%%       https://doi.org/10.1029/2023WR036710                                         %%
%%   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software %%
%%       package: Theory, concepts, and MATLAB implementation, Environmental Modeling %%
%%       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             %%
%%   Vrugt, J.A., C.G.H. Diks, and M.P. Clark (2008), Ensemble Bayesian model         %%
%%       averaging using Markov chain Monte Carlo sampling, Environmental Fluid       %%
%%       Mechanics, 8(5-6), 579-595, doi:10.1007/s10652-008-9106-3                    %%
%%   Vrugt, J.A., and B.A. Robinson (2007), Treatment of uncertainty using ensemble   %%
%%       methods: Comparison of sequential data assimilation and Bayesian model       %%
%%       averaging, Water Resources Research, 43, W01411, doi:10.1029/2005WR004838    %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% Problem settings defined by user
DREAMPar.lik = 2;       % Model returns log-likelihood

%% Initial sampling and parameter ranges
Par_info.initial = 'latin';         % Latin hypercube sampling
Par_info.boundhandling = 'reflect'; % Explicit boundary handling

%% Define name of function (.m file) for posterior exploration
Func_name = 'BMA_lik';

%% Load data from Vrugt and Robinson, WRR, 43, W01411, doi:10.1029/2005WR004838, 2007
data = load('discharge.txt');           % Daily discharge forecasts (mm/d) 
                                        % verifying data 
T_idx = 1:1:3000;                       % Start/end day of training period
D = data(T_idx,1:8);                    % Ensemble forecasts
y = data(T_idx,9);                      % Verifying observations

options.BMA = 'yes';                    % Activate BMA method
PDF = 'normal';                         % Forecast pdf ('normal'/'gamma')
VAR = '4';                              % Variance option ('1'/'2'/'3'/'4')
[DREAMPar,Par_info,D_bc,A,B] = ...      % Setup BMA model (+ bias crrction)
    setup_BMA(DREAMPar, ...
Par_info,D,y,VAR);
plugin.BMA = struct('PDF',PDF, ...      % Structure with BMA info AMALGAM
    'VAR',VAR,'D',D_bc,'y',y, ...
    'K',size(D,2));

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_zs';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;                            % # Markov chains
        DREAMPar.T = 5000;                          % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 3;                             % # Markov chains
        DREAMPar.T = 25000;                         % # generations
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.steps = 1000*ones(1,DREAMPar.d);       % # discrete steps
end

options.save = 'yes';       % Save memory during run (for restart trials)

%% Run DREAM-Suite
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar,...
    Par_info,[],options,[],plugin);

%% Postprocess
parset = genparset(chain);
P = parset(end-19999:end,1:DREAMPar.d);
% Weight normalization - this works but please use MODELAVG toolbox
P(:,1:8) = P(:,1:8)./sum(P(:,1:8),2);
