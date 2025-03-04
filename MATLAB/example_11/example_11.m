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
%% Example 11: Target is a d-variate t-distribution with df degrees of freedom        %%
%%             and correlation matrix R                                               %%
%%                                                                                    %%
%% Check the following papers                                                         %%
%%  Ter Braak, C.J.F., and J.A. Vrugt (2008), Differential Evolution Markov Chain     %%
%%      with snooker updater and fewer chains, Statistics and Computing,              %%
%%      10.1007/s11222-008-9104-9.                                                    %%
%%                                                                                    %%
%% ---------------------------------------------------------------------------------- %%

%% Problem specific parameter settings
DREAMPar.d = 25;                        % Dimension of the problem
DREAMPar.lik = 2;                       % Model output is log-likelihood

%% Provide information parameter space and initial sampling
Par_info.initial = 'latin';                 % Latin hypercube sampling
Par_info.min = -5 * ones(1,DREAMPar.d);     % Lower bound parameter values
Par_info.max = 15 * ones(1,DREAMPar.d);     % Upper bound parameter values

%% Define name of function (.m file) for posterior exploration
Func_name = 'mvt_lik';

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_zs';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 20;                            % # Markov chains
        DREAMPar.T = 5000;                          % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 3;                             % # Markov chains
        DREAMPar.T = 30000;                         % # generations
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.steps = 1000*ones(1,DREAMPar.d);       % # discrete steps
end

%% Run DREAM-Suite
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar, ...
    Par_info);
