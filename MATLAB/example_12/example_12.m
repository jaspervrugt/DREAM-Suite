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
%% Example 12: Model-based geostatistics                                              %%
%%                                                                                    %%
%% Check the following paper                                                          %%
%%  Minasny, B., J.A. Vrugt, and A.B. McBratney (2011), Confronting uncertainty in    %%
%%      model-based geostatistics using Markov chain Monte Carlo simulation,          %%
%%      Geoderma, 163, 150-162, doi:10.1016/j.geoderma.2011.03.011                    %%
%%                                                                                    %%
%% ---------------------------------------------------------------------------------- %%

%% Problem settings defined by user
DREAMPar.d = 5;                         % Dimension of the problem
DREAMPar.lik = 2;                       % Model output is log-likelihood

%% Provide information parameter space and initial sampling
Par_info.initial = 'latin';             % Latin hypercube sampling
Par_info.boundhandling = 'reflect';     % Explicit boundary handling
Par_info.min = [0.0 0.0 0.00 0.00 0.00];  % If 'latin', min values
Par_info.max = [100 100 1000 1000 20.0];  % If 'latin', max values

%% Define name of function (.m file) for posterior exploration
Func_name = 'blpmodel';

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_zs';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;                            % # Markov chains
        DREAMPar.T = 2000;                          % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 3;                             % # Markov chains
        DREAMPar.T = 6000;                          % # generations
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.steps = 1000*ones(1,DREAMPar.d);       % # discrete steps
end

%% Run DREAM-Suite
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar, ...
    Par_info);

%% Postprocess the results to generate some fitting results
ParSet = genparset(chain); postproc_variogram
