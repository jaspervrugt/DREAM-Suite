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
%% Example 22: Nash-Cascade model: Limits of Acceptability                            %%
%%                                                                                    %%
%% Check the following papers                                                         %%
%%   Vrugt, J.A. and K.J. Beven (2018), Embracing equifinality with efficiency:       %%
%%       Limits of Acceptability sampling using the DREAM_{(LOA)} algorithm,  Journal %%
%%       of Hydrology,  559 , pp. 954-971, doi:10.1016/j.jhydrol.2018.02.026          %%
%%   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software %%
%%       package: Theory, concepts, and MATLAB implementation, Environmental Modeling %%
%%       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             %%
%%   Nash, J.E. (1960), A unit hydrograph study with particular reference to British  %%
%%       catchments, Proceedings - Institution of Civil Engineers, 17, 249-282        %%
%%   Nash, J.E., J.V. Sutcliffe (1970), River flow forecasting through conceptual     %%
%%       models part I - A discussion of principles, Journal of Hydrology, 10(3),     %%
%%       282-290                                                                      %%
%%   Beven, K.J., and A.M. Binley (1992), The future of distributed models: Model     %%
%%       calibration and uncertainty prediction. Hydrological Processes, 6, 279–298,  %%
%%       doi: 10.1002/hyp.3360060305                                                  %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% Problem settings defined by user
DREAMPar.d = 2;                         % Dimension of the problem
DREAMPar.lik = 23;                      % Limits of acceptability likelihood

%% Provide information parameter space and initial sampling
Par_info.initial = 'latin';             % Latin hypercube sampling
Par_info.boundhandling = 'reflect';     % Explicit boundary handling
Par_info.names = {'k','n'};             % recession constant, # reservoirs
Par_info.min = [ 1   1  ];              % If 'latin', min values
Par_info.max = [ 10  10 ];              % If 'latin', max values

%% Define name of function (.m file) for posterior exploration
Func_name = 'NC_model';
%% Create the synthetic time series
[Y] = NC_model([4 2]); Meas_info.S = normrnd(Y,0.1*Y);

%% Optional settings
options.modout = 'yes';                 % Return model (function) simulations of samples (yes/no)?
options.epsilon = 0.2*Meas_info.S;
options.save = 'yes';

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_zs';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 8;                             % # Markov chains
        DREAMPar.T = 2500;                          % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 3;                             % # Markov chains
        DREAMPar.T = 7000;                          % # generations
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.steps = 90*ones(1,DREAMPar.d);         % # discrete steps
end

%% Run DREAM-Suite
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar,...
    Par_info,Meas_info,options);
