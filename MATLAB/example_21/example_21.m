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
%% Example 21: Soil water flow modeling: Limits of Acceptability                      %%
%%                                                                                    %%
%% Check the following papers                                                         %%
%%   Vrugt, J.A. and K.J. Beven (2018), Embracing equifinality with efficiency:       %%
%%       Limits of Acceptability sampling using the DREAM_{(LOA)} algorithm,  Journal %%
%%       of Hydrology,  559 , pp. 954-971, doi:10.1016/j.jhydrol.2018.02.026          %%
%%   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software %%
%%       package: Theory, concepts, and MATLAB implementation, Environmental Modeling %%
%%       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             %%
%%   Scharnagl, B., J.A. Vrugt, H. Vereecken, and M. Herbst (2011), Bayesian inverse  %%
%%	     modeling of soil water dynamics at the field scale: using prior information  %%
%%	     on soil hydraulic properties, Hydrology and Earth System Sciences, 15,       %%
%%       3043–3059, doi:10.5194/hess-15-3043-2011                                     %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% Problem settings defined by user
DREAMPar.d = 6;                         % Dimension of the problem
DREAMPar.lik = 23;                      % Limits of acceptability likelihood

%% Provide information parameter space and initial sampling
Par_info.initial = 'latin';             % Initial sample from Latin hypercube sampling
Par_info.boundhandling = 'reflect';     % Explicit boundary handling
Par_info.names = {'\theta_{\rm r}','\theta_{\rm s}','\alpha','n',...
    'K_{\rm s}','h_{\rm bot}'};
%				    1	   2	  3	   4	 5		 6
%				[ thetar thetas	alpha  n     Ks    h_bot
Par_info.min =	[ 0.0     0.30	0.02  1.05	0.01   -500 ];   % If "latin"
Par_info.max =	[ 0.1     0.55	0.50  2.50	100/24 -10  ];   % If "latin"

%% Define name of function (.m file) for posterior exploration
Func_name = 'HYDRUS';

%% Now load the data to define the measured data
data = readmatrix('all_observations.xlsx','Sheet','TDR_data', ...
    'Range','A4:AG64'); data = data(:,[4 6:end]);
% Mean soil moisture at time t as observed summary metrics
Meas_info.S = mean(data,'omitnan')';

% Determine limits of acceptability: epsilon for each entry Meas_info.S
D = sort(data);                     % Sort each data column ascending order
[m,n] = size(D);                    % # rows, columns of data matrix
alfa1 = 0.025; alfa2 = 1 - alfa1;   % conf. levels for sig. level alfa = 0.05;
rng = prctile(D,100*[alfa1 alfa2]); % 95% SM range at each measurement time
% rng = D(round([alfa1 alfa2]*m),1:n);
options.epsilon = .5*diff(rng)';    % Define LOAs: nx1 vector

%% Define other data that needs to be ported to HYDRUS
data_hydrus = Load_data;

%% Optional settings
options.IO = 'yes';                 % Input-output writing of model files (only for parallel!)
options.parallel = 'yes';           % Run chains in parallel
options.save = 'yes';               % Save memory of DREAM during trial
%options.restart = 'yes';           % Restart run (no)
options.modout = 'yes';             % Save model output

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_zs';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;                        % # Markov chains
        DREAMPar.T = 2500;                      % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 3;                         % # Markov chains
        DREAMPar.T = 8000;                      % # generations
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.steps = 500*ones(1,DREAMPar.d);    % # discrete steps
end

%% Run DREAM-Suite
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar,...
    Par_info,Meas_info,options,[],data_hydrus);
