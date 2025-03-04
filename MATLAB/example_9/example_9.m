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
%% Example 9: Application of spectral likelihood function to conceptual               %%
%%            watershed models using measured discharge data                          %%
%%                                                                                    %%
%% Check the following papers                                                         %%
%%   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of %%
%%       distribution-adaptive likelihood functions: Generalized and universal        %%
%%       likelihood functions, scoring rules and multi-criteria ranking, Journal of   %%
%%       Hydrology, 615, Part B, 2022, doi:10.1016/j.jhydrol.2022.128542.             %%
%%       https://www.sciencedirect.com/science/article/pii/S002216942201112X          %%
%%   Schoups, G., and J.A. Vrugt (2010), A formal likelihood function for parameter   %%
%%       and predictive inference of hydrologic models with correlated,               %%
%%       heteroscedastic and non-Gaussian errors, Water Resources Research, 46,       %%
%%       W10531, doi:10.1029/2009WR008933                                             %%
%%   Schoups, G., J.A. Vrugt, F. Fenicia, and N.C. van de Giesen (2010), Corruption   %%
%%       of accuracy and efficiency of Markov Chain Monte Carlo simulation by         %%
%%       inaccurate numerical implementation of conceptual hydrologic models, Water   %%
%%       Resources Research, 46, W10530, doi:10.1029/2009WR008648                     %%
%%   Montanari, A., and E. Toth (2007), Calibration of hydrological models in the     %%
%%       spectral domain: An opportunity for scarcely gauged basins?, Water Resources %%
%%       Research, 43, W05434, doi:10.1029/2006WR005184                               %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% Problem settings defined by user
DREAMPar.d = 7;                         % Dimension of the problem
DREAMPar.lik = 15;                      % Model output is simulation: Spectral (Whittle) likelihood

%% Parameter names of hmodel
par_names = {'I_{\rm max}','S_{\rm u,max}','Q_{\rm s,max}',...
    '\alpha_{\rm e}','\alpha_{\rm f}','K_{\rm f}','K_{\rm s}'};
%% Provide information parameter space and initial sampling
% parname:      Imax  Smax  Qsmax   alE   alF   Kfast  Kslow
Par_info.min = [ 0.5   10     0    1e-6   -10     0      0    ];    % If 'latin', min values
Par_info.max = [ 10   1000   100   100     10     10    150   ];    % If 'latin', max values
Par_info.initial = 'latin';             % Latin hypercube sampling
Par_info.boundhandling = 'reflect';     % Explicit boundary handling

%% Define name of function (.m file) for posterior exploration
Func_name = 'rainfall_runoff';

%% Define the observed streamflow data
daily_data = load('03451500.dly'); Meas_info.Y = daily_data(731:end,6);

%% Optional settings
options.parallel = 'no';               % Run each chain on different core
options.modout = 'yes';                % Store model simulations
options.save = 'yes';

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_zs';

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
