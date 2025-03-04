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
%% Example 13: Nash-Cascade series of linear reservoirs: synthetic study              %%
%%                                                                                    %%
%% Check the following papers                                                         %%
%%   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of %%
%%       distribution-adaptive likelihood functions: Generalized and universal        %%
%%       likelihood functions, scoring rules and multi-criteria ranking, Journal of   %%
%%       Hydrology, 615, Part B, 2022, doi:10.1016/j.jhydrol.2022.128542.             %%
%%       https://www.sciencedirect.com/science/article/pii/S002216942201112X          %%
%%   Nash, J.E. (1960), A unit hydrograph study with particular reference to British  %%
%%       catchments, Proceedings - Institution of Civil Engineers, 17, 249-282        %%
%%   Nash, J.E., J.V. Sutcliffe (1970), River flow forecasting through conceptual     %%
%%       models part I - A discussion of principles, Journal of Hydrology, 10(3),     %%
%%       282-290                                                                      %%
%%   Beven, K.J., and A.M. Binley (1992), The future of distributed models: Model     %%
%%       calibration and uncertainty prediction. Hydrological Processes, 6, 279–298,  %%
%%       doi: 10.1002/hyp.3360060305                                                  %%
%%   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and B.A. Robinson (2008), %%
%%       Treatment of input uncertainty in hydrologic modeling: Doing hydrology       %%
%%       backward with Markov chain Monte Carlo simulation, Water Resources Research, %%
%%       44, W00B09, doi:10.1029/2007WR006720                                         %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% Problem settings defined by user
DREAMPar.d = 1;                         % Dimension of the problem
DREAMPar.lik = 12;                      % Gaussian likelihood with measurement error

%% Provide information parameter space and initial sampling
Par_info.initial = 'latin';             % Latin hypercube sampling
Par_info.boundhandling = 'reflect';     % Explicit boundary handling
Par_info.min =  1;                      % If 'latin', min values
Par_info.max = 100;                     % If 'latin', max values

%% Define name of function (.m file) for posterior exploration
Func_name = 'Nash_Cascade';

%% Create the data: synthetic time series
y = Nash_Cascade(2);
%% Now define heteroscedastic measurement error
Meas_info.Sigma = max( 1/5 * y , 1e-2);
%% Observed data: Model simulation with heteroscedastic measurement error
Meas_info.Y = normrnd(y,Meas_info.Sigma);

%% Optional settings
options.modout = 'yes';                % Return model (function) simulations of samples (yes/no)?

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_dzs';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;                            % # Markov chains
        DREAMPar.T = 2000;                          % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 3;                             % # Markov chains
        DREAMPar.T = 7000;                          % # generations
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.steps = 1000*ones(1,DREAMPar.d);       % # discrete steps
end

%% Run DREAM-Suite
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar,...
    Par_info,Meas_info,options);
