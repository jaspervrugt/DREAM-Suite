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
%% Example 31: AR(2)-model parameter estimation: Test of likelihood functions         %%
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
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% Provide information parameter space and initial sampling
Par_info.initial = 'latin';                 % Latin hypercube sampling
Par_info.boundhandling = 'reflect';         % Explicit boundary handling
[par_names,fpar_mod,parmin_mod,parmax_mod] = deal([ ]);

%% Define name of function (.m file) for posterior exploration
Func_name = 'ar2_model';
%% Determine likelihood function
DREAMPar.lik = 44;

global LV   % Global for likelihood functions 13,14,16,17,44,45
switch DREAMPar.lik
    case 13 %% Normal distribution
        % index:           1   2    3   4
        % parname:        s0  s1  phi1 phi2
        fpar_nuis   =  [  0.1  0    0   0 ];
        parmin_nuis =  [   0   0    0   0 ];
        parmax_nuis =  [   1   1    1   1 ];
        LV.filename = 'Normal'; id_nuis = [3];
    case 14 %% GL (OBSOLETE, DO NOT USE)
        % index:           1    2    3   4   5   6    7    8    9  10  11
        % parname:       std0 std1 beta xi  mu1 phi1 phi2 phi3 phi4 K lambda
        fpar_nuis   =  [  0.1   0    0   1   0   0    0    0    0   0   1  ];
        parmin_nuis =  [   0    0   -1  0.1  0   0    0    0    0   0  0.1 ];
        parmax_nuis =  [   1    1    1  10  100  1    1    1    1   1   1  ];
        LV.filename = 'GL'; id_nuis = [1 2 3 4 6];
    case 16 %% Laplace distribution
        % index:           1   2    3
        % parname:        s0  s1  phi1
        fpar_nuis   =  [  0.1  0    0 ];
        parmin_nuis =  [   0   0    0 ];
        parmax_nuis =  [   1   1    1 ];
        LV.filename = 'Laplace'; id_nuis = [ 3 ];
    case 17 %% SST
        % index:          1   2    3   4    5    6
        % parname:       s0  s1  nu   xi phi1  phi2
        fpar_nuis   =  [ 0.1  0  1e10  1    0    0 ];
        parmin_nuis =  [  0   0    2   0.1  0    0 ];
        parmax_nuis =  [  1   1   100  10   1    1 ];
        LV.filename = 'SL'; id_nuis = [3 4 5];
    case 44 %% GL+
        % index:          1   2    3   4    5    6
        % parname:       s0  s1  beta xi phi1  phi2
        fpar_nuis   =  [  1   0    0   1    0    0 ];
        parmin_nuis =  [  0   0   -1  0.1   0    0 ];
        parmax_nuis =  [  2   1    1   10   1    1 ];
        LV.filename = 'GL_plus'; id_nuis = [3 4 5 6];
    case 45 %% SGT
        % index:          1   2    3   4   5    6    7
        % parname:       s0  s1  labda p   q  phi1 phi2
        fpar_nuis   =  [  1   0    0   2  1e10  0    0 ];
        parmin_nuis =  [  0   0   -1  0.5   2   0    0 ];
        parmax_nuis =  [  2   1    1  100  100  1    1 ];
        LV.filename = 'UL'; id_nuis = [3 4 5 6 7];
end

%% Load data
file_name = 'y.mat';
%% Create autocorrelated data
switch isfile(file_name)
    case 0 % Create measurement data
        M = 5000; phi1 = 0.7; phi2 = 0.2; y(1:2,1:2) = zeros(2,2); %randn(2,2);
        for i = 1:2
            switch i
                case 1 % innovations = 'SEP';
                    beta = -0.5; xi = 5;
                    rnd = SEPrnd(beta,xi,[M 1]);
                case 2 % innovations = 'SGT'
                    lambda = -0.5; p = 1.2; q = 5.7;
                    rnd = SGTrnd(0,1,lambda,p,q,[M 1]);
            end
            for t = 3:M
                y(t,i) = phi1 * y(t-1,i) + phi2 * y(t-2,i) + rnd(t,1);
            end
        end
        save y.mat y M
    case 1
        load(file_name)
end
plugin.N = M;
switch DREAMPar.lik
    case 44
        Meas_info.Y = y(:,1);
    case 45
        Meas_info.Y = y(:,2);
end
Meas_info.sigma2 = 'constant';   % Constant measurement error variance

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
nuis_names = nuis_var_names(DREAMPar);              % Extract names nuis. variables
names = [par_names,nuis_names];                     % Combine with parameter names
n_modpar = numel(fpar_mod);                         % # model parameters
LV.id_vpar = [1:n_modpar , id_nuis + n_modpar ];    % Model parameter + nuisance variable selection
LV.fpar = [fpar_mod fpar_nuis];                     % Merge default values parameters and nuisance variables
parmin = [parmin_mod parmin_nuis];                  % Merge the parameter names and nuisance variables
parmax = [parmax_mod parmax_nuis];
Par_info.min = parmin(LV.id_vpar);                  % Min/max values of parameter selection
Par_info.max = parmax(LV.id_vpar);
Par_info.names = names(LV.id_vpar);                 % Names of parameters and nuisance variables
DREAMPar.d = size(Par_info.min,2);                  % # parameters

%% Optional DREAM settings
options.parallel = 'no';            % This example runs in parallel
options.save = 'yes';               % Save DREAM workspace during run
options.modout = 'no';              % Store model simulations

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_zs';
switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;            % # Markov chains
        DREAMPar.T = 3000;          % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 3;             % # Markov chains
        DREAMPar.T = 10000;         % # generations
end

%% Run DREAM-Suite
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar,...
    Par_info,Meas_info,options,[],plugin);
