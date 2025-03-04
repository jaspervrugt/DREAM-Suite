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
%% Example 34: Inverse modeling of Parlange's semi-analytic infiltration equation     %%
%%             using HYDRUS simulated infiltration curves                             %%
%%                                                                                    %%
%% Check the following papers                                                         %%
%%  Vrugt, J.A., J.W. Hopmans, Y. Gao, M. Rahmati, J. Vanderborght, and               %%
%%      H. Vereecken (2023), The time validity of Philip's two-term infiltration      %%
%%      equation: An elusive theoretical quantity? Vadose Zone Journal, e20309,       %%
%%      pp. 1-25, https://doi.org/10.1002/vzj2.20309                                  %%
%%  Vrugt, J.A. and Y. Gao (2022), On the three-parameter infiltration equation of    %%
%%      Parlange et al. (1982): Numerical solution, experimental design, and          %%
%%      parameter estimation, Vadose Zone Journal, 21:e20167, pp. 1-25,               %%
%%      https://doi.org/10.1002/vzj2.20167                                            %%
%%                                                                                    %%
%% ---------------------------------------------------------------------------------- %%

clc; clear; close all hidden;           % clear workspace and figures
plugin.model_setup = 2;                 % which model setup?

%% Problem settings defined by user
DREAMPar.lik = 11;                      % Gaussian likelihood

%% Provide information parameter space and initial sampling
Par_info.initial = 'latin';             % Latin hypercube sampling
Par_info.boundhandling = 'reflect';     % Explicit boundary handling

%% Define name of function (.m file) for posterior exploration
Func_name = 'Haverkamp_I_patch';

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_zs';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;                % # Markov chains
        DREAMPar.T = 10000;             % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 3;                 % # Markov chains
        DREAMPar.T = 20000;             % # generations
end

switch plugin.model_setup
    case 1
        %          S [cm/h^1/2] Ks [cm/h] beta [-]
        Par_info.min = [ 0  0  0 ];     % Minimum parameter values
        Par_info.max = [ 20 50 2 ];     % Maximum parameter values
    case 2
        %          S [cm/h^1/2] Ks [cm/h] Ki [cm/h] beta [-]
        Par_info.min = [0  0  0 0];     % Minimum parameter values
        Par_info.max = [10 50 2 5];     % Maximum parameter values
end
DREAMPar.d = numel(Par_info.min);       % Dimension of the problem

%% Optional settings
options.modout = 'yes';                % Return model simulations (yes/no)?
options.parallel = 'no';               % Run each chain on a different core
options.save = 'yes';                  % Save workspace DREAM during run
options.print = 'no';                  % No figures printed to screen

%% Define the measured infiltration data
load HYDRUS_1D_Data.mat

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                          PREPROCESS DATA                              %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

n_soil = 12;                            % How many soil_types?
plugin.n = 100;                         % # interpolated (t_meas,I_meas)
interpolation_method = 2;               % Interpolation method
% [1] time, [2] cumulative inf.
I_max = 5;                              % Maximum infiltration in cm
% if interpolation_method = 2
interpolation_scheme = 'linear';        % Interpolation scheme?
% 'square_root'/'linear'/'logarithmic'
[data_new,true_pars] = ...              % Pre-processes data
    preprocess_data(data, ...
    Parameters,n_soil,...
    plugin.n,interpolation_method, ...
    interpolation_scheme,I_max);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                         END PREPROCESS DATA                           %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

p95 = nan(n_soil,5,DREAMPar.d);        % Initialize return matrices
ML = nan(n_soil,DREAMPar.d);

for soil_type = 1:12                    % Do each soil & call DREAM Package
    dat = data_new{soil_type};          % Unpack data
    Meas_info.Y = dat(1:plugin.n,2);    % Assign measurement data in cm
    plugin.t = dat(1:plugin.n,1);       % Assign measurement times in h

    % Run DREAM-Suite
    [chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar,...
        Par_info,Meas_info,options,[],plugin);

    parset = genparset(chain);          % Unpack Markov chains
    N = size(chain,1);                  % # samples each chain
    P = parset(floor(1/2*N):N, ...      % Burn-in to get posterior samples
        1:DREAMPar.d+2); M = size(P,1);
    [~,ii] = max(sum(P(1:M, ...         % Index of "best" solution
        DREAMPar.d+1) + ...
        P(1:M,DREAMPar.d+2),2));
    ML(soil_type,1:3) = P(ii(1),1:3);   % Maximum likelihood solution

    for j = 1:DREAMPar.d                % 95% ranges, µ & median solution
        a = sort(P(:,j));
        p95(soil_type,1:5,j) = ...
            [a(round(0.025*M)) ...
            mean(a) median(a) ...
            a(round(0.975*M)) std(a)];
    end

    evalstr = char(strcat(['rename '... % rename file, for figures later
        'DREAM_Suite.mat'],{' '}, ...
        num2str(soil_type),...
        '_',interpolation_scheme,'_', ...
        num2str(plugin.n),'.mat'));
    if ispc                             % statement is platform dependent
        dos(evalstr);
    elseif isunix
        unix(evalstr);
    else
        mac(evalstr);
    end
end

% Save result
save HYDRUS_resuls.mat ML p95 data_new true_pars Parameters data ...
    I_max interpolation_scheme interpolation_method n_soil
