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
%% Example 35: Inverse modeling of Parlange's semi-analytic infiltration equation     %%
%%             using the SWIG database of measured infiltration curves                %%
%%                                                                                    %%
%% Check the following papers                                                         %%
%%  Vrugt, J.A., J.W. Hopmans, Y. Gao, M. Rahmati, J. Vanderborght, and               %%
%%      H. Vereecken (2023), The time validity of Philip%s two-term infiltration      %%
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

%% Provide information parameter space and initial sampling
Par_info.initial = 'latin';             % Latin hypercube sampling
Par_info.boundhandling = 'reflect';     % Explicit boundary handling

%% Provide information parameter space and initial sampling
Par_info.initial = 'normal';            % Draw from normal distribution
Par_info.boundhandling = 'fold';        % Explicit boundary handling
Par_info.mu = [10 10 1];                % µ of N(µ,Σ)
Par_info.cov = [2 0 0; 0 2 0; 0 0 0.3]; % Σ of N(µ,Σ)

plugin.model_setup = 1;
switch plugin.model_setup
    case 1
        %          S [cm/h^1/2] Ks [cm/h] beta [-]
        Par_info.min = [ 0   0   0 ];   % Minimum parameter values
        Par_info.max = [ 1e5 1e5 2 ];   % Maximum parameter values
    case 2
        %          S [cm/h^1/2] Ks [cm/h] Ki [cm/h] beta [-]
        Par_info.min = [0   0   0   0]; % Minimum parameter values
        Par_info.max = [1e5 1e5 1e2 2]; % Maximum parameter values
end
DREAMPar.d = numel(Par_info.min);       % Dimension of the problem

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

%% Optional settings
options.modout = 'yes';                 % Return model simulations?
options.parallel = 'no';                % Run each chain on different core
options.save = 'yes';                   % Save workspace DREAM during run
options.print = 'no';                   % No figures printed to screen

%% Load data
load SWIG_1D.mat SWIG_1D                % Load SWIG data set
data_SWIG = SWIG_1D;                    % Rename content of file
n_soil = size(data_SWIG,1);             % Unpack the data

method_to_use = 2;
% [1] synthetic data: Paper 1                           [= example_34]
% [2] SWIG data - time & infiltration form: Paper 2     [= example 35]
% [3] SWIG data - time & infiltration form combined     [= example 35]

%% Which implementation?
switch method_to_use
    case 2
        DREAMPar.lik = 11;                      % Gaussian likelihood
        p95_SWIG = nan(n_soil,5,...             % Init. parameter ranges
            DREAMPar.d,2);
        ML_SWIG = nan(n_soil,DREAMPar.d,2);     % Init. max. likelhd values
        % Add cell array with 95% simulation ranges
        [alg_SWIG,FX95_SWIG,P_SWIG] = deal(cell(n_soil,2));
        [conv_DREAM,MRstat_DREAM] = ...
            deal(nan(n_soil,2));
        cur_dir = pwd;
        for app = 1:2                           % Loop over each approach
            for st = 1:n_soil                   % Loop over each SWIG soil
                dat = data_SWIG{st};            % Unpack data
                plugin.n = size(dat,1);         % Define plugin structure
                switch app
                    case 1 % Minimize cum. inf. residuals
                        plugin.t = dat(:,1); Func_name = 'Haverkamp_I';
                        % Measurement data
                        Meas_info.Y = dat(1:plugin.n,2);    % (in cm)
                    case 2 % Minimize time residuals
                        plugin.I = dat(:,2); Func_name = 'Haverkamp_t';
                        % Measurement data
                        Meas_info.Y = dat(1:plugin.n,1);    % (in hours)
                end
                % Run eDREAM package
                [chain,output,FX,Z,logL] = ...
                    DREAM_Suite(method, ...
                    Func_name,DREAMPar, ...
                    Par_info,Meas_info, ...
                    options,[],plugin);
                alg_SWIG{st,app} = ...          % Store conv. diagnostics
                    output.MR_stat;
                parset = genparset(chain);      % Extract samples chains
                N = size(parset,1);             % # samples
                P = parset(floor(2/3*N)+1:N,... % Burn-in posterior samples
                    1:DREAMPar.d+2);
                M = size(P,1);                  % # posterior samples
                P_SWIG{st,app} = P;             % Store posterior samples
                [~,ii] = max(sum(P(1:M, ...     % Index ML parameters
                    DREAMPar.d+1) + ...
                    P(1:M,DREAMPar.d+2),2));
                ML_SWIG(st,1:3,app) = ...       % Store ML parameters
                    P(ii(1),1:3);
                FX = FX(floor(2/3*N)+1:N, ...   % Get posterior simulations
                    1:plugin.n);
                A = sort(FX); clear fx95;       % Get 95% of simulations
                n1 = floor(0.025*M);            % 2.5% pred. limit
                n2 = floor(0.5*M);              % mean of pred. interval
                n3 = floor(0.975*M);            % 97.5% pred. limit
                fx95(:,1) = A(n1,1:plugin.n)';  % Store 2.5%
                fx95(:,2) = A(n2,1:plugin.n)';  % Store mean
                fx95(:,3) = FX(ii(1), ...       % Store ML simulation
                    1:plugin.n)';
                fx95(:,4) = A(n3,1:plugin.n)';  % Store 97.5% pred. limit
                FX95_SWIG{st,app} = fx95;       % Store for soil
                for j = 1:DREAMPar.d            % Do same for parameters
                    a = sort(P(:,j));
                    p95_SWIG(st,1:5,j,app) = ...
                        [a(n1) mean(a) ...
                        median(a) a(n3) ...
                        std(a) ];
                end
                if output.MR_stat(end,2)<1.2    % Chains converged or not
                    conv = 1;
                else
                    conv = 0;
                end
                conv_DREAM(st,app) = conv;      % Store conv.
                MRstat_DREAM(st,app) = ...      % Multv. R_hat diagnostic
                    output.MR_stat(end,2);
                cd mat_files                    % Change directory
                evalstr = strcat('save', ...
                    {' '},num2str(st),'_', ...
                    num2str(app),['.mat st ' ...
                    'ML_SWIG P fx95 ' ...
                    'conv_DREAM MRstat_DREAM']);
                eval(char(evalstr));
                cd(cur_dir);
            end
        end
        % Save optimal values
        save DREAM_SWIG.mat ML_SWIG p95_SWIG P_SWIG FX95_SWIG ...
            data_SWIG alg_SWIG n_soil MRstat_DREAM conv_DREAM

    case 3
        Func_name = 'Haverkamp_It';             % Specify function name
        DREAMPar.lik = 12;                      % Gaussian likelihood
        sigma_I = 0.1; sigma_t = 0.05;          % Specify measurement error
        p95_IT_SWIG = nan(n_soil,5, ...         % Initlze return matrices
            DREAMPar.d,2);
        ML_IT_SWIG = nan(n_soil,DREAMPar.d,2);
        [alg_IT_SWIG,FX_95_IT_SWIG, ...
            P_IT_SWIG] = deal(cell(n_soil,1));
        [conv_IT_DREAM,MRstat_IT_DREAM] = ...
            deal(nan(n_soil,2));
        cur_dir = pwd;
        for st = 1:n_soil                       % Loop over each SWIG soil
            dat = data_SWIG{st};                % Unpack data
            plugin.n = 2 * size(dat,1);         % Define plugin structure
            n_d = size(dat,1);                  % # data points
            Meas_info.Sigma = ...               % Measurement error std.
                [sigma_t * ones(n_d,1) ; ...
                sigma_I * ones(n_d,1)];
            plugin.t = dat(:,1);                % Time (h)
            plugin.I = dat(:,2);                % Cum. inf. (cm)
            Meas_info.Y = [dat(:,1) ; ...       % Assign measurement data
                dat(:,2)];
            % Run DREAM-Suite
            [chain,output,FX,Z,logL] = ...
                DREAM_Suite(method, ...
                Func_name,DREAMPar, ...
                Par_info,Meas_info, ...
                options,[],plugin);
            alg_IT_SWIG{st} = output.MR_stat;   % Store. conv. diagnostics
            parset = genparset(chain);          % Extract samples chains
            N = size(parset,1);                 % # samples
            P_IT = parset(floor(2/3*N)+1:N, ... % Burn-in posterior samples
                1:DREAMPar.d+2);
            M = size(P_IT,1);                   % # posterior samples
            P_IT_SWIG{st} = P_IT;               % Store posterior samples
            [~,ii] = max(sum(P_IT(1:M, ...      % Index ML parameters
                DREAMPar.d+1) + ...
                P_IT(1:M,DREAMPar.d+2),2));
            ML_IT_SWIG(st,1:3) = ...            % Store ML parameters
                P_IT(ii(1),1:3);
            FX = FX(floor(2/3*N)+1:N, ...       % Get posterior simulations
                1:plugin.n_data);
            A = sort(FX); clear fx95_IT;        % Get 95% of simulations
            n1 = floor(0.025*M);                % 2.5% pred. limit
            n2 = floor(0.5*M);                  % mean of pred. interval
            n3 = floor(0.975*M);                % 97.5% pred. limit
            fx95_IT(:,1) = A(n1,1:plugin.n)';   % Store 2.5%
            fx95_IT(:,2) = A(n2,1:plugin.n)';   % Store mean
            fx95_IT(:,3) = FX(ii(1), ...        % Store ML simulation
                1:plugin.n)';
            fx95_IT(:,4) = A(n3,1:plugin.n)';   % Store 97.5% pred. limit
            FX_95_IT_SWIG{st} = fx95_IT;        % Store for soil
            for j = 1:DREAMPar.d                % Do same for parameters
                a = sort(P_IT(:,j));
                p95_IT_SWIG(st,1:5,j,app) = ...
                    [a(n1) mean(a) ...
                    median(a) a(n3) ...
                    std(a) ];
            end
            if output.MR_stat(end,2)<1.2        % Chains converged or not
                conv = 1;
            else
                conv = 0;
            end
            conv_IT_DREAM(st) = conv;           % Store conv.
            MRstat_IT_DREAM(st) = ...           % Multv. R_hat diagnostic
                output.MR_stat(end,2);
            cd mat_files                        % Change directory
            evalstr = strcat('save', ...
                {' '},num2str(st),'_', ...
                num2str(app), ...
                ['.mat st ML_IT_SWIG P_IT fx95_IT ...'
                'conv_IT_DREAM MRstat_IT_DREAM']);
            eval(char(evalstr));
            cd(cur_dir);
        end
        % Save optimal values
        save DREAM_IT_SWIG.mat ML_IT_SWIG p95_IT_SWIG P_IT_SWIG ...
            FX_95_IT_SWIG data_SWIG alg_IT_SWIG n_soil ...
            MRstat_IT_DREAM conv_IT_DREAM
end
