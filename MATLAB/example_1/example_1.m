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
%% Example 2: w1*N(µ1,Σ1) + w2*N(µ2,Σ2) with µ1 = -5, µ2 = 5, Σ1/Σ2 = identity matrix %%
%%                   and w1 = 1/3 and w2 = 2/3                                        %%
%%                                                                                    %%
%% Check the following papers                                                         %%
%%  Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and J.M.    %%
%%      Hyman (2009), Accelerating Markov chain Monte Carlo simulation by             %%
%%      differential evolution with self-adaptive randomized subspace sampling,       %%
%%      International Journal of Nonlinear Sciences and Numerical Simulation, 10(3),  %%
%%      271-288.                                                                      %%
%%   Vrugt, J.A., H.V. Gupta, W. Bouten and S. Sorooshian (2003), A Shuffled Complex  %%
%%      Evolution Metropolis algorithm for optimization and uncertainty assessment of %%
%%      hydrologic model parameters, Water Resour. Res., 39 (8), 1201,                %%
%%      doi:10.1029/2002WR001642.                                                     %%
%%                                                                                    %%
%% ---------------------------------------------------------------------------------- %%

%% Problem settings defined by user
DREAMPar.d = 2;                        % Dimension of the problem
DREAMPar.lik = 2;                      % Model output is log-likelihood

%% Provide information parameter space and initial sampling
Par_info.initial = 'Normal';           % N2(µ,Σ) initial distribution
Par_info.mu = zeros(1,DREAMPar.d);     % if 'normal', µ-mean of distribution
Par_info.cov = 10 * eye(DREAMPar.d);   % if 'normal', Σ-covariance matrix
% Par_info.norm = 1;                    % sample in normalized [0-1] space

%% Define name of function (.m file) for posterior exploration
Func_name = 'banana_lik';
%% Written as a MATLAB Live Script!

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_zs';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;                        % # Markov chains
        DREAMPar.T = 100000;                    % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 3;                         % # Markov chains
        DREAMPar.T = 20000;                     % # generations
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.min = -100 * ones(1,DREAMPar.d);   % Min value for discrete sampling
    Par_info.max = 100 * ones(1,DREAMPar.d);    % Max value for discrete sampling
    Par_info.steps = 2001*ones(1,DREAMPar.d);   % # discrete steps
end

%% Run DREAM-Suite package
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar,...
    Par_info);
