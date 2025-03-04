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
%% Example 10: w1*N(µ1,Σ1) + w2*N(µ2,Σ2) with µ1 = -5, µ2 = 5, Σ1/Σ2 =identity matrix %%
%%                   and w1 = 1/3 and w2 = 2/3. Informative prior                     %%
%%                                                                                    %%
%% Check the following papers                                                         %%
%%  Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and J.M.    %%
%%      Hyman (2009), Accelerating Markov chain Monte Carlo simulation by             %%
%%      differential evolution with self-adaptive randomized subspace sampling,       %%
%%      International Journal of Nonlinear Sciences and Numerical Simulation, 10(3),  %%
%%      271-288.                                                                      %%
%%  Ter Braak, C.J.F., and J.A. Vrugt (2008), Differential Evolution Markov Chain     %%
%%      with snooker updater and fewer chains, Statistics and Computing,              %%
%%      10.1007/s11222-008-9104-9.                                                    %%
%%                                                                                    %%
%% ---------------------------------------------------------------------------------- %%

%% Problem settings defined by user
DREAMPar.d = 2;                         % Dimension of the problem
DREAMPar.lik = 1;                       % Model output is likelihood

%% Provide information parameter space and initial sampling
Par_info.initial = 'normal';            % Multinormal initial sampling distribution
Par_info.mu = zeros(1,DREAMPar.d);      % If 'normal', define mean of distribution
Par_info.cov = eye(DREAMPar.d);         % If 'normal', define covariance matrix
Par_info.prior = @(x,a,b) mvnpdf(x,a,b);% Specify a bi-variate prior distribution
Par_info.a = [-2 -2];                   % Mean of prior distribution
Par_info.b = eye(2);                    % Covariance of prior distribution

%% Define name of function (.m file) for posterior exploration
Func_name = 'mixture_lik';

% So the mixture models has two modes at -5 and 5; with the specified prior
% distribution the mode around 5 should disappear. You can compare the
% theoretical distribution with that derived from DREAM

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_zs';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;                            % # Markov chains
        DREAMPar.T = 2000;                          % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 3;                             % # Markov chains
        DREAMPar.T = 7000;                          % # generations
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.min = [ -20 -20 ];                     % Min value for discrete sampling
    Par_info.max = [  20  20 ];                     % Max value for discrete sampling
    Par_info.steps = 1000*ones(1,DREAMPar.d);       % # discrete steps
end

%% Run DREAM-Suite
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar, ...
    Par_info);
