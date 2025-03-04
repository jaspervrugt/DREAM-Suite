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
%% Example 29: w1*N(µ1,Σ1) + w2*N(µ2,Σ2) with µ1 = -5, µ2 = 5, Σ1/Σ2 = identity mtrix %%
%%                   and w1 = 1/3 and w2 = 2/3. Informative prior on parameters       %%
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
DREAMPar.d = 5;                         % Dimension of the problem
DREAMPar.lik = 2;                       % Model output is log-likelihood

%% Provide information parameter space and initial sampling
Par_info.initial = 'prior';             % Use a prior
%% With univariate uniform prior
% % Par_info.prior = { 'unifpdf(-10,10)',...
% %                    'unifpdf(-10,10)',...
% %                    'unifpdf(-10,10)',...
% %                    'unifpdf(-10,10)',...
% %                    'unifpdf(-10,10)' };
%% With univariate informative prior
Par_info.prior = { 'normpdf(2,1)',...
    'gampdf(1,1)',...
    'unifpdf(-10,10)',...
    'unifpdf(-10,10)',...
    'gppdf(1,2,-4)' };
%% With multivariate normal prior
% Par_info.prior = @(x,mu,Sigma) mvnpdf(x,mu,Sigma);
% Par_info.mu = [-2 -2 -2 -2 -2];
% Par_info.Sigma = 5*eye(5);
%% With multivariate uniform prior
% Par_info.prior = @(x,a,b) mvunifpdf(x,a,b);
% Par_info.a = -5*ones(1,5);
% Par_info.b = 5*ones(1,5);

% TRY ALSO WITH
Par_info.boundhandling = 'fold';         % Boundary handling
Par_info.min = -10 * ones(1,DREAMPar.d); % If 'latin', min values
Par_info.max =  10 * ones(1,DREAMPar.d); % If 'latin', max values

%% Define name of function (.m file) for posterior exploration
Func_name = 'mixture_lik';

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;                            % # Markov chains
        DREAMPar.T = 10000;                         % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 3;                             % # Markov chains
        DREAMPar.T = 30000;                         % # generations
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.steps = 500*ones(1,DREAMPar.d);        % # discrete steps
end

%% Run DREAM-Suite
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar,...
    Par_info);
