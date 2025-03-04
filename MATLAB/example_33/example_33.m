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
%% Example 33: Two-dimensional target distribution: from Haario et al. 1999           %%
%%                                                                                    %%
%% Check the following papers                                                         %%
%%   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software %%
%%       package: Theory, concepts, and MATLAB implementation, Environmental Modeling %%
%%       and Software, 75, pp. 273-316, https://doi.org/10.1016/j.envsoft.2015.08.013 %%
%%  Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and J.M.    %%
%%      Hyman (2009), Accelerating Markov chain Monte Carlo simulation by             %%
%%      differential evolution with self-adaptive randomized subspace sampling,       %%
%%      International Journal of Nonlinear Sciences and Numerical Simulation, 10(3),  %%
%%      271-288                                                                       %%
%%   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and B.A. Robinson (2008), %%
%%       Treatment of input uncertainty in hydrologic modeling: Doing hydrology       %%
%%       backward with Markov chain Monte Carlo simulation, Water Resources Research, %%
%%       44, W00B09, https://doi.org/10.1029/2007WR006720                             %%
%%   Vrugt, J.A., H.V. Gupta, W. Bouten and S. Sorooshian (2003), A Shuffled Complex  %%
%%      Evolution Metropolis algorithm for optimization and uncertainty assessment of %%
%%      hydrologic model parameters, Water Resour. Res., 39 (8), 1201,                %%
%%      https://doi.org/10.1029/2002WR001642                                          %%
%%   Haario, H., E. Saksman, and J. Tamminen (1999), Adaptive proposal distribution   %%
%%      for random walk Metropolis algorithm. Computational Statistics, 14, 375–395,  %%
%%      https://doi.org/10.1007/s001800050022                                         %%
%%                                                                                    %%
%% ---------------------------------------------------------------------------------- %%

%% Problem settings defined by user
DREAMPar.d = 2;                         % Dimension of the problem
DREAMPar.lik = 1;                       % Model output is likelihood

%% Provide information parameter space and initial sampling
Par_info.initial = 'latin';             % latin hypercube initial chain states
Par_info.min = [-18 -3];                % lower bound
Par_info.max = [ 18  3];                % upper bound
Par_info.boundhandling = 'fold';        % enforce bounds using folding
Par_info.norm = 1;                      % sample in normalized space (= test)

%% Define name of function (.m file) for posterior exploration
Func_name = 'rectangle_lik';
%% Written as a MATLAB Live Script!

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_zs';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;            % # Markov chains
        DREAMPar.T = 10000;         % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 3;             % # Markov chains
        DREAMPar.T = 20000;         % # generations
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.steps = [7200 600];   % # discrete steps
end

%% Run DREAM-Suite
[chain,output,FX] = DREAM_Suite(method,Func_name,DREAMPar,...
    Par_info);

%% Now do some postprocessing, specifically for this example
P = genparset(chain); Post = P(end-29999:end,:);

% If sampled in normalized space then back transform to original space
if ~isfield(Par_info,'norm')
    Par_info.norm = 0; % sampled in regular space
end

switch Par_info.norm
    case 1 % must back transform
        Par_info.minun = Par_info.min; Par_info.maxun = Par_info.max;
        Par_info.min = [0 0]; Par_info.max = [1 1];
        Post(:,1:2) = X_unnormalize(Post(:,1:2),Par_info);
    otherwise % do nothing
end

% Look only at first parameter as 2nd is equal to prior (inconsequential)
edges = [-18:.25:-1/2 -0.3:0.2:0.3 0.5:.25:18];
% Make a histogram of the posterior samples
[N,edges] = histcounts(Post(:,1),edges,'normalization','pdf');
% Midpoints of bins
bin = 1/2*(edges(1:end-1)+edges(2:end));
% Does the empirical density integrate to one?
trapz(bin,N)
% Sufficiently close!

% Now plot a histogram of marginal distribution according to MCMC samples
figure(100),bar(bin,N,'b'); hold on
% Now plot true target density
x = linspace(-18,18,3601);
f = @(x) 35*(x>=-0.5 & x<=0.5) + 1;
Z = trapz(x,f(x));
% 1d scaling factor = int_x1 -18:18 but without [-1/2:1/2] = 35 * 1 = 35
%                             -.5:.5 is 36 * 1 = 36. Thus, integral is 71
plot(x,f(x)/Z,'r','linewidth',2);
legend('$\;$ MCMC-estimate',...
    '$\;$ True target','interpreter','latex','fontsize',16);
% Why the small mismatch in the rectangular area of the peak? (DREAM_ZS)
% A result in part of bin choice. but also the probability mass by the MCMC
% method may not be exactly similarly distributed as the target.

% Simple check is that about 50% of the probability mass of the distribution
% is in the tails: 18:-1/2 and 1/2:18 has integral of 35; -.5 to .5 of 36.
% Now check what percentage of posterior samples is outside S: -.5 to .5
M = size(P,1); out_MCMC = 100*(sum(P(:,1)<-0.5 | P(:,1)>.5))/M;
out_true = 100 * 35/(35+36);
disp([out_MCMC out_true])
% Not a reason for concern: can play with # samples, # chains,
% # post. samples, etc. Even try using another DREAM sampler, etc.
