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
%% Example 2: N(µ,Σ) with µ = 0 and Σ = dxd matrix with 1:d on main diagonal and      %%
%%                   0.5 correlation between dimensions                               %%
%% Check: func_normal for details                                                     %%
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
DREAMPar.thinning = 2;                  % Only store each 10th sample
DREAMPar.lik = 2;                       % Model output is log-likelihood

%% Provide information parameter space and initial sampling
Par_info.initial = 'latin';             % Latin hypercube sampling
Par_info.min = -15 *ones(1,DREAMPar.d); % If 'latin', min values
Par_info.max = 15 *ones(1,DREAMPar.d);  % If 'latin', max values
Par_info.boundhandling = 'none';

%% Define name of function (.m file) for posterior exploration
Func_name = 'normal_lik';

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
%method = 'mtdream_zs';
%method = 'dream_d'
method = 'dream_zs'
switch method
    case {'dream','dream_d'}
        DREAMPar.N = DREAMPar.d;                    % # Markov chains
        DREAMPar.T = 10000;                         % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 10;                            % # Markov chains
        DREAMPar.T = 25000;                         % # generations
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.min = -50*ones(1,DREAMPar.d);          % Min value for discrete sampling
    Par_info.max =  50*ones(1,DREAMPar.d);          % Max value for discrete sampling
    Par_info.steps = 1000*ones(1,DREAMPar.d);       % # discrete steps
end

options.save = 'yes';       % save file during trial    
options.print = 'no';       % no visual output to screen

%% Run DREAM-Suite
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar, ...
    Par_info,[],options);

%% Check mean and variance of MCMC samples
P = genparset(chain); P = P(end-24999:end,1:DREAMPar.d+2);

figure(100);
plot(mean(P),'b+','markersize',8,'linewidth',2);
line([0 DREAMPar.d],[0 0],'color','r','linewidth',2);
xlabel('Dimension','fontsize',16);
ylabel('Mean, $\mu_{i}$','fontsize',16);
set(gca,'fontsize',16);
[~,hleg] = legend('\color{blue} estimated','\color{red} true', ...
    'box','off','fontsize',16);
set(hleg(4),'markersize',15,'linewidth',3);
set(hleg(5),'linewidth',3);
axis square
title('Target mean against estimates from MCMC','fontsize',16);

figure(101);
plot(diag(cov(P)),'b+','markersize',8,'linewidth',2);
line([0 DREAMPar.d],[0 DREAMPar.d],'color','r','linewidth',2);
xlabel('Dimension','fontsize',16);
ylabel('Variance, $\sigma^{2}_{i}$','fontsize',16);
set(gca,'fontsize',16);
[~,hleg] = legend('\color{blue} estimated','\color{red} true', ...
    'box','off','fontsize',16);
set(hleg(4),'markersize',15,'linewidth',3);
set(hleg(5),'linewidth',3);
axis square
title('Target variance against estimates from MCMC','fontsize',16);
