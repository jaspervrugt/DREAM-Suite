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
%% Example 15: Approximate Bayesian Computation: Bivariate normal benchmark test      %%
%%                                                                                    %%
%% Check the following papers                                                         %%
%%   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software %%
%%       package: Theory, concepts, and MATLAB implementation, Environmental Modeling %%
%%       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             %%
%%   Sadegh, M., and J.A. Vrugt (2014), Approximate Bayesian computation using Markov %%
%%       chain Monte Carlo simulation: DREAM_(ABC), Water Resources Research,         %%
%%       doi:10.1002/2014WR015386.                                                    %%
%%   Sadegh, M., and J.A. Vrugt (2013), Bridging the gap between GLUE and formal      %%
%%       statistical approaches: approximate Bayesian computation, Hydrology and      %%
%%       Earth System Sciences, 17, 4831–4850                                         %%
%%   Vrugt, J.A., and M. Sadegh (2013), Toward diagnostic model calibration and       %%
%%       evaluation: Approximate Bayesian computation, Water Resources Research, 49,  %%
%%       4335–4345, doi:10.1002/wrcr.20354                                            %%
%%   Turner, B.M., and P.B. Sederberg (2013), Approximate Bayesian computation with   %%
%%       differential evolution, Journal of Mathematical Psychology, In Press.        %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% How many bivariate normal distributions?
Npairs = 10;

%% Problem settings defined by user
DREAMPar.d = 2 * Npairs;    % Dimension of the problem
DREAMPar.lik = 22;          % ABC informal likelihood function

%% Provide information parameter space and initial sampling
Par_info.initial = 'latin';         % Latin hypercube sampling
Par_info.boundhandling = 'fold';    % Explicit boundary handling
Par_info.min = zeros(1,2*Npairs);   % If 'latin', min values
Par_info.max = 10*ones(1,2*Npairs); % If 'latin', max values

%% Define name of function (.m file) for posterior exploration
Func_name = 'ABC_binormal';

%% Lets create the observed summary metrics - the mean (mu) of ten bivariate normals
Meas_info.S = Par_info.min' + rand(DREAMPar.d,1) .* ( Par_info.max - Par_info.min )';

%% Optional settings
%options.rho = inline(' sqrt( 1 / 20 * sum((X - Y).^2)) ');
options.rho = @(X,Y) sqrt( 1 / 20 * sum((X - Y).^2));

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_zs';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;                            % # Markov chains
        DREAMPar.T = 3000;                          % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 3;                             % # Markov chains
        DREAMPar.T = 50000;                         % # generations
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.steps = 100*ones(1,DREAMPar.d);        % # discrete steps
end

%% Run DREAM-Suite
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar,...
    Par_info,Meas_info,options);

%% postprocess results
parset = genparset(chain); parset = parset(end-49999:end,1:2*Npairs);
id_1 = 1:2:2*Npairs; id_2 = id_1 + 1;
figure(100),plot(parset(:,id_1),parset(:,id_2),'r.'); hold on;
plot(Meas_info.S(1:Npairs),Meas_info.S(Npairs+1:2*Npairs),'kx');
axis square

