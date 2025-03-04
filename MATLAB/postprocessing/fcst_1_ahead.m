function varargout = fcst_1_ahead(fx_post,Y_meas,s0,phi1,p_alfa)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This function computes the one-observation ahead forecasts and returns
%% these values for the posterior realizations. Furthermore, the code also
%% returns the LS, CRPS, SS, DDS and IS scoring rules and RLBL, CV, C and W 
%% performance metrics
%%
%% SYNOPSIS: [fx1_mod,S,P] = ...
%%               fcst_1_ahead(fx_post,Y_meas,s0,phi1,p_alfa) 
%%  INPUT: 
%%    fx_post [REQUIRED] nxm matrix posterior simulations (par. unc)
%%    Y_meas  [REQUIRED] nx1 vector of measurements
%%    s0      [REQUIRED] mx1 vector of posterior intercepts
%%    phi1    [REQUIRED] mx1 vector of posterior AR(1)-coefficients
%%    p_alfa  [OPTIONAL] Significance level, 0 < p_alfa < 1
%%  OUTPUT: 
%%    fx1_mod nxm matrix with one-observation ahead forecasts (tot unc)
%%    tot_unc nx2 matrix with lower and upper (1-p_alfa) pred. intervals
%%    S       1x5 vector with LS,CRPS,SS,DDS and IS score rules
%%    P       1x4 vector with RLBL, CV, C and W performance metrics
%% 
%% Written by Jasper A. Vrugt
%% DREAM Package
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    p_alfa = 0.05;
end

calc_method = 'vector';                         % Computation approach    
[n,m] = size(fx_post);                          % Matrix of forecasts of 
                                                % parameter uncertainty
E = Y_meas - fx_post;                           % Calculate residuals
[std_e,~,exitflag] = get_sigma2(s0,E,fx_post);  % Get std_e for all realizations
if exitflag ~= 1, fx1_mod = []; return; end
phi1 = phi1(:)'; std_eps = sqrt(1 - phi1.^2);   % Orient phi and compute std_eps
Eps = normrnd(0,repmat(std_eps,n,1),n,m);       % Draw partial residuals
Err = zeros(n,m);                               % Initialize zeroth matrix
switch calc_method
    case 'loop'
        for t = 2:n
            Err(t,1:m) = std_e(t,1:m) .* ...
                ( phi1 .* E(t-1,1:m)./std_e(t-1,1:m) + Eps(t,1:m) );
        end
    case 'vector'
        Err(2:n,1:m) = std_e(2:n,1:m) .* ( bsxfun(@times,phi1,E(1:n-1,1:m) ./ ...
            std_e(1:n-1,1:m)) + Eps(2:n,1:m) );
end
fx1_mod = fx_post + Err;                        % 1-obs-ahead forecast

%% Now compute prediction limits
p1 = round(p_alfa/2*m); p2 = round((1-p_alfa/2)*m); Ys = sort(fx1_mod,2); 
% Total uncertainty at p_alfa significance level
tot_unc = [Ys(1:n,p1) Ys(1:n,p2)];

%% Now compute score rules
% Logarithmic score
[LS,~,~] = log_score(fx1_mod,Y_meas);
% Continuous Rank Probability Score [Matheson and Winkler, 1976]: minimize
CRPS = crps_jav(fx1_mod,Y_meas);
% Compute the spherical score (zeta = 2)
[SS,~,~] = spherical_rule(fx1_mod,Y_meas);
% Dawid-Sebastiani score: minimize
[DSS,~,~,m_F,s_F] = dawid_sebas(fx1_mod,Y_meas);
% Compute the interval score
[IS,~] = interval_score(tot_unc,Y_meas,p_alfa);
% Return 
S = [LS,CRPS,SS,DSS,IS];

%% Compute summary metrics of forecast density
% Compute p-values (CHECK if < or > than obs)
p_val = p_values(fx1_mod,Y_meas);
% 1. Compute reliability from p-values [Renard et al., 2011]: maximize
[RLBL,~,~] = rlbl(p_val); % --> RLBL: want to maximize
% 2. Compute coefficient of variation [Evin et al., 2013]
CV = 1/n * sum(s_F./m_F);
% 3. Coverage
C = sum(tot_unc(:,1) < Y_meas & tot_unc(:,2) > Y_meas)/n;
% 4. Mean spread/width of 100alfa par/total uncertainty
width_mod = mean(tot_unc(:,2) - tot_unc(:,1));
% Return 
P = [RLBL , CV , C , width_mod ];

% Return arguments
varargout(1:4) = {fx1_mod,tot_unc,S,P};