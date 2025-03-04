function [SS,SS_value,num_inf] = spherical_score(fcst,obs,zeta,method,N)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Computes the spherical score of the observations given the pdf of the forecast     %%
%% ensemble                                                                           %%
%%                                                                                    %%
%% SYNOPSIS: [SS,SS_value,num_zero] = spherical_score(fcst,obs); 	                  %%
%%  where          								                                      %%
%%   fcst      [input]  REQUIRED: n x m matrix of ensemble forecasts 		          %%
%%   obs       [input]  REQUIRED: n x 1 vector of measured data (aka observations)    %%
%%   zeta      [input]  OPTIONAL: norm of pseudospherical score (zeta > 1)            %%
%%   method    [input]  OPTIONAL: method of computation [0] st. [1] ref. [2] best     %%
%%   N         [input]  OPTIONAL: number of integration points                        %%
%%   SS        [output] Mean of non-zero spherical scores                             %%
%%   SS_value  [output] n x 1 vector with spherical scores                            %%
%%   num_inf   [output] Number of spherical scores of infinity                        %%
%%                                                                                    %%
%% Reference:                                                                         %%
%%                                                                                    %%
%% Example: fcst = normrnd(0,1,1000,1000);                                            %%
%%  	    obs = rand(1000,1); 						                              %%
%%          DS = spherical_rule(fcst,obs);    				                          %%
%%                                                                                    %%
%% (c) Written by Jasper A. Vrugt, Jan. 2022                                          %%
%% University of California Irvine 						                              %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

if nargin < 5 % Number of integration points of PDF for method 2
    N = 1000;
end
if nargin < 4 % We integrate PDF based on N > 100 points
    method = 2;
end
if nargin < 3
    zeta = 2; % Spherical score
end

%% Ensemble size
[n,m] = size(fcst);         % Determine the size of the forecast matrix
                            % n measurement times
                            % m ensemble members
if size(obs,1) ~= n
    error('spherical_rule:WrongDimension',...
        'The length of the observation vector does not match number of rows of forecast matrix');
end
if size(obs,2) ~= 1
    error('spherical_rule:WrongDimension','The observation vector should have one column only');
end
if zeta < 1
    error('spherical_rule:WrongValue','The pseudospherical score is not defined for zeta < 1');
end

% determine range, minimum and maximum of fcst (or not)
switch method
    case 0
        n2 = 100; Xf = nan(n,n2);
    case {1,2} 
        r = range(fcst,2); dr = r/5; R = r + 2*dr;
        xf_min = min(fcst,[],2) - dr; xf_max = max(fcst,[],2) + dr;
        % get value of n from N
        n2 = 2^ceil(log2(N)); Xf = nan(n,n2);
        for t = 1:n % create xf values
            Xf(t,1:n2) = linspace(xf_min(t),xf_max(t),n2);
        end
end

SS_value = nan(n,1);       % initialize logarithmic scores for each observation
% Loop over each entry of measurement vector
for t = 1:n
    % Get predictive pdf of forecast; evaluate pdf at observed and simulated values
    switch method
        case 0 % Let ksdensity construct the 100 points
            % Evaluate pdf of fcst and return a hundred points
            [f,Xf(t,1:n2)] = ksdensity(fcst(t,1:m)');
            % Evaluate pdf at observation
            f_obs = ksdensity(fcst(t,1:m)',obs(t));
        case 1 % We determine the number and location of the points
            % Evaluate pdf at obs and xf at same time
            f = ksdensity(fcst(t,1:m)',[obs(t) Xf(t,1:n2)]);
            % Now make sure the probabilities sum to one! (= Savage Representation) 
            f = f/sum(f(2:n2+1)); % --> Inconsequential!!
            % Get pdf at observation
            f_obs = f(1); f = f(2:n2+1);
        case 2 % method 1 but using an external kernel density estimator
            % Evaluate pdf at vector of values
            [~,f] = kde(fcst(t,1:m)',n2,Xf(t,1:n2),R(t));
            % Get pdf at observation
            f_obs = interp1(Xf(t,1:n2),f,obs(t));
    end
    % Now integrate the distribution - in theory, between -infty and infty
    I2 = trapz(Xf(t,1:n2),f.^zeta);
    % Compute spherical score: positive orientation [larger = better]
    SS_value(t) = ( f_obs^(zeta-1) ) / ( I2^((zeta-1)/zeta) );
end
ii = SS_value > -inf;                   % Get index of spherical scores > -infinity
num_inf = n - sum(ii);                  % Number of spherical score of -infinity
SS = mean(SS_value(ii));                % Return the mean of the non-missing LS scores