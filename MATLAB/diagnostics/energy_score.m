function [mean_ES,ES_value,num_nan] = energy_score(fcst,obs,beta)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% The ENERGY SCORE is a generalization of the continuous rank probability score      %%
%%                                                                                    %%
%% SYNOPSIS: [mean_ES,ES_value,num_nan] = energy_score(fcst,obs);                     %%
%%           [mean_ES,ES_value,num_nan] = energy_score(fcst,obs,beta);                %%
%%  where                                                                             %%
%%   fcst      [input]  REQUIRED: n x m matrix of ensemble forecasts                  %%
%%   obs       [input]  REQUIRED: n x 1 vector of measured data (aka observations)    %%
%%   beta      [input]  REQUIRED: scalar with value of beta                           %%
%%   mean_ES   [output] Mean of non-missing ES values                                 %%
%%   ES_value  [output] n x 1 vector with ES values                                   %%
%%   num_nan   [output] Number of missing values of ES                                %%
%%                                                                                    %%
%% Reference:                                                                         %%
%%   T. Gneiting, and A.E. Raftery (2007), Strictly proper scoring rules, prediction, %%
%%      and estimation, Journal of the American Statistical Association,              %%
%%      102 (477), 359-378                                                            %%
%%                                                                                    %%
%% Example: fcst = normrnd(0,1,1000,1000);                                            %%
%%  	    obs = rand(1000,1); 						                              %%
%%          [mean_ES,ES_values,num_nan] = energy_score(fcst,obs);        			  %%
%%                                                                                    %%
%% (c) Written by Jasper A. Vrugt, Aug. 2022                                          %%
%% University of California Irvine 						                              %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% Check input arguments
if nargin < 3
    beta = 1;
end

%% Ensemble size
[n,m] = size(fcst);         % Determine the size of the forecast matrix
% n measurement times
% m ensemble members
if size(obs,1) ~= n
    error('ENERGY_SCORE:WrongDimension',...
        'The length of the observation vector does not match number of rows of forecast matrix');
end
if size(obs,2) ~= 1
    error('ENERGY_SCORE:WrongDimension','The observation vector should have one column only');
end

ES_value = nan(n,1);        % initialize ES values
ys = sort(fcst,2);          % Sort entries in all rows of fcst in increasing order

%% Approach B1: For empirical CDF!
for t = 1:n
    first_term = sum ( (abs(ys(t,1:m) - obs(t) ) ).^beta );
    % second_term = 0;
    % for i = 1:m
    %    for j = 1:m
    %        second_term = second_term + (abs(ys(t,i) - ys(t,j) ) ).^beta;
    %    end
    % end
    second_term = 0;
    for i = 1:m
        second_term = second_term + sum ( (abs(ys(t,1:m) - ys(t,i) ) ).^beta );
    end
    ES_value(t) = 1/2 * 1/m^2 * second_term - 1/m * first_term;
end

mean_ES = mean(ES_value,'omitnan');                   % Compute mean of energy score
num_nan = sum(isnan(ES_value));                       % Return the number of nan values
% Positive orientation [larger = better]