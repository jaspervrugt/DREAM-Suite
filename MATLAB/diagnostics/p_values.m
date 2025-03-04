function [p_val,num_nan] = p_values(fcst,obs,method)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Computes the p-values of the observations given the members of the ensemble        %%
%%                                                                                    %%
%% SYNOPSIS: [p_val,num_nan] = p_values(fcst,obs) 	                                  %%
%%  where          								                                      %%
%%   fcst      [input]  REQUIRED: n x m matrix of ensemble forecasts 		          %%
%%   obs       [input]  REQUIRED: n x 1 vector of measured data (aka observations)    %%
%%   method    [input]  OPTIONAL: scalar with method (1: loop; 2: direct)             %%
%%   p_val     [output] n x 1 vector with p-values of forecasts                       %%
%%   num_nan   [output] Number of missing values of p-values                          %%
%%                                                                                    %%
%% Note: method 2 is more than 15 times faster than method 1 for larger size matrices %%
%% Example: fcst = normrnd(0,1,1000,1000);                                            %%
%%  	    obs = rand(1000,1); 						                              %%
%%          p_val = p_values(fcst,obs);   					                          %%
%%                                                                                    %%
%% (c) Written by Jasper A. Vrugt, July 2021                                          %%
%% University of California Irvine 						                              %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

if nargin < 3
    method = 2;
end
%% Ensemble size
[n,m] = size(fcst);         % Determine the size of the forecast matrix
                            % n measurement times
                            % m ensemble members
if size(obs,1) ~= n
    error('p_values:WrongDimension',...
        'The length of the observation vector does not match number of rows of forecast matrix');
end
if size(obs,2) ~= 1
    error('p_values:WrongDimension','The observation vector should have one column only');
end

switch method
    case 1 % Loop over time
        p_val = nan(n,1);           
        for t = 1:n
            p_val(t) = sum( fcst(t,1:m) < obs(t) )./m;
        end
    case 2 % Vectorized implementation
        p_val = sum(fcst<obs,2)/m; 
end
num_nan = sum(isnan(p_val));            % Return the number of nan values