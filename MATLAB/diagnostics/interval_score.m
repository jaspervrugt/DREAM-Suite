function [IS,IS_value] = interval_score(LU,obs,p_alfa)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Computes the interval score of the observations given the pdf of the forecast      %%
%% ensemble                                                                           %%
%%                                                                                    %%
%% SYNOPSIS: [IS,IS_value] = interval_score(LU,obs); 	                              %%
%%  where          								                                      %%
%%   LU        [input]  REQUIRED: n x 2 matrix of lower and upper forecasts           %%
%%   obs       [input]  REQUIRED: n x 1 vector of measured data (aka observations)    %%
%%   p_alfa    [input]  REQUIRED: significance level (p_alfa in (0,1))                %%
%%   IS        [output] Mean of interval scores                                       %%
%%   IS_value  [output] n x 1 vector of interval scores                               %%
%%                                                                                    %%
%% Reference: T. Gneiting, and A.E. Raftery (2007), Strictly proper scoring rules,    %%
%%     prediction, and estimation, Journal of the American Statistical Association,   %%
%%     102 (477), 359-378                                                             %%
%%                                                                                    %%
%% Example: fcst = normrnd(0,1,1000,1000);                                            %%
%%  	    obs = rand(1000,1); 						                              %%
%%          IS = interval_score(LU,obs);            		                          %%
%%                                                                                    %%
%% (c) Written by Jasper A. Vrugt, July 2021                                          %%
%% University of California Irvine 						                              %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

if nargin < 3
    p_alfa = 0.05;
end

%% Ensemble size
[n,r] = size(LU);           % Determine the size of the forecast matrix
                            % n measurement times
                            % m ensemble members
if r ~= 2
    error('interval_score:WrongDimension',...
        'The LU vector should have two columns with lower and upper quantiles');
end    
if size(obs,1) ~= n
    error('interval_score:WrongDimension',...
        'The length of the observation vector does not match number of rows of forecast matrix');
end    
if size(obs,2) ~= 1
    error('interval_score:WrongDimension',...
        'The observation vector should have one column only');
end
% Extract the lower and upper quantiles
L = LU(1:n,1); U = LU(1:n,2);
% Compute interval score for all entries at once using inner product
IS_value = (U - L) + 2/p_alfa * (L - obs) .* (obs < L) + ...
    2/p_alfa * (obs - U) .* (obs > U);
% Return mean of all IS values [positive orientation, larger=better]
IS = - mean(IS_value);