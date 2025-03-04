function [LS,LS_value,num_zero] = log_score(fcst,obs,mu_F,sigma_F)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Computes the logarithmic score of the observations given the pdf of the forecast   %%
%% ensemble                                                                           %%
%%                                                                                    %%
%% SYNOPSIS: [LS,LS_value,num_zero] = log_score(fcst,obs,mu_F,sigma_F); 	          %%
%%  where          								                                      %%
%%   fcst      [input]  REQUIRED: n x m matrix of ensemble forecasts 		          %%
%%   obs       [input]  REQUIRED: n x 1 vector of measured data (aka observations)    %%
%%   mu_F      [input]  OPTIONAL: n x 1 vector with mean of CDFs                      %%
%%   sigma_F   [input]  OPTIONAL: n x 1 scalar with standard deviation of CDFs        %%
%%   LS        [output] Mean of non-zero logarithmic scores                           %%
%%   LS_value  [output] n x 1 vector with logarithmic scores                          %%
%%   num_zero  [output] Number of logarithmic scores of zero                          %%
%%                                                                                    %%
%% Reference: I.J. Good (1952), Rational decisions, Journal of the Royal Statistical  %%
%%     Society, Series B (Statistical Methodology), 14 (1), 107–114                   %%
%%                                                                                    %%
%% Example: fcst = normrnd(0,1,1000,1000);                                            %%
%%  	    obs = rand(1000,1); 						                              %%
%%          LS = log_score(fcst,obs);    					                          %%
%%                                                                                    %%
%% (c) Written by Jasper A. Vrugt, July 2021                                          %%
%% University of California Irvine 						                              %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

switch nargin < 3
    case 1, G = 0; 
    case 0, G = 1;
end
%% Ensemble size
[n,m] = size(fcst);         % Determine the size of the forecast matrix
                            % n measurement times
                            % m ensemble members
if size(obs,1) ~= n
    error('log_score:WrongDimension',...
        'The length of the observation vector does not match number of rows of forecast matrix');
end    
if size(obs,2) ~= 1
    error('log_score:WrongDimension','The observation vector should have one column only');
end

switch G
    case 0
        LS_value = nan(n,1);       % initialize logarithmic scores for each observation
        num_zero = 0;              % number of zero-densities
        % Loop over each entry of measurement vector
        for t = 1:n 
            % get predictive pdf of forecast; evaluate pdf at observed value = py
            py = ksdensity(fcst(t,1:m)',obs(t)); %,'Support','positive');
            if py < realmin
                py = realmin; num_zero = num_zero + 1;
            end
            % compute the logarithmic score
            LS_value(t) = log(py);
        end
    case 1
        py = normpdf(obs,mu_F,sigma_F);
        ii = py < realmin; py(ii) = realmin;
        num_zero = sum(ii);
        % Compute the logarithmic score
        LS_value = log(py);
end
% Positive orientation: larger is better
LS = mean(LS_value(LS_value > -inf),'omitnan');        % Return the mean of the non-missing LS scores