function [DSS,DS_value,num_nan,mu_F,sigma_F] = dawidsebas_score(fcst,...
    obs,mu_F,sigma_F)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% This function computes the David-Sebastiani score of the observations given the    %% 
%% forecast pdf ( = prediction of ensemble members) 				                  %%
%%                                                                                    %%
%% SYNOPSIS: [DSS,DS_value,num_nan] = dawidsebas_score(fcst,obs); 	                  %%
%%  where          								                                      %%
%%   fcst      [input]  REQUIRED: n x m matrix of ensemble forecasts 		          %%
%%   obs       [input]  REQUIRED: n x 1 vector of measured data (aka observations)    %%
%%   mu_F      [input]  OPTIONAL: n x 1 vector with mean of CDFs                      %%
%%   sigma_F   [input]  OPTIONAL: n x 1 scalar with standard deviation of CDFs        %%
%%   DSS       [output] Mean of non-missing Dawid-Sebastiani scores                   %%
%%   DS_value  [output] n x 1 vector with Dawid-Sebastiani scores                     %%
%%   num_nan   [output] Number of missing values of CRPS                              %%
%%   mu_F      [input]  n x 1 vector with mean of CDFs                                %%
%%   sigma_F   [input]  n x 1 scalar with standard deviation of CDFs                  %%
%%                                                                                    %%
%% Reference: A.P. Dawid, and P. Sebastiani, Coherent dispersion criteria for optimal %%
%%     experimental design, The Annals of Statistics, 27 (1), p. 65-81, 1999          %%
%%                                                                                    %%
%% Example: fcst = normrnd(0,1,1000,1000);                                            %%
%%  	    obs = rand(1000,1); 						                              %%
%%          DSS = dawid_sebas(fcst,obs);  					                          %%
%%                                                                                    %%
%% (c) Written by Jasper A. Vrugt, July 2021                                          %%
%% University of California Irvine 						                              %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Compute mean and standard deviation of predictive pdf of n forecasts
if nargin < 3 
    mu_F = mean(fcst,2); sigma_F = std(fcst,[],2);
end 

% Ensemble size
[n,m] = size(fcst);         % Determine the size of the forecast matrix
                            % n measurement times
                            % m ensemble members
if size(obs,1) ~= n
    error('dawid_sebas:WrongDimension',...
        'Length of observation vector does not match number of rows forecast matrix');
end    
if size(obs,2) ~= 1
    error('dawid_sebas:WrongDimension','Observation vector should have one column');
end

% Compute Dawid-Sebastiani score
DS_value = - 2*log(sigma_F) - ( (obs - mu_F) ./ sigma_F ).^2;

DSS = mean(DS_value,'omitnan');                 % Compute mean of DS values
num_nan = sum(isnan(DS_value));                 % Return the number of nan values