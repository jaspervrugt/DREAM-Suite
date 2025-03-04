function [RLBL,eCDF,uCDF] = rlbl(p_val)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Computes the time-averaged reliability of the forecast distribution of ensemble    %%
%%                                                                                    %%
%% SYNOPSIS: [rlbl,eCDF,uCDF] = rlbl(p_val)      	                                  %%
%%  where          								                                      %%
%%   p_val     [input]  n x 1 vector with p-values of forecasts                       %%
%%   RLBL      [output] scalar with time-averaged reliability of forecasts            %%
%%   eCDF      [output] n x 1 vector with CDF (= sorted p_values) of forecasts        %%
%%   uCDF      [output] n x 1 vector with cdf of uniform distribution                 %%
%%                                                                                    %%
%% Example: fcst = normrnd(0,1,1000,1000);                                            %%
%%  	    obs = rand(1000,1); 						                              %%
%%          p_val = p_values(fcst,obs);   					                          %%
%%          [RLBL,eCDF,uCDF] = rlbl(p_val);                                           %% 
%%                                                                                    %%
%% (c) Written by Jasper A. Vrugt, July 2021                                          %%
%% University of California Irvine 						                              %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Computes reliability of ensemble using p-values

n = numel(p_val);                       % number of p-values ( = number of observations )
eCDF = sort(p_val);                     % empirical CDF (sort p-values in ascending order)
uCDF = (1:n)'./n;                       % cdf of uniform distribution
RLBL = 1 - 2/n * sum(abs(uCDF - eCDF)); % compute reliability