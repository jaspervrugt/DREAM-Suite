function y = exppdf(x,mu)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Probability density of exponential distribution                       %%
%%                                                                       %%
%% Not sure why this was written - maybe to play with exponential prior  %%
%% © Written by Jasper A. Vrugt, Dec. 2016                               %%
%% University of California Irvine                                       %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if x < 0
    y = 0;
else
    % Overwrite of existing function
    y = mu * exp(-mu * x);
end