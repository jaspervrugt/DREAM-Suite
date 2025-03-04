function nuis_names = nuis_var_names(DREAMPar)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Returns the names of nuisance variables selected likelihood function  %%
%%                                                                       %%
%%  SYNOPSIS: nuis_names = nuis_var_names(DREAMPar)                      %%
%% where                                                                 %%
%%  DREAMPar    [input] Structure with algorithmic variables             %%
%%  nuis_names  [outpt] Cell array with names of nuisance variables      %%
%%                                                                       %%
%% © Written by Jasper A. Vrugt, Dec. 2016                               %%
%% DREAM Package                                                         %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

nuis_names = [];        % Return argument, initialize empty

% Switch command to distinguish between likelihood functions
switch DREAMPar.lik
    case 13 % Normal AR(2)-likelihood with (non)constant residual variance
            % Vrugt et al., JOH, 2022 
            % Common in statistical literature
        nuis_names = {'s_{0}','s_{1}','\phi_{1}','\phi_{2}'};
    case 14 % Generalized Likelihood function [= OBSOLETE]
            % Schoups and Vrugt, WRR, 2010
        nuis_names = {'s_{0}','s_{1}','\beta','\xi','\mu_{1}',...
            '\phi_{1}','\phi_{2}','\phi_{3}','\phi_{4}','K','\lambda'};
    case 16 % Laplace AR(1)-likelihood with (non)constant residual variance
            % Vrugt et al., JOH, 2022
            % Relatively common in statistical literature            
        nuis_names = {'s_{0}','s_{1}','\phi_{1}'};
    case 17 % Student AR(2)-likelihood with (non)constant residual variance
            % Scharnagl et al., HESS-D, 2015
            % Vrugt et al., JOH, 2022
        nuis_names = {'s_{0}','s_{1}','\nu','\xi','\phi_{1}','\phi_{2}'};
    case 44 % Generalized likelihood ++, AR(2) & (non)constant residual var.
            % Vrugt et al., JOH, 2022
        nuis_names = {'s_{0}','s_{1}','\beta','\xi','\phi_{1}','\phi_{2}'};
    case 45 % Universal likelihood, AR(2) & (non)constant residual variance
            % Vrugt et al., JOH, 2022
        nuis_names = {'s_{0}','s_{1}','\lambda','p','q','\phi_{1}',...
            '\phi_{2}'};
end

end