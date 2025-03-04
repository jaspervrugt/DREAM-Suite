function [std_e,s1,exflag] = s1_phantom(s0,e,y)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% This function calculates candidate points using discrete proposal distribution     %%
%%                                                                                    %%
%% SYNOPSIS: [std_e,s1] = s1_phantom(s0,e,y)                                          %%
%%  where                                                                             %%
%%   s0         [input] intercept of the measurement error function                   %%
%%   e          [input] nxm matrix of residuals of N parameter vectors                %%
%%   y          [input] measured/simulated y values                                   %%
%%   std_e      [outpt] nxN matrix of measurement error standard deviations           %%
%%   s1         [outpt] 1xN vector of estimated slopes of measurement error function  %%
%%   exflag     [outpt] Exitflag: 1 if converged, otherwise not converged             %%
%%                                                                                    %%
%% More information can be found in the following references                          %%
%%   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of %%
%%       distribution-adaptive likelihood functions: Generalized and universal        %%
%%       likelihood functions, scoring rules and multi-criteria ranking, Journal of   %%
%%       Hydrology, 615, Part B, 2022, doi:10.1016/j.jhydrol.2022.128542.             %%
%%       https://www.sciencedirect.com/science/article/pii/S002216942201112X          %%
%%   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software %%
%%       package: Theory, concepts, and MATLAB implementation, Environmental Modeling %%
%%       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             %%
%%   Schoups, G., and J.A. Vrugt (2010), A formal likelihood function for parameter   %%
%%       and predictive inference of hydrologic models with correlated,               %%
%%       heteroscedastic and non-Gaussian errors, Water Resources Research, 46,       %%
%%       W10531, doi:10.1029/2009WR008933                                             %%
%%                                                                                    %%
%% © Written by Jasper A. Vrugt, April 2017                                           %%
%% University of California Irvine                                                    %%
%%                                                                                    %%
%% Inspired in part by discussion with Mario Hernandez and the following manuscript   %%
%%   Hernández-López, M.R., and F. Francés (2017), Bayesian joint inference of        %%
%%       hydrological and generalized error models with the enforcement of total laws %%
%%       Hydrology and Earth System Sciences Discussions, pp. 1-40,                   %%
%%       doi:10.5194/hess-2017-9                                                      %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

[n,N] = size(e);                        % # elements, # residual vectors
s1 = nan(1,N); std_e = nan(n,N);        % Initialize slope and std_e
if size(y,2) == 1
    if N > 1, y = repmat(y,1,N); end        % Make N copies of y
end
for ii = 1:N                            % Do each residual vector
    [s1(ii),~,exflag] = fzero(@(x) ...      % Root finding: Newton method
        eval_func(x,s0(ii),e(:,ii),...      % so that studentized residuals
        y(:,ii)),0.2,...                    % have a unit variance
        optimset('Display','off'));     
    if (s1(ii) < 0 || exflag < 0)           % If slope s1 is improper then
        [s1(ii),~,exflag] = ...                 % minimize err: Nelder Mead
            fminsearch(@(x) get_opt(x,...       % algorithm using quadratic   
            s0(ii),e(:,ii),y(:,ii)),0.2,...     % objective function with a  
            optimset('Display','off'));         % penalty term for s1 < 0
    end
    if s1(ii) < 0, warning('s1 < 0'); end   % Warning s1 < 0 (= unlikely)
    std_e(:,ii) = s1(ii)*y(:,ii) + s0(ii);  % nxN matrix of meas. err. stds.
end                                     % End of minimization loop

end
% <><><><><><><><><><><><> End of primary function <><><><><><><><><><><><>

% <><><><><><><><><><><><><> Secondary functions <><><><><><><><><><><><><>
% SUBROUTINE 1
function err = eval_func(s1,s0,e,y)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Compute error deviation as 1 minus variance of studentized residuals  %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

std_e = s1*y + s0;              % measurement error function (proposed)
e_n = e./std_e;                 % normalized residuals (proposed)
err = 1 - var(e_n,'omitnan');   % error deviation from unit variance

end

% <><><><><><><><><><><><><> Secondary functions <><><><><><><><><><><><><>
% SUBROUTINE 2
function err = get_opt(s1,s0,e,y)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% More sophisticated variant of eval_func with penalty for magnitude s1 %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

w = 100;                        % set weight for regularization (= penalty)
std_e = s1*y + s0;              % measurement error function (proposed)
e_n = e./std_e;                 % normalized residuals (proposed)
err = 1 - var(e_n,'omitnan');   % error (distance) from unit variance
pen = w*(s1 < 0)*abs(s1);       % compute penalty
err = err.^2 + pen;             % sum of squared error (not residual)
                                % as target of "1" is exactly known!!
end
% <><><><><><><><><><><><> End secondary functions <><><><><><><><><><><><>
