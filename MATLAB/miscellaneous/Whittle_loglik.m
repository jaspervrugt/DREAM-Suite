function loglik = Whittle_loglik(fx,Meas_info)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Function computes Whittle's log-likelihood function using spectral densities of    %%
%% measurements and model output                                                      %%
%%                                                                                    %%
%% SYNOPSIS: loglik = Whittle_loglik(fx,Meas_info)                                    %%
%%  where                                                                             %%
%%   fx        [input]  REQUIRED: n x 1 vector of simulated model output              %%
%%   Meas_info [input]  REQUIRED: Measurement structure with measured values          %%
%%   loglik    [outpt]  Whittle quasi-log-likelihood                                  %%
%%                                                                                    %%
%% For more information please check                                                  %%
%%  Whittle, P. (1957), Curve and Periodogram Smoothing, Journal of the Royal         %%
%%      Statistical Society, Ser. B, 19, 38-63                                        %%
%%  Whittle, P. (1962), Gaussian Estimation in Stationary Time Series, Bulletin of    %%
%%      the International Statistical Institute, 39, 105-129                          %%
%%                                                                                    %%
%% © Written by Jasper A. Vrugt, Feb. 2007                                            %%
%% Los Alamos National Laboratory 			        	                              %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

e = Meas_info.Y - fx;                       % n x 1 residual vector
[phi,s2_pe] = armcov(e,1); phi = phi(2);    % AR(1)-coef. phi and part. variance, s2_pe  
% e_(t) = phi*e_(t-1) + pe_(t)  where pe ~ N(0,s2_pe) is partial residual 

%% The next lines can be computed up front to speed up calculations
n_half = floor((Meas_info.n - 1)/2);
% Compute periodogram of measured data (spectral density)
per_meas = abs(fft(Meas_info.Y)).^(2 / ( 2 * pi * Meas_info.n));
% Now take a certain set of points from per_obs
per_meas = per_meas ( 2: ( n_half + 1 ) );
% Now work with autoregressive error model
id = (1:n_half)' * (( 2 * pi ) / Meas_info.n);
% Now calculate sine and cosine of idx
sin_id = sin(id); cos_id = cos(id);
%% 

% Compute periodogram of simulated values (spectral density)
per_sim = abs(fft(fx)).^(2 / ( 2 * pi * Meas_info.n));
% Now take a certain set of points from per_sim
per_sim = per_sim(2:(n_half + 1));
% Now define I_ar and R_ar
I_ar = phi * sin_id; R_ar = phi * cos_id;
% Now calculate f_ar
f_ar = (1 - R_ar).^2 + I_ar.^2;
% Calculate f_spec
f_spec = (1./f_ar).* s2_pe / (2 * pi);
% Add spectral density AR(1) to spectral density of model
per_tot = per_sim + f_spec;
% Compute ratio of spectral densities of data and joint model
y_f = per_meas./per_tot;
% Which elements of per_tot are larger than zero?
idx = per_tot > 0;

% Whittle's log-likelihood
loglik = - sum(log(per_tot(idx))) - sum(y_f(idx));
%loglik = - log(per_tot(idx)) - y_f(idx);

end