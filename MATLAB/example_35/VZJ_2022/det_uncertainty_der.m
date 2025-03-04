function [std_opt,r90,r95] = det_uncertainty_der(I_meas,t_meas,eta_opt,SSR_opt,eta_min,eta_max,n_data,approach)
% Determine uncertainty at optimum

% Approach one or two
switch approach
    case 1 % Infiltration form
        plugin.t = t_meas; 
    case 2 % Time form
        plugin.I = I_meas;
end

% Now compute nxp Jacobian matrix at current iterate
Jf = jac_Haverkamp_der_anal(eta_opt,plugin,approach);
% Number of parameters
d = 3;
% Determine uncertainty
C = SSR_opt/(n_data-d) * inv(Jf'*Jf + eye(d,d)*1e-10);
% And confidence intervals
std_opt = sqrt(diag(C));
% And 95% intervals
I90 = std_opt * tinv(0.95,n_data-d);
% Now ranges of all three
r90 = [ eta_opt' - I90 eta_opt' + I90 ];
% Now check whether in bound - lower 
r90(1:d,1) = max(r90(1:d,1),eta_min);
% Now check whether in bound - upper 
r90(1:d,2) = min(r90(1:d,2),eta_max);
% And 95% intervals
I95 = std_opt * tinv(0.975,n_data-d);
% Now ranges of all three
r95 = [ eta_opt' - I95 eta_opt' + I95 ];
% Now check whether in bound - lower 
r95(1:d,1) = max(r95(1:d,1),eta_min);
% Now check whether in bound - upper 
r95(1:d,2) = min(r95(1:d,2),eta_max);