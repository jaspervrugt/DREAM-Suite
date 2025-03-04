function [LM_res_SWIG,opt_LM_SWIG,SSR_LM_SWIG,n_LM_SWIG,min90,max90,min95,max95] = par_SWIG(Bfixed,approach)
% For parallel implementation

% Load data
load SWIG_1D
% Rename content of file
data_SWIG = SWIG_1D;
% Unpack the data
n_soil = size(data_SWIG,1);
% Summarize data
for soil_type = 1:n_soil
    dat = data_SWIG{soil_type};
    Data(soil_type,1:4) = [ soil_type size(dat,1) max(dat(:,1)) max(dat(:,2)) ];
end
eta_min = [ 1e-3 1e-3 ]'; eta_max = [ inf inf ]';
% Number of parameters
d = 2;
% How many trials with LM from different starting points?
Ntrials = 10; nu = 2; %1.3; % nu = 10;
warning off
for soil_type = 1:n_soil
    % Unpack data
    dat = data_SWIG{soil_type};
    % Define plugin structure
    n_data = size(dat,1);
    % Store time and infiltration
    t_meas = dat(1:n_data,1); I_meas = dat(1:n_data,2);
    % Execute Levenberg Marquardt method
    [eta_opt,SSR_opt,K] = least_squares_haverkamp(t_meas,I_meas,approach,Bfixed,eta_min,eta_max,nu,Ntrials);
    % How many trials converge to same optimum?
    [nr,Km,Kstd,Kopt] = analyze_opt(eta_opt,SSR_opt,K);
    % Store
    LM_res_SWIG(soil_type,1:4) = [ nr Km Kstd Kopt ];
    %% Store output of least_squares
    zz = find(SSR_opt == min(SSR_opt)); zz = zz(1);
    % Least squares values of LM method
    opt_LM_SWIG(soil_type,1:d) = eta_opt(zz,1:d);
    % SSR of least squares values of LM method
    SSR_LM_SWIG(soil_type) = SSR_opt(zz);
    % Store number of data points
    n_LM_SWIG(soil_type) = n_data;
    % Compute uncertainty at optimum
    [std_opt,r90,r95] = det_uncertainty(I_meas,t_meas,Bfixed,eta_opt(zz,1:d),SSR_opt(zz),eta_min,eta_max,n_data,approach);
    % Now store 95%
    min95(soil_type,1:d) = r95(1:d,1)';
    max95(soil_type,1:d) = r95(1:d,2)';
    min90(soil_type,1:d) = r90(1:d,1)';
    max90(soil_type,1:d) = r90(1:d,2)';
end