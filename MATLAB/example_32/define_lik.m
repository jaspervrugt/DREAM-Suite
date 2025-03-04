function [DREAMPar,Par_info,idx_vpar,fpar,fname] = define_lik(DREAMPar,...
    Par_info,par_names,fpar_mod,parmin_mod,parmax_mod)
% This function returns likelihood variables

switch DREAMPar.lik
    case 13 %% Normal distribution
        % index:           1   2    3   4
        % parname:        s0  s1  phi1 phi2
        fpar_nuis   =  [ 1e-4  0    0   0 ];
        parmin_nuis =  [   0   0    0   0 ];
        parmax_nuis =  [   1   1    1   1 ];
        fname = 'Normal'; idx_nuis = [1 2 3 4];
    case 14 %% GL
        % index:           1    2    3   4   5   6    7    8    9  10  11
        % parname:       std0 std1 beta xi  mu1 phi1 phi2 phi3 phi4 K lambda
        fpar_nuis   =  [  0.1   0    0   1   0   0    0    0    0   0   1  ];
        parmin_nuis =  [   0    0   -1  0.1  0   0    0    0    0   0  0.1 ];
        parmax_nuis =  [   1    1    1  10  100  1    1    1    1   1   1  ];
        fname = 'GL';
    case 17 %% SST
        % index:          1   2    3   4    5    6
        % parname:       s0  s1   nu   xi  phi1  phi2
        fpar_nuis   =  [ 1e-4 0  1e10  1    0    0 ];
        parmin_nuis =  [  0   0   2   0.1   0    0 ];
        parmax_nuis =  [  1   1  100  10    1    1 ];
        fname = 'SST'; idx_nuis = [3 4 5];
    case 44 %% GL_plus
        % index:          1   2    3   4    5    6
        % parname:       s0  s1  beta xi phi1  phi2
        %fpar_nuis   =  [ 0.01 0    0   1    0    0 ];
        fpar_nuis   =  [ 1e-4 0    0   1    0    0 ];
        parmin_nuis =  [  0   0   -1  0.1   0    0 ];
        parmax_nuis =  [  1   1    1   10   1    1 ];
        fname = 'GL_plus'; idx_nuis = [3 4 5];
    case 45 %% SGT
        % index:          1    2    3   4   5    6    7
        % parname:        s0   s1 labda p   q  phi1 phi2
        fpar_nuis   =  [ 1e-4  0    0   2  1e10  0    0 ];
        parmin_nuis =  [  0    0   -1  0.5   2   0    0 ];
        parmax_nuis =  [  1    1    1  100  100  1    1 ];
        fname = 'SGT'; idx_nuis = [3 4 5 6];
    otherwise
        [idx_nuis,fname,fpar_nuis,parmin_nuis,parmax_nuis] = deal([]);
end

nuis_names = nuis_var_names(DREAMPar);          % Extract names nuis. variables
names = [par_names,nuis_names];                 % Combine with parameter names
n_modpar = numel(fpar_mod);                     % Number of model parameters
idx_vpar = [1:n_modpar , idx_nuis + n_modpar ]; % Model parameter + nuisance variable selection
fpar = [fpar_mod fpar_nuis];                    % Merge default values parameters and nuisance variables
parmin = [parmin_mod parmin_nuis];              % Merge the parameter names and nuisance variables   
parmax = [parmax_mod parmax_nuis];   
Par_info.min = parmin(idx_vpar);                % Min/max values of parameter selection
Par_info.max = parmax(idx_vpar); 
Par_info.names = names(idx_vpar);               % Names of parameters and nuisance variables

end