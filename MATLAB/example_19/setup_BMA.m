function [DREAMPar,Par_info,D_bc,A,B] = setup_BMA(DREAMPar, ...
    Par_info,D,y,VAR)
% Setup BMA mixture model estimation
% I recommend using instead the MODELAVG toolbox. This provides a much more
% advanced implementation of the BMA model, including confidence/prediction
% intervals, scoring rules and other metrics of the BMA distribution
% forecast

[n,K] = size(D);                            % # forecasts, # ensmble mmbers
adjust = 'true';                            % Linear bias correction
[D_bc,A,B] = ComputeAB(D,y,n,K,adjust);     % Do bias correction

par_name = cell(1,K);
for z = 1:K, par_name{z} = strcat('\beta_{',num2str(z),'}'); end

switch VAR
    case '1'    % 1: common constant variance
        DREAMPar.d = K + 1;
        Par_info.max = [ ones(1,K) 2 * std(y) ];
        par_name{K+1} = '\sigma';
        for z = 1:K, par_name{z} = strcat('\beta_{',num2str(z),'}'); end
    case '2'    % 2: individual constant variance
        DREAMPar.d = 2 * K;
        Par_info.max = [ ones(1,K) 2 * std(y)*ones(1,K) ];
        for z = 1:K, par_name{K+z} = strcat('\sigma_{',num2str(z),'}'); end
    case {'3'}  % 3: common non-constant variance
        DREAMPar.d = K + 1;
        Par_info.max = [ ones(1,K) 2 ];
        par_name{K+1} = 'c';
    case {'4'}  % 4: individual non-constant variance
        DREAMPar.d = 2 * K;
        Par_info.max = [ ones(1,K) 2*ones(1,K) ];
        for z = 1:K, par_name{K+z} = strcat('c_{',num2str(z),'}'); end
    otherwise
        error(['AMALGAM:setup_BMA:Unknown variance option; ' ...
            'choose between ''1''/''2''/''3''/''4'' ']);
end
Par_info.names = par_name;                  % Parameter names
Par_info.unit_simplex = 'yes';              % Weights on unit Simplex
Par_info.min = zeros(1,DREAMPar.d);       % Min. values BMA weights/vars

end
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
% Secondary functions
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% 1: ComputeAB
function [D_bc,A,B] = ComputeAB(D,y,n,K,adjust)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%ComputeAB: Function performs linear bias correction if indicated
%
%  SYNOPSIS
%   [D_bc,A,B] = ComputeAB(D,y,n,K)
%   [D_bc,A,B] = ComputeAB(D,y,n,K,adjust)
%  where
%   D       [input] NxK matrix of ensemble forecasts of K models
%   y       [input] Nx1 vector of verifying data
%   n       [input] Number of forecasts
%   K       [input] Number of ensemble members
%   adjust  [input] OPT: linear bias correction ('true') or not
%   D_bc    [outpt] NxK matrix with bias-corrected ensemble forecasts
%   A       [outpt] Intercepts of linear bias-correction
%   B       [outpt] Slopes of linear bias-correction
%
% Written by Jasper A. Vrugt
% Los Alamos National Laboratory
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

if nargin < 5, adjust = 'true'; end

if strcmpi(adjust,'true')
    % Bias correction terms (linear regression)
    B = zeros(1,K); A = zeros(1,K);
    for k = 1:K
        T = cov(D(:,k),y); T = T(1,2);      % Temporary variables
        B(k) = T/var(D(:,k));               % Slope linear regrs. func.
        A(k) = mean(y) - B(k)*mean(D(:,k)); % Intercept linear regrs. func.
    end
else
    B = ones(1,K); A = zeros(1,K);
end

D_bc = nan(size(D));                        % Init. bias-correctd forecasts
for k = 1:K
    D_bc(1:n,k) = B(k) * D(1:n,k) + A(k);   % Bias-corrected forecasts
end

end
