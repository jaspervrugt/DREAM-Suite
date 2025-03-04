function logpdf = GAME_mixture_density(gmix,X,cov_type,shared)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
% GGGG   AA   M   M  EEEE  SSSS   AA   M   M  PPPP  L     I  N    N  GGGG %
% G  G  A  A  M   M  E     S     A  A  M   M  P  P  L     I  NN   N  G  G %
% G  G  A  A  MM MM  EEE   S     A  A  MM MM  P  P  L     I  N N  N  G  G %
% GGGG  AAAA  M M M  EEE - SSSS  A  A  M M M  PPPP  L     I  N  N N  GGGG %
%    G  A  A  M   M  E        S  AAAA  M   M  P     L     I  N   NN     G %
% GGGG  A  A  M   M  EEEE  SSSS  A  A  M   M  P     LLLL  I  N    N  GGGG %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Computes log density of normal mixture importance distribtion           %
%                                                                         %
% SYNOPSIS:                                                               %
%  logpdf = GAME_mixture_density(gmix,X,cov_type)                         %
% WHERE                                                                   %
%  gmix        [input] Structure trained normal mixture importance dist.  %
%  .n              # samples (= R) of X                                   %
%  .d              # dimensions (= parameters) of X                       %
%  .k              # components of mixture distribution                   %
%  .p              # parameters of mixture distribution                   %
%  .w              maximum likelihood weights of mixture components       %
%  .mu             jxd matrix of mean values each component               %
%  .Sigma          dxdxj array of covariance matrices each component      %
%  .I              integral of pdf each component before applying weight  %
%  .loglik         log-likelihood of normal mixture                       %
%  .AIC            Akaike information criterion of normal mixture         %
%  .BIC            Bayesian information criterion of normal mixture       %
%  X          [input] Nxd matrix of samples                               %
%  cov_type   [input] OPT: Covariance type                                %
%   = 1            Diagonal covariance matrix [diagonal entries only]     %
%   = 2            Full covariance matrix [= all entries are estimated]   %
%  shared     [input] OPT: shared covariance matrix Σ mixture components  %
%   = 0            Mixture components share same covariance matrix Σ      %
%   = 1            Mixture components have different covariance matrices  %
%  logpdf     [outpt] Nx1 vector logarithmic density mixture at X samples %
%                                                                         %
% © Written by Jasper A. Vrugt, Jan 2015                                  %
%   Adapted from wdensity                                                 %
% University of California Irvine                                         %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

if nargin < 4, shared = 0; end              % Each component own Σ matrix
if nargin < 3, cov_type = 2; end            % Covariance type

[mu,Sigma,w,known] = ...                    % Unpack the structure mix
    deal(gmix.mu,gmix.Sigma,gmix.w,0);
log_pr = log(w);                            % Log value component weights
[n,d] = size(X);                            % # samples, # parameters
K = gmix.k;                                 % # mixture components
ell = nan(n,K);                             % log-likelihood mixture at X
save gmix.mat gmix

for k = 1:K
    switch shared
        case 1      % k components same Σ matrix
            z = 1; if k > 1, known = 1; end
        otherwise
            z = k;  % k components own Σ matrix
    end
    if cov_type == 1 && known == 0          % diagonal Σ matrix, 
        L = sqrt(Sigma(:,:,z));             % can take sqrt, Sigma = L'*L;    
        if any(L < eps(max(L)) * d)         % check L entries
            error(['GAME_mixture_density: ' ...
                'Ill-conditioned covariance']);
        end
        logdetSigma = ...
            sum(log(Sigma(:,:,z)));         % compute log(|Sigma(:,:,z)|)
    elseif cov_type == 2 && known == 0      % full Σ matrix
        [L,f] = chol(Sigma(:,:,z));         % Cholesky fact., Sigma = L'*L
        diagL = diag(L);                    % diagonal L entries
        if (f ~= 0) || any(abs(diagL) < ... % check diagL entries
            eps(max(abs(diagL))) * size(L,1))
            error(['GAME_mixture_density: ' ...
                'Ill-conditioned covariance']);
        end
        logdetSigma = 2*sum(log(diagL));    % compute log(|Sigma(:,:,z)|)
    else
        error(['GAME_mixture_density: ' ...
            'unknown covariance type']);
    end
    Xc = bsxfun(@minus,X,mu(k,:));          % Center samples
    if cov_type == 2
        XcRinv = Xc/L;
    else
        XcRinv = bsxfun(@times,Xc,(1./L));
    end
    d_M = sum(XcRinv.^2, 2);                % nx1 vector Mahalanobis dist.
                                            % = distance measure between a
                                            % point and a distribution
    ell(1:n,k) = - .5 * d_M ...             % Loglik each sample each comp.
        - .5 * logdetSigma + log_pr(k) ...
        - .5 * d * log(2*pi) ...
        - log(gmix.I(k));
    % ell(:,j) = -.5 * d_M -.5 *logDetSigma + log_pr(j) - d*log(2*pi)/2;
    % pdf = p/norm_cnst; log(pdf) = log(p) - log(norm_cnst); subtract 
    % The last term on the right hand side accounts for bounded case when
    % jth component mixture does not integrate to w(j) * 1, j = 1,...,k
    % In unbounded case, integral, gmix.I(j) = 1 for all components, thus
    % we yield w(j)
    % pdf = w * mvnpdf(X,µ,Σ); then log(pdf) = log(w) + log(mvnpdf(X,µ,Σ))
    % Thus, for 1 sample, we add log(w); for N samples effect is N*log(w)
end

% OLD CODE: LONGER
% % for j = 1:k
% %     switch shared
% %         case 1      % k normal components have same Σ matrix
% %             if j == 1
% %                 if cov_type == 1        % diagonal Σ matrix
% %                     L = sqrt(Sigma);
% %                     if any(L < eps(max(L)) * d)
% %                         error(['GAME_sampling ERROR: weighted_density: '...
% %                             'Ill-conditioned covariance']);
% %                     end
% %                     logdetSigma = sum(log(Sigma));
% %                 elseif cov_type == 2    % full Σ matrix
% %                     [L,f] = chol(Sigma); diagL = diag(L);
% %                     if (f ~= 0) || any(abs(diagL) < ...
% %                             eps(max(abs(diagL))) * size(L,1))
% %                         error(['GAME_sampling ERROR: weighted_density: '...
% %                             'Ill-conditioned covariance']);
% %                     end
% %                     logdetSigma = 2*sum(log(diagL));
% %                 else
% %                     error(['GAME_sampling ERROR: unknown covariance' ...
% %                         ' type']);
% %                 end
% %             end
% %         otherwise   % k normal components have own Σ matrix
% %             if cov_type == 1        % diagonal covariance
% %                 L = sqrt(Sigma(:,:,j)); % a vector
% %                 if any(L < eps(max(L))*d)
% %                     error(['GAME_sampling ERROR: weighted_density: ' ...
% %                         'Ill-conditioned covariance']);
% %                 end
% %                 logdetSigma = sum(log(Sigma(:,:,j)));
% %             elseif cov_type == 2    % full Σ matrix
% %                 [L,f] = chol(Sigma(:,:,j)); diagL = diag(L);
% %                 if (f ~= 0) || any(abs(diagL) < ...
% %                         eps(max(abs(diagL))) * size(L,1))
% %                     error(['GAME_sampling ERROR: weighted_density: ' ...
% %                         'Ill-conditioned covariance']);
% %                 end
% %                 logdetSigma = 2*sum(log(diagL));
% %             else
% %                 error(['GAME_sampling ERROR: unknown covariance' ...
% %                     ' type']);
% %             end
% %     end
% %     Xc = bsxfun(@minus,X,mu(j,:));      % Center samples
% %     if cov_type == 2
% %         XcRinv = Xc/L;
% %     else
% %         XcRinv = bsxfun(@times,Xc,(1./L));
% %     end
% %     d_M = sum(XcRinv.^2, 2);            % nx1 vector Mahalanobis distances
% %     % = distance measure between a
% %     % point and a distribution
% %     ell(1:n,j) = -.5 * d_M ...          % Log-liklhd each sample each comp.
% %         - .5 * logdetSigma ...
% %         + log_pr(j) ...
% %         - d*log(2*pi)/2 ...
% %         - log(gmix.I(j));
% %     % ell(:,j) = -.5 * d_M -.5 *logDetSigma + log_pr(j) - d*log(2*pi)/2;
% %     % The last term on the right hand side accounts for bounded case
% %     % In unbounded case, integral, gmix.I(j) = 1 for all components, j in k
% %
% % end

ell_max = max(ell,[],2);                % Max. log-likelihood each sample 
p = exp(bsxfun(@minus,ell,ell_max));    % minus ell_max to avoid underflow
pdf = sum(p,2);                         % pdf(i) = Σ_j w_j P(x_i|Θ_j) / ...
                                        %     exp(ell_max(i)) 
% p = bsxfun(@rdivide, p, pdf);         % normalize densities
logpdf = log(pdf) + ell_max;            % Calculate log-df

end
