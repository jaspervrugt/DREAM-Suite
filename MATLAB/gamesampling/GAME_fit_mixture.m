function gmix = GAME_fit_mixture(X,logP,metric,K,lower,upper)
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
% This function determines the optimal Gaussian mixture distribution      %
% for Nxd matrix of posterior samples, X, derived from (e)DREAM Package.  %
% We try mixture distributions consisting of k = 1,...,K components. Each %
% of these mixtures is trained using the expectation maximization method. %
% Then, the optimal mixture distribution is determined using model        %
% selection criteria such as the AIC, BIC and variance ratio, 'var', of   %
% the target (= posterior) and Gaussian mixture PDF.                      %
%                                                                         %
% The normal mixture PDF of x \in X [= 1xd vector] is written as          %
%          PDF(x) = Σ_{j=1}^{k} w_{j} f_{N_{d}}(x,µ_{j},Σ_{j})            %
%                                                                         %
% where f_{N_{d}}(x,µ_{j},Σ_{j}) is the d-variate normal pdf with mean µ  %
% [= 1xd vector] and dxd covariance matrix Σ and the weights w_{1},...,   %
% w_{k} must sum to 1. The weights, w_{j}, mean vector, µ_{j}, and dxd    %
% covariance matrix, Σ_{j}, of each mixture component, j = 1,...,k, are   %
% determined using Expectation-Maximization. The default option is that   %
% each mixture component has its own and full covariance matrix. Then,    %
% the EM algorithm has to train p1 = d - 1 weights, namely, w_{1},...,    %
% w_{k-1}, p2 = k*d mean entries, µ_{1},...,µ_{k}, and p3 = k*d*(d+1)/2   %
% elements of the covariance matrices, Σ_{1},...,Σ_{k}. The total # of    %
% parameters of the normal mixture equals p = p1+p2+p3. In MATLAB, we     %
% yield that f_{N_{d}}(x,µ_{j},Σ_{j}) = mvnpdf(x,µ_{j},Σ_{j})             %
%                                                                         %
% SYNOPSIS:                                                               %
%  gmix = GAME_fit_mixture(X,logP,metric,K,lower,upper)                   %
% WHERE                                                                   %
%  X           [input] Rxd matrix posterior samples DREAM                 %
%  logP        [input] Rx1 vector of logarithmic values posterior density %
%  metric      [input] string (name) metric used for mixture selection    %
%   = 'bic'        Bayesian information criterion        DEFault          %
%   = 'var'        Variance reduction                                     %
%  K           [input] Maximum # components mixture      DEF: 5           %
%                      [= importance] distribution                        %
%  lower       [input] 1xd vector lower bounds parameters                 %
%  upper       [input] 1xd vector upper bounds parameters                 %
%  gmix        [outpt] Structure trained normal mixture importance dist.  %
%  .n              # samples (= R) of X                                   %
%  .d              # dimensions (= parameters) of X                       %
%  .K              # components of mixture distribution                   %
%  .p              # parameters of mixture distribution                   %
%  .w              maximum likelihood weights of mixture components       %
%  .mu             jxd matrix of mean values each component               %
%  .Sigma          dxdxj array of covariance matrices each component      %
%  .I              integral of pdf each component before applying weight  %
%  .loglik         log-likelihood of normal mixture                       %
%  .AIC            Akaike information criterion of normal mixture         %
%  .BIC            Bayesian information criterion of normal mixture       %
%                                                                         %
% © Written by Jasper A. Vrugt, Jan 2015                                  %
% University of California Irvine                                         %
%                                                                         %
% Version 1:    June 2012       Initial setup and definition              %
% Version 1.1:  Jan. 2015       Major overhaul of code                    %
% Version 1.2:  Feb. 2015       VAR method: min. variance ratio q and p   %
%				Oct. 2015       Add bounded Gaussian mixture distribution %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

var_ratios = inf(1,K);      % variance of the ratio between q and p
bic = inf(1,K);             % initialize mixture BIC
% Fit a mixture distribution through this sample --> # modes determined 
% by maximizing log-density. This mixture of Gaussians will act as 
% importance distribution. We use a "full" covariance matrix - that is
% with correlated dimensions and each component havings its own properties
for k = 1:K
    % Try fitting j components - if that does not work go with j-1 mixtures
    try
        % Fit a Gaussian mixture using built-in MATLAB function
        gmix(k) = emgm(X',k,lower,upper);                           %#ok
        logQ = GAME_mixture_density(gmix(k),X);
        var_ratios(k) = var(logP - logQ);
        bic(k) = gmix(k).BIC;
    catch
        % Display that fitting crashed
        fprintf(['GAME_fit_mixture: Cannot determine a suitable ' ...
            '%d-component mixture distribution'],k) 
        fprintf(['GAME_fit_mixture: We continue with a %d component ' ...
             'Gaussian mixture distribution\n'],k+1);
    end
end
% Sort BIC in ascending order to find best mixture [= lowest BIC]
[~,j1] = sort(bic(1:K)); j1 = j1(1); 
[~,j2] = sort(var_ratios(1:K)); j2 = j2(1);

% Now return mixture based on bic or var metric
switch metric
    case {'bic'}
        gmix = gmix(j1);
    case {'var'}
        gmix = gmix(j2);
end

if (j1 == j2)
    fprintf(['GAME_fit_mixture: BIC and VAR metric lead to same ' ...
        'optimal mixture distribution\n']);
else
    fprintf(['GAME_fit_mixture: BIC and VAR metric do not yield ' ...
        'same optimal mixture distribution. Check results with ' ...
        'GAMEoptions.metric = ''var'' \n']);
end

end

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
% SECONDARY FUNCTIONS
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
function [gmix,m_loglik] = emgm(X,k,lower,upper)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%EMGM Training of normal mixture model using the Expectation-Maximization %
% algorithm                                                               %
%                                                                         %
% SYNOPSIS:                                                               %
%  [gmix,loglik] = emgm(X,k,lower,upper)                                  %
% WHERE                                                                   %
%  X           [input] Rxd matrix posterior samples DREAM                 %
%  k           [input] scalar with # components normal mixture dist.      %
%  lower       [input] 1xd vector lower bounds parameters                 %
%  upper       [input] 1xd vector upper bounds parameters                 %
%  gmix        [outpt] Structure trained normal mixture importance dist.  %
%  .d              # dimensions (= parameters) of X                       %
%  .n              # samples (= R) of X                                   %
%  .k              # components of mixture distribution                   %
%  .p              # parameters of mixture distribution                   %
%  .w              maximum likelihood weights of mixture components       %
%  .mu             jxd matrix of mean values each component               %
%  .Sigma          dxdxj array of covariance matrices each component      %
%  .I              integral of each component pdf before applying weight  %
%  .loglik         log-likelihood of normal mixture                       %
%  .AIC            Akaike information criterion of normal mixture         %
%  .BIC            Bayesian information criterion of normal mixture       %
%  m_loglik    [outpt] 1xit vector of mean log-likelihoods mixture dist.  %
%                                                                         %
% © Modified by Jasper A. Vrugt, Jan 2015                                 %
%   based on code of Michael Chen (sth4nth@gmail.com)                     %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

%
% X: d x n data matrix
%   init: j (1 x 1) or label (1 x n, 1<=label(i)<=j) or center (d x j)
% Written by Michael Chen (sth4nth@gmail.com).

% initialization
fprintf(['GAME_fit_mixture:emgm: Fitting of Gaussian ' ...
    'mixture distribution ... \n']);

R = initialization(X,k);
[~,label(1,:)] = max(R,[],2);
R = R(:,unique(label));

tol = 1e-10;
%maxiter = 500;
maxiter = 2500;
m_loglik = -inf(1,maxiter);
converged = false;
it = 1;
while ~converged && it < maxiter        % Expectation-Maximization
    it = it + 1;
    gmix = maximization(X,R);           % Maximization Step
    [R,m_loglik(it)] = ...              % Expectation Step
        expectation(X,gmix);     
    [~,label(:)] = max(R,[],2);
    u = unique(label);                  % non-empty components
    if size(R,2) ~= size(u,2)
        R = R(:,u);                     % remove empty components
    else
        converged = m_loglik(it) - ...  % EM algorithm converged?
            m_loglik(it - 1) < ...
            tol * abs(m_loglik(it));
    end

end
m_loglik = m_loglik(2:it);
if converged
    fprintf(['GAME_fit_mixture:emgm: ' ...
        'Parameters of mixture ' ...
        'distribution have converged ' ...
        'in %d steps.\n'],it - 1);
else
    fprintf(['GAME_fit_mixture:emgm: ' ...
        'Parameters of mixture ' ...
        'distribution did not converge ' ...
        'in %d steps.\n'],maxiter);
end

[gmix.d,gmix.n] = size(X);              % # data points 
gmix.k = k;                             % # components of mixture
gmix.p = k-1 ...                        % # parameters mixture: k-1 weights 
    + k*gmix.d ...                      %   & d*j mean values 
    + k*gmix.d*(gmix.d+1)/2;            %   & j*(d+1)*d/2 entries of Sigma
gmix.mu = gmix.mu';                     % mean of each mixture component 
gmix.loglik = ...                       % total log-likelihood
    gmix.n * m_loglik(it - 1); 
gmix.AIC = -2 * gmix.loglik ...         % Akaike's information criterion
    + 2 * gmix.p;                       %  → want to minimize
gmix.BIC = -2 * gmix.loglik ...         % Bayes information criterion
    + gmix.p * log(gmix.n);             %  → want to minimize

% Compute integral of components' pdf in bounded case (lower/upper bound)
% We use the function qsimvn developed by Alan Genz. This function uses 
% an algorithm presented in the following paper:                     
%   Genz, Alan (1992), Numerical Computation of Multivariate Normal       
%       Probabilities, Journal of Computational and Graphical Statistics, 
%       1, pp. 141-149. WSU Math, PO Box 643113, Pullman, WA 99164-3113   
%       Email: AlanGenz@wsu.edu

gmix.I = nan(1,k);
switch gmix.d > 1
    case 1
        [low,up] = deal(nan(1,gmix.d)); precision = 1e4;
        for k = 1:gmix.k
            for j = 1:gmix.d
                low(j) = lower(j) - gmix.mu(k,j);
            end
            for j = 1:gmix.d
                up(j) = upper(j) - gmix.mu(k,j);
            end
            gmix.I(k) = qsimvn(precision,gmix.Sigma(:,:,k),low',up');
        end
    otherwise
        for k = 1:gmix.k
            gmix.I(k) = mvncdf(upper(:),gmix.mu(k,:), ...
                gmix.Sigma(:,:,k)) - mvncdf(lower(:), ...
                gmix.mu(k,:),gmix.Sigma(:,:,k));
        end
end

end

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
% SECONDARY FUNCTIONS
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
function R = initialization(X,init)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%INITIALIZATION Initialization of the mixture                             %
%                                                                         %
% SYNOPSIS:                                                               %
%  R = initialization(X,init)                                             %
%                                                                         %
% © Written by Michael Chen (sth4nth@gmail.com)                           %
%   Modified by Jasper A. Vrugt                                           %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

[d,n] = size(X);
if isstruct(init)                               % initialize with a model
    R = expectation(X,init);
elseif length(init) == 1                        % random initialization
    k = init;
    idx = randsample(n,k);
    m = X(:,idx);
    [~,label] = max(bsxfun(@minus,m'*X, ...
        dot(m,m,1)'/2),[],1);
    [u,~,label] = unique(label);
    while k ~= numel(u)
        idx = randsample(n,k);
        m = X(:,idx);
        [~,label] = max(bsxfun(@minus,m'*X, ...
            dot(m,m,1)'/2),[],1);
        [u,~,label] = unique(label);
    end
    R = full(sparse(1:n,label,1,n,k,n));
elseif size(init,1) == 1 && size(init,2) == n   % initialize with labels
    label = init;
    k = max(label);
    R = full(sparse(1:n,label,1,n,k,n));
elseif size(init,1) == d                        % initialize with centers
    k = size(init,2);
    m = init;
    [~,label] = max(bsxfun(@minus,m'*X, ...
        dot(m,m,1)'/2),[],1);
    R = full(sparse(1:n,label,1,n,k,n));
else
    error(['GAME_fit_mixture:emgm:initialization: ' ...
        'Invalid initialization.']);
end

end

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
% SECONDARY FUNCTIONS
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
function gmix = maximization(X,R)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%MAXIMIZATION Maximization step of EM algorithm                           %
%                                                                         %
% SYNOPSIS:                                                               %
%  gmix = maximization(X,R)                                               %
%                                                                         %
% © Written by Michael Chen (sth4nth@gmail.com)                           %
%   Modified by Jasper A. Vrugt                                           %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

[d,n] = size(X);                            % # parameters, # samples of X
k = size(R,2);                              % # components mixture
nk = sum(R,1);

w = nk/n;                                   % component weights; Σw = 1
mu = bsxfun(@times,X*R,1./nk);              % components mean
Sigma = zeros(d,d,k);                       % components covariance mtrx, Σ 
sqrtR = sqrt(R);
for j = 1:k
    Xo = bsxfun(@minus,X,mu(:,j));
    Xo = bsxfun(@times,Xo,sqrtR(:,j)');
    Sigma(:,:,j) = Xo*Xo'/nk(j);
    Sigma(:,:,j) = Sigma(:,:,j) ...
        + 1e-6 * eye(d);                    % numerical stability (P.D.)
end
gmix = struct('mu',mu,'Sigma',Sigma,'w',w); % Return Gaussian mixture

end

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
% SECONDARY FUNCTIONS
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
function [R,loglik] = expectation(X,gmix)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%EXPECTATION Expectation step of EM algorithm                             %
%                                                                         %
% SYNOPSIS:                                                               %
%  [R,loglik] = expectation(X,gmix)                                       %
%                                                                         %
% © Written by Michael Chen (sth4nth@gmail.com)                           %
%   Modified by Jasper A. Vrugt                                           %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

[~,n] = size(X);                    % # samples of X
k = size(gmix.mu,2);                % # components of mixture distribution
logrho = zeros(n,k);                % initialize logRho
for j = 1:k                         % compute pdf of each mixture comp
    logrho(1:n,j) = loggausspdf( ...
        X,gmix.mu(:,j), ...
        gmix.Sigma(:,:,j));
end
logrho = bsxfun(@plus,logrho, ...   % add log(weights) to each component 
    log(gmix.w));
T = logsumexp(logrho,2);            % take sum of component pdfs
loglik = sum(T)/n;                  % log-density of mixture distribution
logR = bsxfun(@minus,logrho,T);     % average log-density of all X
R = exp(logR);                      % scaled density all X

end

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
% SECONDARY FUNCTIONS
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
function logpdf = loggausspdf(X,mu,Sigma)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%LOGGAUSSPDF Logarithmic density of Nd(µ,Σ) at N samples of X             %
% = log(mvnpdf(X,mu,Sigma))'                                              %
% but with protection against numerical underflow                         %
%                                                                         %
% SYNOPSIS:                                                               %
%  logpdf = loggausspdf(X,mu,Sigma)                                       %
%                                                                         %
% © Written by Michael Chen (sth4nth@gmail.com)                           %
%   Modified by Jasper A. Vrugt                                           %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

d = size(X,1); 
X = bsxfun(@minus,X,mu);
[U,p] = chol(Sigma);
if p ~= 0
    error(['GAME_fit_mixture:emgm:' ...
        'expectation:loggausspdf: ' ...
        'Sigma is not positive definite.']);
end
Q = U'\X;
q = dot(Q,Q,1);                         % quadratic term (M distance)
c = d*log(2*pi)+2*sum(log(diag(U)));    % normalization constant
logpdf = -(c+q)/2;

end
