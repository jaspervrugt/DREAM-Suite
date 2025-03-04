function loglik = mvt_lik(x)
% Multivariate student t distribution with df degrees of freedom and
% correlation matrix R

% Store local variables in memory
persistent d df R logNumer_A logSqrtDetC logDenom

% Define covariance matrix
if isempty(d)
    
    % Define dimensionality of target distribution
    d = size(x,2);
    % Degrees of freedom of student target distribution?
    df = 60;
    % Define mean of t-distribution
    mu = zeros ( 1 , d );
    % Construct the d x d covariance matrix
    A = 0.5 * eye( d ) + 0.5 * ones( d );
    % Rescale to variance-covariance matrix of interest
    for i = 1 : d
        for j = 1 : d
            C(i,j) = A(i,j) * sqrt(i * j);
        end
    end
    % Standardize C to correlation if necessary.  This does NOT standardize X.
    s = sqrt(diag(C));
    if (any(s~=1))
        C = C ./ (s * s');
    end
    % Make sure C is a valid covariance matrix
    R = cholcov(C,0);
    % Define logNumer_A
    logNumer_A = gammaln(( df + d )/2 ) - gammaln(df/2);
    % Define logSqrtDetC
    logSqrtDetC = sum(log(diag(R)));
    % Define logDenom
    logDenom = logSqrtDetC + (d/2) * log(df*pi);

end
disp(R)
% Normalize x with R
Z = x / R;
disp(Z)
% Define logNumer_B
logNumer_B = -((df + d)/2) .* log(1+sum(Z.^2, 2)./df);
% Calculate log-density of multivariate t-distribution, for t_pdf(x, mu, and df)
loglik = logNumer_A + logNumer_B - logDenom;

end