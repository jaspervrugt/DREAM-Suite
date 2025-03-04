function S = ABC_binormal(x)
% Generate samples from bivariate normal distribution

% Store local variables in memory
persistent n_pairs sigma id_1 id_2

% Load the data and define local variables - only once
if isempty(n_pairs)
    % Number of x,y pairs
    n_pairs = 10;
    % Define the sigma of the bivariate distribution
    sigma = 0.01^2 * eye(2);
    % Now reorganize the x-vector so that we get bivariate mu values
    id_1 = 1 : 2 : 2 * n_pairs - 1; id_2 = 2 : 2 : 2 * n_pairs;
end

% Check whether delta has been specified
if size(x,2) == 2 * n_pairs
    % Set delta to zero
    x = [ x 0 ];
end

% Define mu
mu = [ x(id_1)' x(id_2)' ];

% Now loop over each bivariate normal distribution
for i = 1 : n_pairs
    % Draw 50 points from each bivariate normal distribution
    X(i,1:2) = mean ( mvnrnd ( mu(i,1:2) , sigma, 50 ) + ...
        normrnd( 0 , x(end) , 50 , 2 ) );
end

% Now return vector with model simulated values
S = X(:);