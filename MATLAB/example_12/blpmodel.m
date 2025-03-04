function loglik = blpmodel(p)
% p contains the parameters of the linear model and variogram function best linear predictor

% Store model input variables locally in memory
persistent np nt M H z 

% Load data and define variables - only once
if isempty(np)
    
    % Load the data
    data = load('forest.txt'); 
    % Define coordinates
    x = data(:,1); y = data(:,2); z = data(:,3); 
    % Create matrix of coordinates
    X = [x , y];  
    % Create trend vector
    M = ones(length(z),1);
    % distance matrix between all observations
    H = distmat(X);   
    % Define np
    np = 1;
    % Define nt
    nt = 5;

end

% parameter of the linear model
beta = p ( 1 : np );

% the linear model
S =  M * beta';

% parameter of the variogram
theta = p ( np + 1 : nt );   

% compute covariance structure of the data
C = cova_matern(theta(1),theta(2),theta(3),theta(4),H);    

% Compute Log likelihood
dC = det(C) + realmin;

% Calculate the L1 norm
L1 = ( z - S )'* inv(C) *( z - S );

% Calculate the log-likelihood
loglik  = -0.5 * ( log(dC) + L1 );

end