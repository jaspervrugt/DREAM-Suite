function [pred] = BMA_quantile(BMA,X,w,sigma,alpha,b,c);
% This function finds the quantiles of the BMA model
% assumption is that conditional pdf is a normal distribution
%
% -------------------------------------------------------------------------
% Input:
%  X            matrix with (bias-corrected) model forecasts
%  sigma        sigma value from EM/MCMC algorithm (single value or vector)
%  w            weights from EM/MCMC algorithm
%  alpha        Confidence interval
%
% Output:
%  pred         BMA prediction corresponding to requested confidence level
% -------------------------------------------------------------------------

if nargin < 6,
    b = []; c = [];
end

if nargin == 6,
    c = [];
end

switch BMA.PDF,
    
    case 'normal',
        
        % How many sigma values do we have?
        n_sigma = numel(sigma);
        
        % Check whether we have the same sigma for all models of ensemble
        if n_sigma < BMA.k,
            % Create copies of sigma
            sigma = sigma * ones(1,BMA.k);
        end
        
    case 'heteroscedastic'
        
        % Calculate measurement error of data
        sigma = abs(b*X);
        
end

% Calculate mean of std
mean_std = mean(std(X(1:end,1:BMA.k)));

% Now calculate min_X and max_X
min_X = max ( 1e-3 , min(X(1:end,1:BMA.k),[],2) - 1 * mean_std ); max_X = max(X(1:end,1:BMA.k),[],2) + 1*mean_std;

% Now loop over each observation and determine BMA quantiles
for i = 1 : size(X,1),
    
    % Define x values
    x = [ min_X(i) : (max_X(i) - min_X(i))/9999 : max_X(i) ];
    
    % Initialize pdf to be zero
    pdf = 0;
    
    % Loop over each model
    switch BMA.PDF
        
        case 'normal'
            
            % Loop over all models of ensemble
            for j = 1 : BMA.k,
                
                % Now calculate the pdf for each x value
                pdf = pdf + w(j) * normpdf(x,X(i,j),sigma(j));
                
            end
            
        case 'heteroscedastic'
            
            % Loop over all models of ensemble
            for j = 1 : BMA.k,
                
                % Now calculate the pdf for each x value
                pdf = pdf + w(j) * normpdf(x,X(i,j),sigma(i,j));
                
            end
            
        case {'gamma'}
            
            % Loop over all models of ensemble
            for j = 1 : BMA.k,
                
                mu = abs(X(i,j));                       % Derive mean of gamma distribution
                var = abs(c + b(j) * X(i,j));           % Derive variance of gamma distribution
                A = mu.^2./var; B = var./mu;            % Derive A and B of gamma distribution
                prob = gampdf(x,A,B);
                % Check for problems with inf
                idx = ( find ( prob == Inf ) ); prob(idx) = prob(min((idx+1),size(x,2)));
                pdf = pdf + w(j) * prob;       % Now calculate the pdf for each x value
                
            end
            
    end
    
    % Determine cdf
    cdf = cumsum(pdf)/sum(pdf);
    % Remove duplicate y-values!
    
    idx = find ( cdf < alpha ( 1 ) );
    if isempty(idx),
        x = [ -1 x ]; cdf = [ 0 cdf ];
    else
        idx = idx(end); x = x(idx:end); cdf = cdf ( idx:end );
    end
    
    % Same for alpha ( 2 );
    idx = find ( cdf > alpha ( 2 ) );
    if isempty(idx),
        x = [ x x(end) + 1 ]; cdf = [cdf 1];
    else
        idx = idx(1); x = x(1:idx); cdf = cdf ( 1 : idx);
    end
    
    % Now determine the appropriate quantiles
    for m = 1 : numel(alpha),
        % Store in matrix "pred" --> interpolate cdf at alpha value
        pred(i,m) = interp1(cdf,x,alpha(m));
    end
end