function varargout = Calc_likelihood(Xp,FXp,DREAMPar,Par_info,Meas_info,...
    Lik_info,options,MAP_info)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% This function computes the log-likelihood of the candidate points Xp               %%
%%                                                                                    %%
%% SYNOPSIS: logPR_Xp = Calc_likelihood(Xp,FXp,DREAMPar,Par_info,Meas_info,...        %%
%%               Lik_info,options,MAP_info)                                           %%
%%  where                                                                             %%
%%   Xp        [input] N x d matrix of candidate points                               %%
%%   FXp       [input] n x N matrix with simulated values candidate points            %%
%%   DREAMPar  [input] Algorithmic structure with variables MCMC method               %%
%%   Par_info  [input] Parameter structure                                            %%
%%   Meas_info [input] Measurement data structure                                     %%
%%   Lik_info  [input] Structure with information likelihood function                 %%
%%   options   [input] Structure with computational settings/options                  %%
%%   MAP_info  [input] OPT: Information MAP solution for sandwich estimator           %%
%%   loglik_Xp [outpt] Nx1 vector with logarithmic value of likelihood of Xp          %%
%%   E         [outpt] nxN matrix with residuals of Xp                                %%
%%   std_m     [outpt] nxN matrix with measurement errors of Xp for Kalman proposal   %%
%%   eps_n     [outpt] nxN matrix of standardized partial residuals of Xp             %%
%%   f_eps_n   [outpt] nxN matrix of density standardized partial residuals of Xp     %%
%%   ell       [outpt] nxN matrix of log-likelihoods of Xp                            %%
%%                                                                                    %%
%% © Written by Jasper A. Vrugt, Feb 2007                                             %%
%% Los Alamos National Laboratory                                                     %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

if nargin < 8 % Sandwich correction not relevant
    MAP_info = struct;
end

n = Meas_info.n;        % Relabel number of measurements
n_S = Meas_info.n_S;    % Relabel number of summary metrics
N = size(Xp,1);         % # candidate vectors? [can differ from DREAMPar.N]
loglik_Xp = nan(N,1);   % Initialize log-likelihood candidate points
E = [];                 % Initialize nxN residual
if n > 0                % FXp is a simulation
    E = bsxfun(@minus,Meas_info.Y,FXp);     
end
ell = nan(max(max(n_S,n),1),N);

% User specified measurement error?
if isfield(Meas_info,'Sigma')
    std_e = Meas_info.Sigma;
    % Now make n copies of std_e
    if numel(std_e) == 1, std_e = std_e * ones(n,1); end
    % Initialize std_m
    std_m = repmat(std_e,1,N);
else
    std_m = nan(Meas_info.n,N);
end

for ii = 1 : N  % Loop over each candidate point
    
    switch DREAMPar.lik
        case 1  % Func_name returns a likelihood
            ell(:,ii) = log(FXp(:,ii)); loglik_Xp(ii,1) = sum(ell(:,ii));
            
        case 2  % Func_name returns a log-likelihood
            ell(:,ii) = FXp(:,ii); 
            loglik_Xp(ii,1) = sum(ell(:,ii)); 

        case 11 % Measurement error integrated out
            loglik_Xp(ii,1) = -(n/2)*log(sum(abs(E(1:n,ii)).^2));
            % --> This is unnormalized posterior density [= uniform prior]
            % --> not eliptical, that means log(sum) not equal to sum(log)
            ell(1,ii) = loglik_Xp(ii,1);
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % NOTES
            % 1: 1/std(e) * f(e_n,0,1) = f(e,0,std(e));
            % 2: RMSE = sqrt( sum(abs(E(:,ii)).^2) / n )
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        case 12 % Normal likelihood with homos/heteroscedastic meas. error
            e_n = E(1:n,ii)./std_e;  % Standardize residuals
            ell(1:n,ii) = - (1/2)*log(2*pi) - log(std_e) - 1/2 * e_n.^2;
            loglik_Xp(ii,1) = sum(ell(1:n,ii),'omitnan');
            % loglik_Xp(ii,1) = -(n/2)*log(2*pi) - ...
            %     sum(log(std_e),'omitnan') - 1/2 * sum(e_n.^2);
            % == sum(log(normpdf(E,0,std_e)))
            %disp(Xp)
            %disp(loglik_Xp)
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % NOTES
            % 1. Similar to following matrix implementation
            %    loglik_Xp(ii,1) = - (n/2)*log(2*pi) - sum(log(std_e)) ...
            %                    - 1/2*E(:,ii)'*diag(1./std_e.^2)*E(:,ii);
            % 2. 1/std(e) * f(e_n,0,1) = f(e,0,std(e));
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        case {13,14,16,17,44,45} % NL/GL/LAPL/SSTL/GL+/SGTL function
            par = Lik_info.fpar;                            % fixed pars
            par(Lik_info.id_vpar) = Xp(ii,1:DREAMPar.d);    % variable pars
            nuisvar = par(Lik_info.id_nuisvar);             %#ok iso nuis. vars 
            eval(Lik_info.stringL);                         % Likelihood
            ell(:,ii) = loglik;                             % store all logliks
            std_m(:,ii) = std_e;                            % for Kalman proposals
            loglik_Xp(ii,1) = sum(loglik,'omitnan');        % loglik is vector

        case 15 % Whittle quasi-likelihood function (Whittle 1957, 1962)
            loglik_Xp(ii,1) = Whittle_loglik(FXp(1:n,ii),Meas_info);
            ell(1,ii) = loglik_Xp(ii,1);                    % store all logliks
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % NOTES
            % 1: only log-likelihoods for spectral densities > 0 included
            % 2: that is why we cannot return individual log-likelihoods 
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        case 21 % Approximate Bayesian Computation (some estimate options.epsilon along)
            phi1 = options.rho(FXp(1:n_S,ii) , Meas_info.S(:)) + ...
                normrnd(0,options.epsilon);
            ell(:,ii) = - (1/2)*log(2*pi) - log(options.epsilon) - ...
                    1/2 * (phi1./options.epsilon).^2;
            loglik_Xp(ii,1) = sum(ell(:,ii));
           % loglik_Xp(ii,1) = - ( n_S / 2) * log(2 * pi) - ...
           %     sum ( log ( options.epsilon ) ) - ...
           %     1/2 * sum ( ( phi1./options.epsilon ).^2 );
            
        case 22 % Approximate Bayesian Computation (alternative to continuous kernel)
            % Now calculate log-density (not a true log-density! - epsilon has multiple values)
% %             loglik_Xp(ii,1) = min(options.epsilon - ...
% %                 abs(options.rho(FXp(1:n_S,ii) , Meas_info.S(:))));
            loglik_Xp(ii,1) = min(options.epsilon - ...
                options.rho(FXp(1:n_S,ii) , Meas_info.S(:)));
            ell(1,ii) = loglik_Xp(ii,1);
            
        case 23 % Approximate Bayesian Computation (limits of acceptability)
            % log-likelihood = number of points that satisfy LOAs
            ell(:,ii) = abs(FXp(1:n_S,ii) - Meas_info.S(:)) <= ...
                options.epsilon;
            loglik_Xp(ii,1) = sum(ell(:,ii));
            
        case 31 % GLUE: Table 1: Option a in Beven and Freer, 2001
            loglik_Xp(ii,1) = - DREAMPar.GLUE * log ( var(E(:,ii)) );
            
        case 32 % GLUE: Table 1: Option b in Beven and Freer, 2001
            a = var(E(:,ii)); b = var(Meas_info.Y) ;
            if ( a < b )
                loglik_Xp(ii,1) = DREAMPar.GLUE * log ( 1 - a / b );
            else
                % Note: Careful with this implementation! Guiding DREAM
                % to solution of a < b --> so this is never in posterior!
                loglik_Xp(ii,1) = -(a-b) * 1000 * DREAMPar.GLUE;
            end
            
        case 33 % GLUE: Table 1: Option c in Beven and Freer, 2001
            loglik_Xp(ii,1) = - DREAMPar.GLUE * ( var(E(:,ii)) );
            
        case 34 % GLUE: Beven and Binley, 1992 (last option, Page 284)
            loglik_Xp(ii,1) = - log ( sum ( abs ( E(:,ii) ) ) );
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        case 52 % Matrix implementation of Gaussian likelihood: GLS form
            if isfield(Meas_info,'C')
                loglik_Xp(ii,1) = - ( n/2 ) * log(2 * pi) - ...
                    1/2 * log(Meas_info.detC) - ...
                    1/2 * E(:,ii)' * Meas_info.invC * E(:,ii);
                % This formulation needs to be adapted so that
                % C = sigma2_e * V and user specifies V
            else % sigma is estimated along with the model parameters
                % form only if diagonal covariance matrix - zeros elsewhere
                loglik_Xp(ii,1) = - ( n/2 ) * log(2 * pi) - ...
                    sum(log(std_e)) - ...
                    1/2 * E(:,ii)' * diag(1./std_e.^2) * E(:,ii);
                % Equivalent to matrix form
                % Sigma2_e = diag(std_e);
                % logL_Xp(ii,1) = - ( n/2 ) * log(2 * pi) - ...
                %     log(det(Sigma2_e)) - 1/2 * E(:,ii)' * ...
                %     inv(Sigma2_e) * E(:,ii);
            end
            
        case 61 % Laplace likelihood with variable learning rate
            lambda = Xp(ii,end-2);
            b = std_e/sqrt(2); % Check LAPrnd
            f_lambda = 1./(2*b).^lambda .* exp(-abs(E(:,ii)./b)).^lambda;
            % Normalizing constant, lambda > 0, prior [0,2]?
            Z_lambda = 1./(lambda * (2*b).^(lambda-1));
            % Compute log-likelihood
            % logL_Xp(ii,1) = sum(log(f_lambda/Z_lambda));
            loglik_Xp(ii,1) = sum(log(f_lambda)) - sum(log(Z_lambda));
            
        case 62 % Normal likelihood with variable learning rate
            lambda = Xp(ii,end-2);
            f_lambda = 1./(std_e*sqrt(2*pi)).^lambda .* ...
                exp(-1/2*(E(:,ii)./std_e).^2).^lambda;
            % Normalizing constant, lambda > 0, prior [0,2]?
            Z_lambda = 1./(sqrt(lambda) * (std_e * sqrt(2*pi)).^(lambda-1));
            % Compute log-likelihood
            % logL_Xp(ii,1) = sum(log(f_lambda/Z_lambda));
            loglik_Xp(ii,1) = sum(log(f_lambda)) - sum(log(Z_lambda));
            
    end
    
end

% Sandwich adjustment of log-likelihood
if isfield(MAP_info,'map')
    [a,loglik_adj] = deal(nan(N,1));
    % Check whether we work with normalized values [0-1] or not
    switch Par_info.norm 
        case 1
            Xp_n = X_normalize(Xp,Par_info);
        otherwise
            Xp_n = Xp;
    end
    for ii = 1:N
        dx = Xp_n(ii,1:DREAMPar.d) - MAP_info.map; % 1 x p vector
        a(ii) = (dx*MAP_info.Godambe*dx') / (dx*MAP_info.Fisher*dx');
        loglik_adj(ii,1) = a(ii) * (loglik_Xp(ii,1) - MAP_info.loglik);
    end
    if any(a < 0)
        fprintf(['DREAM_PACKAGE WARNING: Negative value of a(ii) in ',...
            'Sandwich correction\n']);
    end
    loglik_Xp = loglik_adj;
end
% End of sandwich adjustment

% Now define return arguments
switch nargout
    case 1  % Return only the log-likelihood
        varargout(1) = {loglik_Xp};
    case 3  % Return default variables
        varargout(1:3) = {loglik_Xp,E,std_m};
    case 5  % Return normalized residuals + density of them as well
        % Note: normpdf(e,0,std(e)) == 1/std_e * normpdf(e_n,0,1);
        switch DREAMPar.lik
            case {1,2,21,22,23} % User returns likelihood/log-likelihood directly
                % Or ABC method - MAP not defined for this
                e_n = []; f_e_n = [];
            case {11,15,31,32,33,34} % Normal likelihood with sigma integrated out
                % Whitte Likelihood, Informal likelihoods of Beven etc.
                std_e = std(E(1:n,1)); % Get std from residuals
                e_n = E(1:n,1)/std_e; f_e_n = normpdf(e_n,0,1);
            case 12 % Normal likelihood with homos/heteroscedastic meas. error
                % Obsolete with likelihood 13
                f_e_n = normpdf(e_n,0,1);
            case {13,14,16,17,44,45}
                % do nothing; eps_n and f_eps_n should be known
            case 52
                if isfield(Meas_info,'C')
                    % Need to implement weight matrix, e_n = W*E(1:n,1)
                    W = sqrtm(Meas_info.invC); e_n = W*E(1:n,1);
                else
                    e_n = E(1:n,1)./std_e;
                end
                f_e_n = normpdf(e_n,0,1);
            otherwise % All other likelihoods: normalize residuals and
                e_n = E(1:n,1)./std_e; f_e_n = normpdf(e_n,0,1);
        end
        switch DREAMPar.lik
            case {13,14,16,17,44,45} % Possibly treatment residual corr.
                varargout(1:5) = {loglik_Xp,E,std_m,eps_n,f_eps_n}; %#ok 
            otherwise
                varargout(1:5) = {loglik_Xp,E,std_m,e_n,f_e_n};
        end
    case 6  % Return log-likelihood vector as well
        eps_n = []; f_eps_n = [];
        varargout(1:6) = {loglik_Xp,E,std_m,eps_n,f_eps_n,ell};
end

end