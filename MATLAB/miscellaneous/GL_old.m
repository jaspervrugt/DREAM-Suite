function [logL,std_e,a,fa,ExpY,SimY] = GL_old(iflag,par,ModY,ObsY,Nrep)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Generalized likelihood function regression models with correlated, heteroscedastic %%
%%  and non-Gaussian errors                                                           %%
%%                                                                                    %%
%%  Regression model:                                                                 %%
%%  1. Expected values: simulated using external function, optional bias correction   %%
%%  2. Residual standard deviations: modeled as a linear function of expected values  %%
%%  3. Residual higher-order moments: modeled by Skew Exponential Power distribution  %%
%%  4. Residual correlations: modeled by autoregressive time-series (up to order 4)   %%
%%                                                                                    %%
%%  SYNOPSIS: [logL,std_e,a,fa,ExpY,SimY] = GL_old(iflag,par,ModY,ObsY,Nrep)          %%
%%  where                                                                             %%
%%   iflag      [input] Estimation ('est') or simulation ('sim')                      %%
%%   par        [input] Row vector deterministic & nuisance variables (fixed/estimtd) %%
%%    s0:  nuisvar(1) Intercept of linear heteroscedastic model                       %%
%%    s1:  nuisvar(2) Slope of linear heteroscedastic model                           %%
%%    ba:  nuisvar(3) Kurtosis (-1: uniform, 0: normal; 1: Laplace)                   %%
%%    xi:  nuisvar(4) Skewness (1: symmetric; <1: negative skew; >1: positive skew)   %%
%%    mu1: nuisvar(5) Bias correction parameter                                       %%
%%    fi1: nuisvar(6) First-order AR coefficient (0,1)                                %%
%%    fi2: nuisvar(7) Second-order AR coefficient (0,1)                               %%
%%    fi3: nuisvar(8) Third-order AR coefficient (0,1: check book)                    %%
%%    fi4: nuisvar(9) Fourth-order AR coefficient (0,1: check book)                   %%
%%    lbd: nuisvar(10) Box-Cox transformation parameter (skewness)                    %%
%%    K:   nuisvar(11) Box-Cox transformation parameter (heteroscedasticity)          %%
%%   ModY       [input] nx1 vector of modeled response variables deterministic model  %%
%%   ObsY       [input[ nx1 vector of observed response values                        %%
%%   Nrep       [input] Number of replicates - resampling                             %%
%%   logL       [outpt] nx1 vector of log-likelihood values                           %%
%%   std_e      [outpt] nx1 vector standard deviation of raw residuals                %%
%%   a          [outpt] nx1 vector stndrdzd decorrelated residuals ~SEP(0,1,beta,xi)  %%
%%   fa         [outpt] nx1 vector of SEP densities of a                              %%
%%   ExpY       [outpt] nx1 vector bias-corrected modeled response variables          %%
%%   SimY       [outpt] nxN matrix of replicate simulations (for Bayes_pdf)           %%
%%                                                                                    %%
%%  Reference:                                                                        %%
%%   Schoups, G., and J. A. Vrugt (2010), A formal likelihood function for parameter  %%
%%       and predictive inference of hydrologic models with correlated,               %%
%%       heteroscedastic, and non-Gaussian errors, Water Resources Research, 46,      %%
%%       W10531, doi:10.1029/2009WR008933                                             %%
%%                                                                                    %%
%%  Notes: 1. This is a conditional likelihood function: y_(-1) and y_0 assumed zero  %%
%%         2. Measurement error computed from simulated data                          %%
%%         3. This function is OBSOLETE/ERRONEOUS, please use GL_plus                 %%
%%                                                                                    %%
%% © Written by Gerrit Schoups - modified by Jasper A. Vrugt                          %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

if nargin < 5, Nrep = 1; end
SimY = [];

nuisvar = par(end-10:end);  % Separate deterministic & statistical parameters
s0 = nuisvar(1);            % intercept of linear heteroscedastic model
s1 = nuisvar(2);            % slope of linear heteroscedastic model
ba = nuisvar(3);            % kurtosis (-1: uniform, 0: normal; 1: Laplace)
xi = nuisvar(4);            % skewness (1: symmetric; <1: negative skew; >1: positive skew)
mu1 = nuisvar(5);           % parameter for bias correction
fi1 = nuisvar(6);           % first-order autoregressive coefficient
fi2 = nuisvar(7);           % second-order autoregressive coefficient
fi3 = nuisvar(8);           % third-order autoregressive coefficient
fi4 = nuisvar(9);           % fourth-order autoregressive coefficient
K = nuisvar(10);            % Box-Cox transformation parameter (skewness)
lba = nuisvar(11);          % Box-Cox transformation parameter (heteroscedasticity)

% EXPECTED VALUES
n = size(ModY,1);                       % # entries of simulated record
ExpY = ModY.*min(10,exp(mu1.*ModY));    % bias corrected response variable

% STANDARD DEVIATION
std_e = max(1e-6,s1.*ExpY + s0);    % Measurement error std. response variable

% KURTOSIS AND SKEWNESS
A1 = gamma(3*(1+ba)/2);
A2 = gamma((1+ba)/2);
Cb = (A1/A2)^(1/(1+ba));
Wb = sqrt(A1)/((1+ba)*(A2^(1.5)));
M1 = gamma(1+ba)/sqrt(A1*A2);
M2 = 1; mu_xi = M1*(xi-1/xi);
sig_xi = sqrt((M2-M1^2)*(xi^2 + 1/xi^2) + 2*M1^2 - M2);

% CORRELATION
phi_p = [1 -fi1 -fi2 -fi3 -fi4];        % Coefficients of AR polynomial

% ESTIMATION
e = ((ObsY+K).^lba - 1)./lba ...        % BC transformed residuals
    - ((ExpY+K).^lba - 1)./lba;
a = filter(phi_p,1,e);                  % i.i.d. errors (a)
a = a./std_e;                           % Studentized a values
a_xi = (mu_xi + sig_xi.*a)./...         % SEP transformed a values
    (xi.^sign(mu_xi + sig_xi.*a));  
fa = (2*sig_xi/(xi + 1/xi)) * Wb .* ... % Density of a
    exp(-Cb.*(abs(a_xi).^(2/(1+ba))));
logL = n.*log(Wb*2*sig_xi/(xi+1/xi))... % Log-likelihood
    - sum(log(std_e)) - Cb ...
    .*(sum(abs(a_xi).^(2./(1+ba))));
logL = logL + (lba-1) * ...             % SEP log-likelihood
    sum(log(ObsY + K));

% SIMULATION: Generate response variables (= SimY)
switch iflag
    case 'sim'
        % rand('seed',sum(100*clock));    % Random number generator
        % Generate N i.i.d. errors (a) from SEP(0,1,xi,beta)
        % Step 1 - Generate N random variates from gamma distribution with 
        % shape parameter 1/p and scale parameter 1
        p = 2/(1+ba); SimY = nan(n,Nrep);
        for ii = 1:Nrep
            grnd = gamrnd(1/p,ones(n,1));
            % Step 2: n random signs (+1 or -1) with equal probability
            signrnd = sign(rand(n,1)-0.5);
            % Step 3: n random variates from EP(0,1,beta)
            EP_rnd = signrnd.*(abs(grnd).^(1/p)).*...
                sqrt(gamma(1/p))./sqrt(gamma(3/p));
            % Step 4: n random signs (+1 or -1) with probability 1-w and w
            w = xi/(xi+1/xi);
            signrndw = sign(rand(n,1)-w);
            % Step 5: n random variates from SEP(mu_xi,sig_xi,xi,beta)
            SEP_rnd = -signrndw.*abs(EP_rnd)./(xi.^signrndw);
            % Step 6: Normalize for n random variates from SEP(0,1,xi,beta)
            a = (SEP_rnd-mu_xi)./sig_xi;
            % nx1 vector of raw BC residuals
            e = filter(1,phi_p,std_e.*a);
            % Replicate simulated values
            SimY(1:n,ii) = (lba.*(((ExpY+K).^lba - 1)./lba + e) + 1).^(1/lba) - K;
            % Assumption: E[g^-1(Y)] = g^-1(E[Y]) where g = Box-Cox transformation
        end
end

end