function e = EnKF_Lorenz63(pars,plugin)
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
% Lorenz 1963 model. After:
% Lorenz, E.N. (1963), Deterministic non-periodic flow, Journal of
%   Atmospheric Sciences, 20, 130–141, 1963
% 
% SODA method. After:
% Vrugt, J.A., C.G.H. Diks, H.V. Gupta, W. Bouten, and J.M. Verstraten
%   (2005), Improved treatment of uncertainty in hydrologic modeling: 
%   Combining the strengths of global optimization and data assimilation, 
%   Water Resources Research, 41, W01017, doi:10.1029/2004WR003059
% © Jasper A. Vrugt, University of Amsterdam, 2003
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

persistent m N R C
if isempty(m)
    m = plugin.m; N = plugin.N; R = plugin.R; C = plugin.C;
end

% Data assimilation of the Lorenz model using the Ensemble Kalman Filter.
% The parameters of this model are assumed known and their values are
% provided as input to this function, along with other variables that are
% stored in the structure "plugin". This variable allows user to port
% information from example_26 to this function.
%
% This example is equivalent to case study 1 in this published article and
% uses the Lorenz model - developed by Edward Lorenz:
% Lorenz, E.N. (1963), Deterministic non-periodic flow, Journal of
%   Atmospheric Sciences, 20, 130–141, 1963
%
% The approach used in this case study is identifical to the SODA approach
% advocated by Vrugt et al. (2005):
% Vrugt, J.A., C.G.H. Diks, H.V. Gupta, W. Bouten, and J.M. Verstraten (2005),
%   Improved treatment of uncertainty in hydrologic modeling: Combining the
%   strengths of global optimization and data assimilation, Water Resources
%   Research, 41, W01017, doi:10.1029/2004WR003059
%
% I recommend that you check the following lecture to understand in more 
% detail the concept, theory, and implementation of data assimilation: 
% https://www.youtube.com/watch?v=l5cb4ONIB-c&feature=youtu.be&disable_polymer=true
% including methods such as SODA. 
% This and other lectures can also be found on 
% my website: faculty.sites.uci.edu/jasper/teaching
%
% In short: SODA uses an inner layer of state estimation with the ENKF 
% embedded within an outer layer of MCMC simulation with DREAM. This 
% approach I originally developed to get a much better understanding of
% epistemic (model structural) errors - and to minimize the impact of input
% data and epistemic errors on the optimized parameter values. What is
% more, detailed analysis of the state updates can reveal new insights 
% into model strucxtural errors

% Unpack parameters of Lorenz model and parameter AR_rho EnKF model error 
p = num2cell(pars); [sigma,rho,beta,AR_rho] = deal(p{:});
% Lorenz equations using specified parameter values
Lorenz_1963 = @(t,u) [ - sigma*u(1) + sigma*u(2); ...
                         rho*u(1) - u(2) - u(1)*u(3); ...
                       - beta*u(3) + u(1)*u(2) ];

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%               Now execute ensemble Kalman Filter
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% Set the first mean simulated value equal to observed value
E(1,1:m) = plugin.Y(1:m,1);
% Draw m x N matrix of measurement errors from N_m(0,R)
Meas_err = mvnrnd(zeros(1,m),R,N)';
% m x N matrix of initial states of ensemble members: min of zero
X_a = bsxfun(@plus,plugin.Y(1:m,1),Meas_err);
% m x N matrix of state forecasts of ensemble members
X_f = nan(m,N);
% Draw m x N matrix with model errors
Mod_err = mvnrnd(zeros(1,m),C,N)';
% 
% Now run the model one-observation ahead
for t = 1 : plugin.NT - 1     % plugin.NT - 1 --> we start at t = 0
    % Generate state forecast for each ensemble member
    for ii = 1 : N
        % Run the model from one time step to next for the nth state vector
        [~,state_forecast] = ode45(Lorenz_1963,[(t-1) * plugin.dt t * ...
            plugin.dt],X_a(1:m,ii),plugin.options);
        % Store the simulated state forecast at final time
        X_f(1:m,ii) = state_forecast(end,1:m);
        % Compute m x n matrix of model errors at time t using AR(1) scheme
        Mod_err(1:m,ii) = AR_rho * Mod_err(1:m,ii) + ...
            sqrt(1 - AR_rho^2) * mvnrnd(zeros(1,m),C)';
    end
    % m x N matrix of state forecasts
    X_f = X_f + Mod_err;
    % Compute m x 1 vector of mean state forecast and store in matrix E
    e = mean(X_f,2); E(t+1,1:m) = e';
    % Draw m x N matrix of measurement errors at time t
    Meas_err = mvnrnd(zeros(1,m),R,N)';
    % Generate m x N matrix of measured data replicates at time t + 1 
    D = bsxfun(@plus,plugin.Y(1:m,t+1),Meas_err);
    %<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    % Intermittent block (not needed for filter but for plotting)
    %<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    
    % Now determine 95% ranges of forecasted state
    % T = sort(X_f,2); a1 = ceil(N*(alpha/2)); a2 = N - a1;
    % 2.5/97.5 percentile (column m+1 to 3*m of S)
    % S(t+1,m+1:2*m) = T(1:m,a1)'; S(t+1,2*m+1:3*m) = T(1:m,a2)';
    
    %<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    % End intermittent block (not needed for filter but for plotting)
    %<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    
    % Step 3: Subtract mean state from each vector (ensemble member) of Xf
    A = bsxfun(@minus,X_f,e);
    % Step 4: Compute m x m sample covariance matrix C of state forecasts
    C = 1/(N-1) * (A*A'); % which is similar to C = cov(A');
    % Step 5: Compute the m x m Kalman gain matrix 
    % K = C * inv(C + R); 
    K = C/(C+R);
    % Step 6: Compute m x N matrix of analysis states
    X_a = X_f + K * ( D - X_f );
end

% Return mean state forecast at each time as a single vector
e = E(:);

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%                   End of ensemble Kalman Filter
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

end