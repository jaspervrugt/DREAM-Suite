function [mean_crps,crps,num_nan] = CRP_score(fcst,obs,type_cdf,method)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% The CRPS is a quadratic measure of the difference between the forecast cdf and     %%
%% the empirical cdf of the observation.                                              %%
%%                                                                                    %%
%% SYNOPSIS: [mean_crps,crps,num_nan] = CRP_score(fcst,obs);                          %%
%%           [mean_crps,crps,num_nan] = CRP_score(fcst,obs,type_cdf);                 %%
%%           [mean_crps,crps,num_nan] = CRP_score(fcst,obs,type_cdf,method);          %%
%%  where                                                                             %%
%%   fcst      [input]  REQUIRED: n x m matrix of ensemble forecasts                  %%
%%   obs       [input]  REQUIRED: n x 1 vector of measured data (aka observations)    %%
%%   type_cdf  [input]  OPTIONAL: Empirical cumulative distribution function          %%
%%      Name    Equation                Description                                   %%
%%     'ecdf'   i/n                     Linear interpolation (MATLAB default)         %%
%%     'weib'   i/(n+1)                 Unbiased exceedance probabilities             %%
%%     'med'    (i-0.3175)/(n+0.365)    Median exceedance probabilities               %%
%%     'apl'    (i-0.35)/n              Used with PWMs                                %%
%%     'blom'   (i-0.375)/(n+0.25)      Unbiased normal quantiles                     %%
%%     'cunn'   (i-0.4)/(n+0.2)         Approx. quantile-unbiased                     %%
%%     'grin'   (i-0.44)/(n+0.12)       Optimised for Gumbel distribution             %%
%%     'hazen'  (i-0.5)/n               A traditional choice                          %%
%%     'bern'   (i-0.3)/(n + 0.4)       Benard and Bos-Levenbach, 1953                %%
%%     [p1 p2]  (i-p1)/(n+p2)           Two element vector defining custom            %%
%%   method    [input]  OPTIONAL: Which computation method; 1/2/3 (default: 3)        %%
%%   mean_crps [output] Mean of non-missing CRPS values                               %%
%%   crps_val  [output] n x 1 vector with CRPS values                                 %%
%%   num_nan   [output] Number of missing values of CRPS                              %%
%%                                                                                    %%
%% Reference:                                                                         %% 
%%   J.E. Matheson, and R.L. Winkler (1976), Scoring rules for continuous             %%
%%     probability distributions, Management Science, 22, 1087–1096.                  %%
%%   J.R. Stedinger, R.M. Vogel, and E. Foufoula-Georgiou (1995), Frequency analysis  %%
%%     of extreme events" in Maidment, D. (ed.) Handbook of Hydrology.                %%
%%     New York: McGraw-Hill.                                                         %%
%%                                                                                    %%
%% Example: fcst = normrnd(0,1,1000,1000);                                            %%
%%  	    obs = rand(1000,1); 						                              %%
%%          [meanCRPS,crps_values,num_nan] = crps_jav(fcst,obs,'weib');  	          %%
%% Note: the CRPS is positively-oriented!! (larger is better)                         %%
%%                                                                                    %%
%% (c) Written by Jasper A. Vrugt, July 2021                                          %%
%% University of California Irvine 						                               %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% Check input arguments
if nargin < 4
    method = 4;
end
if nargin < 3 || isempty(type_cdf)
    type_cdf = 'ecdf';
end
type_cdf = lower(type_cdf);     

%% Plotting position
if ischar(type_cdf)
    switch type_cdf
        case 'ecdf'
            p1 = 0; p2 = 0;
        case 'weib' % Weibull
            p1 = 0; p2 = 1;
        case 'med'  % Median
            p1 = 0.3175; p2 = 0.365;
        case 'apl'
            p1 = 0.35; p2 = 0;
        case 'blom'
            p1 = 0.375; p2 = 0.25;
        case 'cunn' % Cunnane
            p1 = 0.4; p2 = 0.2;
        case 'grin' % Gringorten
            p1 = 0.44; p2 = 0.12;
        case 'hazen'
            p1 = 0.5; p2 = 0;
        case 'bern'
            p1 = 0.3; p2 = 0.4;
        otherwise
            error('CRPS_JAV:Wrong_Empirical_CDF','Unknown empirical CDF');
    end
else
    if numel(type_cdf)~=2
        error('CRPS_JAV:Wrong_Empirical_CDF',...
            'Empirical CDF should have 2 elements, p1 and p2');
    end
    if max(type_cdf) > 1 || min(type_cdf) < 0
        error('CRPS_JAV:Wrong_Empirical_CDF',...
            'Parameters p1 and p2 of empirical CDF must be between zero and one');
    end
    p1 = type_cdf(1); p2 = type_cdf(2);
end

%% Ensemble size
[n,m] = size(fcst);         % Determine the size of the forecast matrix
% n measurement times
% m ensemble members
if size(obs,1) ~= n
    error('CRPS_JAV:WrongDimension',...
        'The length of the observation vector does not match number of rows of forecast matrix');
end
if size(obs,2) ~= 1
    error('CRPS_JAV:WrongDimension','The observation vector should have one column only');
end

crps = nan(n,1);            % initialize crps values
ys = sort(fcst,2);          % Sort entries in all rows of fcst in increasing order

%% Note: the CRPS is negatively-oriented in each computation method
%%       this is reversed into a positively oriented score at bottom of code
switch method %% Different computation methods (3rd is fastest)

    case 1
        % Compute probabilty
        p = ((1:m)-p1)./(m+p2);
        r1 = p.^2;  %ith probability squared
        r2 = (1.0-p).^2; %
        % Loop through all number of observations
        for t = 1:n
            % check if there are any missing values
            missingFcst = any(isnan(fcst(t,:)));
            missingObs = isnan(obs(t));
            if  ~missingFcst && ~missingObs
                ind = find(fcst(t,:)<=obs(t),1,'last');
                if ~isempty(ind) % i.e. obs(i) >  fcst(i,1)
                    crpsLeft = 0;
                    if ind>1 % left of the observation
                        fcstLeft = fcst(t,1:ind);
                        dxLeft = diff(fcstLeft);
                        pLeft = r1(1:ind-1);
                        crpsLeft = pLeft*dxLeft';
                    end
                    if obs(t) < fcst(t,end) % right of the observation
                        fcstRight = fcst(t,ind+1:end);
                        dxRight = diff(fcstRight);
                        if isempty(dxRight)
                            crpsRight = 0;
                        else
                            pRight = r2(ind+1:end-1);
                            crpsRight = pRight*dxRight';
                        end
                        % when the cdf crosses the observation left part
                        crpsCentreLeft = r1(ind).*(obs(t)-fcst(t,ind));
                        % right part
                        crpsCentreRight = r2(ind).*(fcst(t,ind+1)-obs(t));
                        crps(t) = crpsLeft + crpsRight + ...
                            crpsCentreLeft + crpsCentreRight;
                    else  % if observation lies right of all members of 
                          % forecast (ie. obs>fcst(end))
                        crps_right_outside = 1.0^2*(obs(t)-fcst(t,end));
                        crps(t)  =  crps_right_outside + crpsLeft;
                    end
                else  % observation lies left of the all member of 
                      % forecast (ie. obs(i) < fcst(i,1))
                    dxRight = diff(fcst(t,:));
                    pRight = r2(1:end-1);
                    crpsRight = pRight*dxRight';
                    crps_left_outside = 1.0^2*(fcst(t,1)-obs(t));
                    crps(t)  =  crps_left_outside + crpsRight;
                end
            end % not missing
        end % all forecasts

    case 2 %% Method above but vectorized
        lg = 0<1;                   % zeroth value to skip first column
        %% Incorporate the plotting position
        p = ((0:m)'-p1)./(m + p2);  % Compute probabilities: column for inner product
        r1 = p.^2;                  % Squared probabilities
        r2 = (1-p).^2;              % (1 - probabilities)^2
        i1 = 1:m-1; i2 = 1+i1;      % Integers for alfa and beta vectors
        %% Do end members, i = 0 and i = m outside measurement loop
        [alfa0m,beta0m] = deal(nan(n,1));
        % Only needs to be done for i = 0 as i = m gives same result
        idx0 = obs(1:n) <= ys(1:n,m); % (changed < to <=)
        alfa0m(idx0) = 0; alfa0m(~idx0) = obs(~idx0) - ys(~idx0,m);
        idx0 = obs(1:n) <= ys(1:n,1); % (changed < to <=)
        beta0m(idx0) = ys(idx0,1) - obs(idx0); beta0m(~idx0) = 0;
        % --> Could initialize alfa0m/beta0m at zero so as to remove two = zero statements
        %% Loop over each entry of measurement vector
        for t = 1:n
            [alfa,beta] = deal(nan(1,m+1));                         % initialize alfa and beta
            alfa([1 m+1]) = alfa0m(t); beta([1 m+1]) = beta0m(t);   % assign alfa/beta: i = 0 and i = m
            ii1 = obs(t) <= ys(t,i1);                               % first condition (changed < to <=)
            ii2 = (ys(t,i1) < obs(t)) & ( obs(t) <= ys(t,i2));      % second condition (changed < to <=)
            ii3 = ys(t,i2) < obs(t);                                % third condition
            alfa([lg ii1]) = 0;                                     % alfa: i=2:m; 1st cond.
            alfa([lg ii2]) = [ 0, obs(t) - ys(t,ii2) ];             % alfa: i=2:m; 2nd cond.
            alfa([lg ii3]) = [ 0, ys(t,i2(ii3)) - ys(t,i1(ii3)) ];  % alfa: i=2:m; 3rd cond.
            beta([lg ii1]) = [ 0, ys(t,i2(ii1)) - ys(t,i1(ii1)) ];  % beta: i=2:m; 1st cond.
            beta([lg ii2]) = [ 0, ys(t,i2(ii2)) - obs(t) ];         % beta: i=2:m; 2nd cond.
            beta([lg ii3]) = 0;                                     % beta: i=2:m; 3rd cond.
            crps(t) = alfa * r1 + beta * r2;                  % CRPS score as inner product
        end

    case 3 %% Approach B1: For empirical CDF!
        for t = 1:n
            crps(t) = 2/m^2 * sum( (ys(t,1:m) - obs(t)) ...
                .* ( m * ( obs(t) < ys(t,1:m) ) - (1:m) + 1/2 ) );
        end

    case 4 %% Approach B2: For empirical CDF! Vectorized approach of for loop above
        crps = 2/m^2 * sum ( ( ys(1:n,1:m) - obs(1:n,1) ) .* ...
            ( m * ( obs(1:n,1) < ys(1:n,1:m) ) - (1:m) + 1/2 ) , 2);

    case 5 %% Quantile formulation - see also BMA_crps
        % Can implement later
        
end

crps = -crps;                       % From negative to positive orientation (larger is better)   
mean_crps = mean(crps,'omitnan');   % Compute mean of CRPS
num_nan = sum(isnan(crps));         % Return the number of nan values