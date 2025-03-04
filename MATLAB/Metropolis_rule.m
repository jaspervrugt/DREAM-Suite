function [accept,id_acc] = Metropolis_rule(DREAMPar,Meas_info,...
    log_alfa_sn,logLPR_Xp,logLPR_X,options)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Metropolis rule for acceptance or rejection of candidate points       %%
%%                                                                       %%
%% SYNOPSIS: [accept,id_acc] = Metropolis_rule(DREAMPar,Meas_info,...    %%
%%                log_alfa_sn,logLPR_Xp,logLPR_X,options)                %%
%%                                                                       %%
%% © Written by Jasper A. Vrugt, June 2006 & modified June 2024          %%
%% Los Alamos National Laboratory                                        %%
%% University of California Irvine                                       %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

logPR_Xp = logLPR_Xp(:,1);      % Log-prior density of candidate points
logL_Xp = logLPR_Xp(:,2);       % Log-likelihood of candidate points
logPR_X = logLPR_X(:,1);        % Log-prior density of current chain states
logL_X = logLPR_X(:,2);         % Log-likelihood of current chain states

a_L = exp(logL_Xp - logL_X);    % Likelihood ratio
Z = rand(DREAMPar.N,1);         % Draw standard uniform labels
accept = nan(DREAMPar.N,1);     % Initialize vector with accept id

switch options.ABC
    case {'no'}     % NO ABC --> REGULAR MCMC WITH PRIOR AND LIKELIHOOD 
                    % --> POSSIBLY DIAGNOSTIC BAYES (no/yes)
        if strcmp(options.DB,'no')      % NO DIAGNOSTIC BAYES
            a_PR = exp(logPR_Xp - logPR_X); % Prior ratio
            a_S = exp(log_alfa_sn);         % Snooker correction
            alfa = a_S .* a_L .* a_PR;      % Product of ratios
            accept = alfa >= Z;             % Accept if alfa ≥ Z
        elseif strcmp(options.DB,'yes') % DIAGNOSTIC BAYES
            for z = 1:DREAMPar.N            % Check pairwise logPR_Xp and logPR_X
                if (logPR_Xp(z) >= logPR_X(z))  % If proposal better than current chain state
                    if logPR_X(z) < 0               % And if current chain state outside epsilon
                        accept(z) = 1;                  % Then, accept proposal
                    else                            % If current chain state inside epsilon
                        accept(z) = a_L(z) >= Z(z);     % Then, accept with Metropolis probability
                    end
                else                            % If proposal worse than current chain state
                    if logPR_Xp(z) < 0              % And if proposal outside epsilon
                        accept(z) = 0;                  % Then, reject proposal
                    else                            % If proposal inside epsilon
                        accept(z) = a_L(z) >= Z(z);     % Accept with Metropolis probability 
                    end
                end
            end
        end
    case {'yes'}    % YES TO ABC OR LIMITS OF ACCEPTABILITY
                    % --> SUMMARY METRICS AS FITNESS MEASURE (LIKELIHOOD)
        switch DREAMPar.lik
            case 21 % Turner approach
                a_PR = exp(logPR_Xp - logPR_X); % Prior ratio
                a_S = exp(log_alfa_sn);         % Snooker correction
                alfa = a_S .* a_L .* a_PR;      % Product of ratios
                accept = alfa >= Z;             % Accept if alfa ≥ Z
            case 22 % Sadegh and Vrugt (2014) --> epsilon is one value
                    % accept = (logL_Xp ≤ logL_X) | (logL_Xp ≤ options.epsilon);
                accept = (logL_Xp >= logL_X) ...    % epsilon multiple values
                            | (logL_Xp >= 0); 
            case 23 % Limits of acceptability --> Vrugt and Beven, 2018
                accept = (logL_Xp >= logL_X) ...    
                            | (logL_Xp == Meas_info.n_S); 
        end
    otherwise
        error(['DREAM_PACKAGE ERROR:Metropolis_rule: I do not know ',...
            'this option']);
end

id_acc = find(accept > 0);  % Index proposals to accept [= row numbers Xp]

end
