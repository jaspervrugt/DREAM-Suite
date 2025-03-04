function [tot_unc,fx_mod] = Bayes_pdf(Up,Ufx,idU,Nr,...
    RMSE_map,DREAMPar,Meas_info,Lik_info,p_alfa)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Computes prediction uncertainty of Bayes posterior for posterior parameter samples %%
%%                                                                                    %%
%% SYNOPSIS:                                                                          %%
%%  [tot_unc,fx_mod] = Bayes_pdf(Up,Ufx,idU,Nr,RMSE_map,DREAMPar,...                  %%
%%                         Meas_info,Lik_info,p_alfa)                                 %%
%% WHERE                                                                              %%
%%  Up          [input] Mxd matrix of M unique posterior parameter vectors            %%
%%  Ufx         [input] Mxn matrix of corresponding simulions                         %%
%%  idU         [input] Vector with indices for unique parameter vectors              %%
%%  Nr          [input] Vector with # replicates for each unique parameter vector     %%
%%  RMSE_map    [input] Root Mean Square Error of map solution                        %%
%%  DREAMPar    [input] Structure with algorithmic settings of MCMC algorithm         %%
%%  Meas_info   [input] Structure with measurement information (for inference)        %%
%%  Lik_info    [input] Structure with information about likelihood function          %%
%%  p_alfa      [input] Significance level (e.g: 0.01, 0.05 and 0.1)                  %%
%%  tot_unc     [outpt] Matrix with low/upper alfa prediction interval total          %%
%%  fx_mod      [outpt] Matrix with simulations of total uncertainty                  %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                                                                                    %%
%% The different componments/algorithms of DREAM Package have been described in       %%
%%   Vrugt, J.A., R. de Punder, and P. Grünwald, A sandwich with water:               %%
%%       Bayesian/Frequentist uncertainty quantification under model                  %%
%%       misspecification, Submitted to Water Resources Research, May 2024            %%
%%       https://essopenarchive.org/users/597576/articles/937008-a-sandwich-with-...  %%
%%       water-bayesian-frequentist-uncertainty-quantification-under-model-...        %%
%%       misspecification                                                             %%
%%   Vrugt, J.A. (2024), Distribution-Based Model Evaluation and Diagnostics:         %%
%%       Elicitability, Propriety, and Scoring Rules for Hydrograph Functionals,      %%
%%       Water Resources Research, 60, e2023WR036710,                                 %%
%%       https://doi.org/10.1029/2023WR036710                                         %%
%%   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of %%
%%       distribution-adaptive likelihood functions: Generalized and universal        %%
%%       likelihood functions, scoring rules and multi-criteria ranking, Journal of   %%
%%       Hydrology, 615, Part B, 2022, doi:10.1016/j.jhydrol.2022.128542.             %%
%%       https://www.sciencedirect.com/science/article/pii/S002216942201112X          %%
%%   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software %%
%%       package: Theory, concepts, and MATLAB implementation, Environmental Modeling %%
%%       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             %%
%%   Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and        %%
%%       J.M. Hyman (2009), Accelerating Markov chain Monte Carlo simulation by       %%
%%       differential evolution with self-adaptive randomized subspace sampling,      %%
%%       International Journal of Nonlinear Sciences and Numerical Simulation, 10(3), %%
%%       271-288                                                                      %%
%%   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and B.A. Robinson (2008), %%
%%       Treatment of input uncertainty in hydrologic modeling: Doing hydrology       %%
%%       backward with Markov chain Monte Carlo simulation, Water Resources Research, %%
%%       44, W00B09, doi:10.1029/2007WR006720                                         %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                                                                                    %%
%% © Written by Jasper A. Vrugt, Feb 2007                                             %%
%% Revised, June 2020                                                                 %%
%% Los Alamos National Laboratory                                                     %%
%% University of California Irvine                                                    %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

if nargin < 9, p_alfa = 0.05; end

alfa1 = 100*p_alfa/2;               % Lower percentile
alfa2 = 100*(1-p_alfa/2);           % Upper percentile
rng(1+round(100*rand),'twister');   % Random seed

npars = size(Up,1);                 % How many parameter vectors?
Ns = sum(Nr);                       % # simulations we expect
fx_mod = nan(Ns,Meas_info.n);       % Initial model output uncertainty

% General code
switch DREAMPar.lik
    case {13,14,16,17,44,45} % NL/GL/LAPL/SL/GL+/UL functions
        par = Lik_info.fpar;        % Initialize fixed values
        for ii = 1:npars
            par(Lik_info.id_vpar) = Up(ii,1:DREAMPar.d);    % Variable parameters
            nuisvar = par(Lik_info.id_nuisvar);             % Nuisance variables
            eval(Lik_info.stringF);                         % Generate replicates
            ii_id = (idU == ii);                            % Indices of replicates
            fx_mod(ii_id,1:Meas_info.n) = fx';              % Store repetitions
        end
        switch npars
            case 1
                fprintf(['DREAM PACKAGE: MAP simulation only was used ' ...
                    'to calculate prediction uncertainty with %s ' ...
                    'function!!!\n'],Lik_info.name_lik_func);
            otherwise
                fprintf(['DREAM PACKAGE: %d simulations of %d unique' ...
                    ' parameter vectors were used to compute ' ...
                    'prediction uncertainty with %s function!!!\n'], ...
                    Ns,npars,Lik_info.name_lik_func);
        end
    otherwise
        if isfield(Meas_info,'Sigma')
            if ~isempty(Meas_info.Sigma)
                std_e = Meas_info.Sigma;
            else
                warning('WARNING: Bayes_PDF: Field ''Sigma'' of structure Meas_info is empty: I use the RMSE_map')
                std_e = RMSE_map;
            end
        else
            std_e = RMSE_map;
        end
        % Now make n copies of Sigma
        if numel(std_e) == 1, std_e = std_e * ones(Meas_info.n,1); end
        % Generate replicates
        for ii = 1:npars
            % fi1 = 0; fi_p = [1 , -fi1];                % Coefficients of AR polynomial
            % eps = normrnd(0,sqrt(1-fi1^2),...          % nxNr(zz) matrix of random variates
            %           Meas_info.n,Nr(zz));
            % e_n = filter(1,fi_p,eps);                  % Stand. raw residuals with var(e_n) = 1
            % Equivalent to
            e_n = normrnd(0,1,Meas_info.n,Nr(ii));          % Draw stand. normal residuals with var(e_n) = 1
            e_r = bsxfun(@times,std_e,e_n);                 % Raw residuals
            fx = bsxfun(@plus,Ufx(ii,1:Meas_info.n)',e_r);  % Replicate simulations
            ii_id = (idU == ii);                            % Indices of replicates
            fx_mod(ii_id,1:Meas_info.n) = fx';              % Store replicates
        end
end

tot_unc = prctile(fx_mod,[alfa1 alfa2])';       % Compute [p_alfa/2 , (1-p_alfa/2)] prediction limits
if Meas_info.n == 1, tot_unc = tot_unc'; end    % Unusual case of only a single measurement

end

% OLD CODE: July 28, 2024
% % function [par_unc,tot_unc,fx_mod] = Bayes_pdf(pars_post,fx_post,...
% %     RMSE_MAP,DREAMPar,Meas_info,Lik_info,p_alfa)
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% % %% Computes confidence and prediction uncertainty of Bayes posterior                  %%
% % %%                                                                                    %%
% % %% SYNOPSIS:                                                                          %%
% % %%  [par_unc,tot_unc,fx_mod] = Bayes_pdf(pars_post,fx_post,RMSE_MAP,DREAMPar,...      %%
% % %%                                 Meas_info,p_alfa)                                  %%
% % %% WHERE                                                                              %%
% % %%  pars_post   [input] Matrix with posterior parameter values                        %%
% % %%  fx_post     [input] Matrix with posterior simulations (parameter uncertainty)     %%
% % %%  RMSE_MAP    [input] Root Mean Square Error of MAP solution                        %%
% % %%  DREAMPar    [input] Structure with algorithmic settings of MCMC algorithm         %%
% % %%  Meas_info   [input] Structure with measurement information (for inference)        %%
% % %%  Lik_info    [input] Structure with information about likelihood function          %%
% % %%  p_alfa      [input] Significance level (e.g: 0.01, 0.05 and 0.1)                  %%
% % %%  par_unc     [outpt] Matrix with low/upper alfa confidence interval parameters     %%
% % %%  tot_unc     [outpt] Matrix with low/upper alfa prediction interval total          %%
% % %%  fx_mod      [outpt] Matrix with simulations of total uncertainty                  %%
% % %%                                                                                    %%
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% % %%                                                                                    %%
% % %% The different componments/algorithms of DREAM Package have been described in       %%
% % %%   Vrugt, J.A., R. de Punder, and P. Grünwald, A sandwich with water:               %%
% % %%       Bayesian/Frequentist uncertainty quantification under model                  %%
% % %%       misspecification, Submitted to Water Resources Research, May 2024            %%
% % %%       https://essopenarchive.org/users/597576/articles/937008-a-sandwich-with-...  %%
% % %%       water-bayesian-frequentist-uncertainty-quantification-under-model-...        %%
% % %%       misspecification                                                             %%
% % %%   Vrugt, J.A. (2024), Distribution-Based Model Evaluation and Diagnostics:         %%
% % %%       Elicitability, Propriety, and Scoring Rules for Hydrograph Functionals,      %%
% % %%       Water Resources Research, 60, e2023WR036710,                                 %%
% % %%       https://doi.org/10.1029/2023WR036710                                         %%
% % %%   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of %%
% % %%       distribution-adaptive likelihood functions: Generalized and universal        %%
% % %%       likelihood functions, scoring rules and multi-criteria ranking, Journal of   %%
% % %%       Hydrology, 615, Part B, 2022, doi:10.1016/j.jhydrol.2022.128542.             %%
% % %%       https://www.sciencedirect.com/science/article/pii/S002216942201112X          %%
% % %%   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software %%
% % %%       package: Theory, concepts, and MATLAB implementation, Environmental Modeling %%
% % %%       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             %%
% % %%   Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and        %%
% % %%       J.M. Hyman (2009), Accelerating Markov chain Monte Carlo simulation by       %%
% % %%       differential evolution with self-adaptive randomized subspace sampling,      %%
% % %%       International Journal of Nonlinear Sciences and Numerical Simulation, 10(3), %%
% % %%       271-288                                                                      %%
% % %%   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and B.A. Robinson (2008), %%
% % %%       Treatment of input uncertainty in hydrologic modeling: Doing hydrology       %%
% % %%       backward with Markov chain Monte Carlo simulation, Water Resources Research, %%
% % %%       44, W00B09, doi:10.1029/2007WR006720                                         %%
% % %%                                                                                    %%
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% % %%                                                                                    %%
% % %% © Written by Jasper A. Vrugt, Feb 2007                                             %%
% % %% Revised, June 2020                                                                 %%
% % %% Los Alamos National Laboratory                                                     %%
% % %% University of California Irvine                                                    %%
% % %%                                                                                    %%
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% %
% % if nargin < 7, p_alfa = 0.05; end
% % alfa1 = 100*p_alfa/2;                       % Lower percentile
% % alfa2 = 100*(1-p_alfa/2);                   % Upper percentile
% % rng(1+round(100*rand),'twister');           % Random seed
% %                                             % legacy: randn('state',sum(100*clock));
% % npars = size(pars_post,1);                  % # parameter vectors
% % par_unc = prctile(fx_post,[alfa1 alfa2])';  % p_alfa/2 & (1-p_alfa/2) confidence limits
% %
% % %% How to generate total uncertainty?
% % %% [1] Sample posterior parameters and add structural uncertainty to this
% % %%     --> does not guarantee that 95% contains 95%
% % %% [2] Take MAP parameter values and add structural uncertainty to this
% % %%     --> parameter uncertainty explains part of total uncertainty
% % method = 2;
% %
% % switch method
% %     case 1 % Sample parameters and add structural uncertainty to this
% %         if DREAMPar.d > 1 % Identify the unique parameter vectors
% %             [R_pars,iiR,idR] = unique(pars_post(1:npars,1:DREAMPar.d), ...
% %                 'rows','stable');
% %         else
% %             [R_pars,iiR,idR] = unique(pars_post(1:npars,1:DREAMPar.d), ...
% %                 'stable');
% %         end
% %         R_fx = fx_post(iiR,1:Meas_info.n);              % Simulations of replicates - and in right order
% %         ntot = size(R_pars,1); Nr = nan(ntot,1);        % How many unique rows in R?
% %         for ii = 1:ntot, Nr(ii) = sum(idR == ii); end   % How many replicates each row?
% %         Ns = sum(Nr);                                   % # simulations we expect
% %         fx_mod = nan(Ns,Meas_info.n);                   % Initial model output uncertainty
% %
% %         % PROBLEM: If Ns is large then this may create memory problems!!
% %         % % nmax = 10000;                               % Total simulations for predictive uncertainty
% %         % % if ( Ns * Meas_info.n ) > 1e8
% %         % %     % We need to shorten output, but how; replicates do not quarantee a
% %         % %     % continous output of fx_mod as indices may not be successive integers
% %         % %     id_max = floor(1e8/Meas_info.n);
% %         % %     ii_row = find(cumsum(N_rep)<=id_max,1,'last');
% %         % %     N_rep = N_rep(1:ii_row); R_fx = R_fx(1:ii_row,1:Meas_info.n);
% %         % %     idR = idR(idR<=ii_row); ntot = ii_row;
% %         % %     Ns = sum(N_rep); fx_mod = nan(Ns,Meas_info.n);
% %         % %     % Now get posterior parameter vectors
% %         % %     R_pars = R_pars(1:ii_row,1:DREAMPar.d+2);
% %         % %     % And recreate pars_post - not right;
% %         % %     pars_post2 = R_pars(idR,:);
% %         % % end
% %
% %         % General code
% %         switch DREAMPar.lik
% %             case {13,14,16,17,44,45} % NL/GL/LAPL/SL/GL+/UL functions
% %                 par = Lik_info.fpar;                                    % Initialize fixed values
% %                 for ii = 1:ntot
% %                     par(Lik_info.id_vpar) = R_pars(ii,1:DREAMPar.d);    % Variable parameters
% %                     nuisvar = par(Lik_info.id_nuisvar);                 % Nuisance variables
% %                     eval(Lik_info.stringF);                             % Generate replicates
% %                     ii_id = (idR == ii);                               % Indices of replicates
% %                     fx_mod(ii_id,1:Meas_info.n) = fx';                  % Store repetitions
% %                 end
% %                 fprintf(['DREAM PACKAGE: %d simulations of %d unique parameter' ...
% %                     ' vectors to calculate total simulation uncertainty ' ...
% %                     'with %s function!!!\n'],Ns,ntot,Lik_info.name_lik_func);
% %             otherwise
% %                 if isfield(Meas_info,'Sigma')
% %                     std_e = Meas_info.Sigma;
% %                 else
% %                     std_e = RMSE_MAP;
% %                 end
% %                 % Now make n copies of Sigma
% %                 if numel(std_e) == 1, std_e = std_e * ones(Meas_info.n,1); end
% %                 % Generate replicates
% %                 for ii = 1:ntot
% %                     % % fi1 = 0;
% %                     % % fi_p = [1 , -fi1];    % Coefficients of AR polynomial
% %                     % % eps = normrnd(0,sqrt(1-fi1^2),Meas_info.n,N_rep(zz));
% %                     % % e_n = filter(1,fi_p,eps);                   % Stand. raw residuals with var(e_n) = 1
% %                     % Equivalent to
% %                     e_n = normrnd(0,1,Meas_info.n,Nr(ii));          % Draw stand. normal residuals with var(e_n) = 1
% %                     e_r = bsxfun(@times,std_e,e_n);                 % Raw residuals
% %                     fx = bsxfun(@plus,R_fx(ii,1:Meas_info.n)',e_r); % Replicate simulations
% %                     ii_id = (idR == ii);                            % Indices of replicates
% %                     fx_mod(ii_id,1:Meas_info.n) = fx';              % Store replicates
% %                 end
% %         end
% %
% %     case 2 % Take MAP parameter values and add structural uncertainty to this
% %
% %         ntot = size(pars_post,1); Nr = ntot;
% %         logP = sum(pars_post(1:ntot,DREAMPar.d+1:DREAMPar.d+2),2);
% %         id_MAP = find(logP == max(logP)); id_MAP = id_MAP(1);
% %         ii = 1; R_fx = fx_post(id_MAP,1:Meas_info.n);               % For stringF
% %         % PROBLEM: If Nr is large then this may create memory problems!!
% %
% %         % General code
% %         switch DREAMPar.lik
% %             case {13,14,16,17,44,45} % NL/GL/LAPL/SL/GL+/UL functions
% %                 par = Lik_info.fpar;                                    % Initialize fixed values
% %                 par(Lik_info.id_vpar) = pars_post(id_MAP,1:DREAMPar.d); % Variable parameters
% %                 nuisvar = par(Lik_info.id_nuisvar);                     % Nuisance variables
% %                 eval(Lik_info.stringF);                                 % Generate replicates
% %                 fx_mod(1:Nr,1:Meas_info.n) = fx';                       % Store repetitions
% %                 fprintf(['DREAM PACKAGE: MAP parameter vector used' ...
% %                     ' to compute prediction uncertainty with' ...
% %                     ' %s function \n'],Lik_info.name_lik_func);
% %             otherwise
% %                 if ~isempty(Meas_info.Sigma)
% %                     std_e = Meas_info.Sigma;
% %                 else
% %                     std_e = RMSE_MAP;
% %                 end
% %                 % Now make n copies of Sigma
% %                 if numel(std_e) == 1, std_e = std_e * ones(Meas_info.n,1); end
% %                 % Generate replicates
% %                 e_n = normrnd(0,1,Meas_info.n,Nr);  % Draw stand. normal residuals with var(e_n) = 1
% %                 e_r = bsxfun(@times,std_e,e_n);     % Raw residuals
% %                 fx = bsxfun(@plus,R_fx',e_r);       % Replicate simulations
% %                 fx_mod(1:Nr,1:Meas_info.n) = fx';   % Store replicates
% %         end
% %
% % end
% %
% % tot_unc = prctile(fx_mod,[alfa1 alfa2])';   % p_alfa/2 & (1-p_alfa/2) prediction limits
% %
% % % In the unusual case that we only have a single measurement
% % if Meas_info.n == 1
% %     tot_unc = tot_unc'; par_unc = par_unc';
% % end
% %
% % end