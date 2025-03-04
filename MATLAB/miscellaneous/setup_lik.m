function [Lik_info,DREAMPar,Par_info] = setup_lik(Func_name,DREAMPar,...
    Par_info,Meas_info)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Function compiles the necessary information for the distribution adaptive          %%
%% likelihood functions (including NL and LAPL)                                       %%
%%                                                                                    %%
%% SYNOPSIS: [Lik_info,DREAMPar,Par_info] = setup_lik(Func_name,DREAMPar,...          %%
%%              Par_info,Meas_info)                                                   %%
%%  where                                                                             %%
%%  Func_name   [input] Function (string) which returns (log)lik or simulated values  %%
%%  DREAMPar    [input] Structure with algorithmic variables                          %%
%%  Par_info    [input] Parameter structure: Ranges, initial/prior and bnd treatmnt   %%
%%  Meas_info   [input] Structure with measurement information (for inference)        %%
%%  Lik_info    [outpt] Structure with information about likelihood function          %%
%%  DREAMPar    [outpt] Structure with algorithmic variables                          %%
%%  Par_info    [outpt] Parameter structure: Ranges, initial/prior and bnd treatmnt   %%
%%                                                                                    %%
%% © Written by Jasper A. Vrugt, June 2016                                            %%
%% University of California Irvine                                                    %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Name of likelihood fu-nction
switch DREAMPar.lik
    case 1 % User supplies likelihood
        name_lik = 'Own likelihood'; name_lik_func = Func_name;
    case 2 % User supplies log-likelihood
        name_lik = 'Own log-likelihood'; name_lik_func = Func_name;
    case 11 % Box and Tiao (1953) variant of log-likelihood
        name_lik = 'Gaussian likelihood without sigma';
        name_lik_func = 'Normal';
    case 12 % Normal likelihood with homos/heteroscedastic meas. error
        name_lik = 'Gaussian likelihood'; name_lik_func = 'Normal';
    case 13 % Normal likelihood with AR(2) & (non)constant residual var.
        name_lik = 'Normal likelihood with AR(2)-process';
        name_lik_func = 'NL';
        lik_names = {'s_{0}','s_{1}','\phi_{1}','\phi_{2}'};
        ar_id = [3 4];
    case 14 % Generalized likelihood with AR(2) & (non)constant res. var.
        name_lik = 'Generalized Likelihood';
        name_lik_func = 'GL';
        lik_names = {'s_{0}','s_{1}','\beta','\xi','\mu_{1}',...
            '\phi_{1}','\phi_{2}','\phi_{3}','\phi_{4}','K',...
            '\lambda'};
        ar_id = [6 7 8 9];
    case 15 % Spectral likelihood of Whittle, 1952, 1953
        name_lik = 'Whittle (spectral) likelihood';
        name_lik_func = 'Whittle';
    case 16 % Laplace likelihood with AR(1) & (non)constant residual var.
        name_lik = 'Laplacian likelihood with AR(1) process';
        name_lik_func = 'LAPL';
        lik_names = {'s_{0}','s_{1}','\phi_{1}'};
        ar_id = 3;
    case 17 % Student t likelihood with AR(2) & (non)constant residual var.
        name_lik = 'Skewed Student T likelihood';
        name_lik_func = 'SL';
        lik_names = {'s_{0}','s_{1}','\nu','\xi','\phi_{1}',...
            '\phi_{2}'};
        ar_id = [5 6];
    case 21 % Approximate Bayesian Computation: Turner formulation
        name_lik = 'Approximate Bayesian Computation';
        name_lik_func = 'ABC';
    case 22 % Approximate Bayesian Computation: Sadegh and Vrugt, 2015
        name_lik = ['Approximate Bayesian Computation ' ...
            '- alternative kernel']; name_lik_func = 'ABC-kernel';
    case 23 % Approximate Bayesian Computation: Vrugt and Beven, 2018
        name_lik = ['Approximate Bayesian Computation - ' ...
            'limits of acceptability']; name_lik_func = 'ABC-LOA';
    case 31 % Informal likelihood of Beven and Freer, 2001
        name_lik = 'Informal likelihood: Beven and Freer, 2001';
        name_lik_func = 'Informal_BFa';
    case 32 % Informal likelihood of Beven and Freer, 2001
        name_lik = ['Informal likelihood: Option b in ' ...
            'Beven and Freer, 2001']; name_lik_func = 'Informal_BFb';
    case 33 % Informal likelihood of Beven and Freer, 2001
        name_lik = ['Informal likelihood: Option c in ' ...
            'Beven and Freer, 2001']; name_lik_func = 'Informal_BFc';
    case 34 % Informal likelihood of Beven and Binley, 1992
        name_lik = ['Informal likelihood: Last option, ' ...
            'Page 284 of Beven and Binley, 1992'];
        name_lik_func = 'Informal_BFd';
    case 44 % Generalized likelihood ++, AR(2) & (non)constant residual var.
        name_lik = 'Generalized Likelihood ++';
        name_lik_func = 'GL_plus';
        lik_names = {'s_{0}','s_{1}','\beta','\xi','\phi_{1}',...
            '\phi_{2}'};
        ar_id = [5 6];
    case 45 % Universal likelihood ++, AR(2) & (non)constant residual var.
        name_lik = 'Universal likelihood';
        name_lik_func = 'UL';
        lik_names = {'s_{0}','s_{1}','\lambda','p','q',...
            '\phi_{1}','\phi_{2}'};
        ar_id = [6 7];
    case 52 % Generalized least squares form of likelihood function
        name_lik = ['Matrix implementation of Gaussian likelihood ' ...
            'using GLS form']; name_lik_func = 'Normal_GLS';
        % POWER LIKELIHOODS
    case 61 % Laplace power likelihood with unit integral
        name_lik = 'Power likelihood: Laplace distribution';
        name_lik_func = 'Power_laplace';
    case 62 % Normal power likelihood with unit integral
        name_lik = 'Power likelihood: Normal distribution';
        name_lik_func = 'Power_normal';
    case 63 % Laplace power likelihood with fixed lambda
        name_lik = 'Power likelihood: Laplace distribution - fixed lambda';
        name_lik_func = 'Power_laplace';
    case 64 % Normal power likelihood with fixed lambda
        name_lik = 'Power likelihood: Normal distribution - fixed lambda';
        name_lik_func = 'Power_normal';
end

% Now determine properties of likelihood function
switch DREAMPar.lik
    case {13,14,16,17,44,45} % NL/GL/LAPL/SL/GL+/UL functions
        global LV                               % request global variables
        nest = numel(LV.id_vpar);               % number of estimeable parameters
        ntot = numel(LV.fpar);                  % total number of parameters
        nf = ntot - nest;                       % number of fixed nuisance variables
        nmod = nest - numel(lik_names) + nf;    % number of model parameters
        % check whether s0 and s1 are selected
        if numel(LV.id_vpar) >= nmod+1
            switch LV.id_vpar(nmod+1)
                case nmod+1 % so is estimated
                    s0 = 1;
                    if numel(LV.id_vpar) >= nmod+2
                        if LV.id_vpar(nmod+2) == nmod+2
                            s1 = 1;
                        else
                            s1 = 0;
                        end
                    else
                        s1 = 0;
                    end
                case nmod+2 % s0 not estimated but s1 is
                    s0 = 0; s1 = 1;
                otherwise % s0 and s1 not estimated
                    s0 = 0; s1 = 0;
            end
        else
            s0 = 0; s1 = 0;
        end
        % Now check whether user specified Sigma or not
        if ~isempty(Meas_info.Sigma)
            method = 0;                     % no s0/s1 needed: Sigma is defined by user!
        else
            % Now check treatment (method) of measurement error variance
            switch Meas_info.sigma2
                case 'constant'
                    if (s0 == 0)                % CONSTANT: s0 not selected
                        if (s1 == 0)
                            method = 1;         % s0: fixed
                        elseif (s1 == 1)        % CONSTANT: s0 not selected but s1
                            method = 2;         % s0: phantom variable [s1 inconsequential]
                        end
                    elseif (s0 == 1)            % CONSTANT: s0 selected
                        if (s1 == 0)
                            method = 3;         % s0: estimated
                        elseif (s1 == 1)        % CONSTANT: s0 and s1 selected
                            method = 4;         % s0: estimated [s1 inconsequential]
                        end
                    end
                case 'nonconstant'
                    if (s0 == 0)                % NONCSTNT: s0 not selected
                        if (s1 == 0)
                            method = 5;         % s0: fixed, s1: phantom variable
                        elseif (s1 == 1)        % NONCSTNT: s0 not selected but s1
                            method = 6;         % s0: fixed, s1: estimated
                        end
                    end
                    if (s0 == 1)                % NONCSTNT: s0 selected
                        if (s1 == 0)
                            method = 7;         % s0: estimated, s1: phantom variable
                        elseif (s1 == 1)        % NONCONSTANT: s0 and s1 selected
                            method = 8;         % s0: estimated, s1: estimated
                        end
                    end
            end
        end
        % Thus with a constant variance, method = 2 and method = 4 should
        % lead to removal of s1 from the list of parameters
        id_rem = [];
        switch method
            case 0
                if s0 == 0
                    if s1 == 1                        
                        id_rem = nmod + 1;
                    end
                elseif s0 == 1
                    id_rem = nmod + 1;
                    if s1 == 1                        
                        id_rem(2) = nmod + 2;
                    end 
                end
            case 2
                id_rem = nmod + 1;
            case 4
                id_rem = nmod + 2;
        end
        % Remove id_rem entry from variable parameters & DREAMPar
        LV.id_vpar(id_rem) = []; DREAMPar.d = DREAMPar.d - numel(id_rem);
        % Remove from Par_info
        Par_info.min(id_rem) = []; Par_info.max(id_rem) = [];
        % Remove from Par_info.names
        if isfield(Par_info,'names')
            Par_info.names(id_rem) = [];
        end
        % Remove from Par_info.steps
        if isfield(Par_info,'steps')
            Par_info.steps(id_rem) = [];
            Par_info.step_size(id_rem) = [];
        end
        % Now continue with regular program
        fpar = LV.fpar; id_vpar = LV.id_vpar;   % extract variables
        par = fpar;                             % duplicate copy
        par(id_vpar) = nan;                     % variable parameters
        nest = numel(id_vpar);                  % number of estimeable parameters
        nf = ntot - nest;                       % number of fixed nuisance variables
        id_nuisvar = nmod+1:ntot;               % index of nuisvar
        nuisvar = par(id_nuisvar);              % isolate nuisance variables
        id_ar = nmod + ar_id;                   % index of AR coefficients
        id_par = zeros(1,ntot);                 % vector of variable pars
        id_par(id_vpar) = 1;                    % set variable pars to 1
        t1 = 1 + sum(id_par(id_ar));            % 1st measurement AR model
        str_nuis = lik_names(isnan(nuisvar));   % get the likelihood parameters
        % Add $$ signs in front and end of str_lik - for latex print
        str_nuis = insertBefore(str_nuis,1,'$');
        for i = 1:numel(str_nuis)
            str_nuis(i) = insertAfter(str_nuis(i),...
                numel(char(str_nuis(i))),'$');
        end
        % Below = new implementation which returns nx1 vector of log-likelihoods
        lik_strL = strcat('[loglik,std_e,eps_n,f_eps_n] = ',{' '},...
            name_lik_func,'(''est'',nuisvar,',num2str(method),...
            ',FXp(1:n,ii),Meas_info.Y,Meas_info.Sigma);');
        lik_strF = char(strcat('[~,~,~,~,fx] = ',{' '},name_lik_func,...
            ['(''sim'',nuisvar,',num2str(method),',Ufx(ii,1:Meas_info.n)'',' ...
            'Meas_info.Y,Meas_info.Sigma,Nr(ii));']));
        lik_strL = char(lik_strL); lik_strE = []; %lik_strE = char(lik_strE);
        lik_strF = char(lik_strF); str_sigma = [];
        if isfield(LV,'filename'), filename = 1; end
    otherwise
        str_sigma = []; str_nuis = []; nmod = 0; t1 = 1;
        lik_strL = []; lik_strE = []; lik_strF = []; fpar = [];
        id_vpar = []; id_ar = []; id_nuisvar = []; filename = 0; method = 0;
end

% Now combine str_lik with str_sigma to get all nuisance variables
str_nuis = [str_sigma,str_nuis];
% Now assign fields to Lik_info
Lik_info = struct('name_lik',name_lik,'name_lik_func',name_lik_func,...
    'fpar',fpar,'id_vpar',id_vpar,'index',nmod,'stringL',...
    lik_strL,'stringE',lik_strE,'stringF',lik_strF,'str_nuis',{str_nuis},...
    't_start',t1,'id_ar',id_ar,'id_nuisvar',id_nuisvar,'method',method,...
    'filename',filename);
% NOTE: str_nuis between braces otherwise multiple copies of Lik_info!!

end