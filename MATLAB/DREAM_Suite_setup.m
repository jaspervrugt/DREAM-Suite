function [DREAMPar,Par_info,Meas_info,Lik_info,options] = ...
    DREAM_Suite_setup(method,Func_name,DREAMPar,Par_info, ...
    Meas_info,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%             DDDDD     RRRRR     EEEEEEE     AAA     MM    MM            %
%             DDDDDD    RRRRRR    EEEEEEE    AAAAA    MM    MM            %
%             DD  DD    RR   RR   EE        AA   AA   MMM  MMM            %
%             DD   DD   RR  RR    EEEE      AA   AA   MMMMMMMM            %
%             DD   DD   RRRRR     EEEE      AAAAAAA   MMM  MMM            %
%             DD  DD    RR RR     EE        AAAAAAA   MM    MM            %
%             DDDDDD    RR  RR    EEEEEEE   AA   AA   MM    MM            %
%             DDDDD     RR   RR   EEEEEEE   AA   AA   MM    MM            %
%                                                                         %
%              SSSSSSSS  UU    UU   II   TTTTTTTTTT   EEEEEEE             %
%              SSSSSSS   UU    UU   II   TTTTTTTTTT   EEEEEEE             %
%              SS        UU    UU   II       TT       EE                  %
%              SSSS      UU    UU   II       TT       EEEE                %
%                 SSSS   UU    UU   II       TT       EEEE                %
%                   SS   UU    UU   II       TT       EE                  %
%               SSSSSS   UUUUUUUU   II       TT       EEEEEEE             %
%              SSSSSSS   UUUUUUUU   II       TT       EEEEEEE             %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% This function initializes the main variables used in DREAM-Suite        %
%                                                                         %
% SYNOPSIS: [DREAMPar,Par_info,Meas_info,Lik_info,options] = ...          %
%               DREAM_Suite_setup(method,Func_name,DREAMPar, ...          %
%               Par_info,Meas_info,options)                               %
%                                                                         %
% © Written by Jasper A. Vrugt, Feb 2007                                  %
% Los Alamos National Laboratory 			                              %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Random seed (legacy: randn('state',sum(100*clock)); )
rng(1+round(100*rand),'twister');

switch method
    case {'dream','dream_d'}
        % Name variable
        name = {'nCR','delta','steps','lambda','zeta','p_unit_gamma',...
            'adapt_pCR','thinning','beta0','GLUE','outlier', ...
            'pparallel','psnooker','pkalman','mt'};
        % Default values algorithmic variables DREAM - if not specified
        value = {'3','3',num2str(max(max(floor(DREAMPar.T/50),1),50)),...
            '0.05','1e-12','0.2','''yes''','1','1','10','''iqr''',...
            '1','0','0','1'};
    case {'dream_zs','dream_dzs'}
        % Name variable
        name = {'nCR','delta','steps','lambda','zeta','p_unit_gamma',...
            'adapt_pCR','thinning','beta0','GLUE','N','k','psnooker',...
            'm0','pparallel','psnooker','pkalman','mt'};
        % Default values algorithmic variables DREAM - if not specified
        value = {'3','3',num2str(max(max(floor(DREAMPar.T/50),1),50)),...
            '0.05','1e-12','0.2','''yes''','1','1','10','3','10','0.1',...
            ['max(10*DREAMPar.d,max(3*DREAMPar.N,' ...
            '2*DREAMPar.delta*DREAMPar.N))'],'0.9','0.1','0','1'};
    case {'dream_kzs'}
        % Name variable
        name = {'nCR','delta','steps','lambda','zeta','p_unit_gamma',...
            'adapt_pCR','thinning','beta0','GLUE','N','k','psnooker',...
            'm0','pparallel','psnooker','pkalman','M','a_1','a_2','mt'};
        % Default values algorithmic variables DREAM - if not specified
        value = {'3','3',num2str(max(max(floor(DREAMPar.T/50),1),50)),...
            '0.05','1e-12','0.2','''yes''','1','1','10','3','10','0.1',...
            ['max(10*DREAMPar.d,max(3*DREAMPar.N,' ...
            '2*DREAMPar.delta*DREAMPar.N))'],'0.5','0.1','0.4',...
            ['min(20,max(10*DREAMPar.d,max(3*DREAMPar.N,' ...
            '2*DREAMPar.delta*DREAMPar.N)))'],'0.1','0.2','1'};
    case {'mtdream_zs'}
        % Check whether snooker probability is defined by user
        if isfield(DREAMPar,'psnooker')
            DREAMPar.pparallel = 1 - DREAMPar.psnooker;
        end
        % Otherwise, parallel direction probability
        if isfield(DREAMPar,'pparallel')
            DREAMPar.psnooker = 1 - DREAMPar.pparallel;
        end
        % Name variable
        name = {'N','delta','k','mt','psnooker','pparallel','nCR', ...
            'steps','lambda','zeta','p_unit_gamma','adapt_pCR', ...
            'thinning','beta0','GLUE'};
        % Set default values algorithmic variables DREAM - if not specified
        value = {'3','1','10','5','0.1','0.9','3', ...
            num2str(max(max(floor(DREAMPar.T/50),1),50)),'0.05', ...
            '1e-12','0.2','''yes''','1','1','10'};
end

% Now check
for j = 1 : numel(name)
    if ~isfield(DREAMPar,name(j))
        % Set variable of DREAMPar
        evalstr = strcat('DREAMPar.',char(name(j)),'=',value(j),';');
        eval(char(evalstr));
    end
end

% Set m = m0 if relevant method
if isfield(DREAMPar,'m0')
    % Set DREAMPar.m
    DREAMPar.m = DREAMPar.m0;
elseif strcmpi(method,'mtdream_zs')
    % Assign m0 if not specified
    DREAMPar.m0 = max ( 10 * DREAMPar.d , max ( 3 * DREAMPar.N , ...
        2 * DREAMPar.delta * DREAMPar.N ) );
    % Set DREAMPar.m
    DREAMPar.m = DREAMPar.m0;
end

% Set default value to 'No' if not specified
name = {'parallel','IO','modout','save','restart','DB','epsilon',...
    'diagnostics','print','burnin'};
% Set default values algorithmic variables DREAM - if not specified
value = {'''no''','''no''','''no''','''no''','''no''','''no''','0.025',...
    '''yes''','''yes''','50'};

% Set to "No" those that are not specified
for j = 1 : numel(name)
    if ~isfield(options,name(j))
        % Set variable of DREAMPar to "No"
        evalstr = strcat('options.',char(name(j)),'=',value(j),';');
        eval(char(evalstr));
    end
end

% Approximate Bayesian computation
if ( DREAMPar.lik > 20 ) && ( DREAMPar.lik < 24 )
    % We use ABC approach
    options.ABC = 'yes';
    % Now check rho
    if ~isfield(options,'rho')
        % Specify the distance function as anonymous function handle
        options.rho = @(X,Y) abs ( X - Y );
    end
else
    options.ABC = 'no';
end

switch method
    case {'dream','dream_d'}
        % Matrix DREAMPar.R: Store for each chain (as row) the index of
        % all other chains available for DE
        for i = 1:DREAMPar.N, z = 1:DREAMPar.N; z(i) = [];
            DREAMPar.R(i,1:DREAMPar.N-1) = z; end
        % Define psnooker and pkalman as zero
        DREAMPar.psnooker = 0; DREAMPar.pkalman = 0;
        % Calculate selection probability of parallel direction jump
        DREAMPar.pparallel = 1 - DREAMPar.psnooker - DREAMPar.pkalman;
    case {'dream_zs','dream_dzs','dream_kzs'}
        % Determine sample indices for proposal generation
        if DREAMPar.delta == 1
            % Define DREAMPar select
            DREAMPar.select =  3 * DREAMPar.N;
            % Fixed indices randomly chosen samples of Z
            DREAMPar.R = reshape ( 1 : DREAMPar.select , 3 , DREAMPar.N )';
        elseif DREAMPar.delta > 1
            % Define DREAMPar select
            DREAMPar.select =  2 * DREAMPar.delta * DREAMPar.N;
            % Fixed indices randomly chosen samples of Z
            DREAMPar.R = reshape ( 1 : DREAMPar.select , ...
                2 * DREAMPar.delta , DREAMPar.N )';
        end
        % If likelihood function larger than 21 (22 or 23) or DB then no
        % snooker update allowed!
        if ( ismember(DREAMPar.lik,[22 23]) || ...
                strcmp(options.DB,'yes') ), DREAMPar.psnooker = 0; end
        % Set Kalman jump probability to zero in current case group
        % for DREAM_ZS or DREAM_DZS
        if any(strcmp(method,{'dream_zs','dream_dzs'}))
            DREAMPar.pkalman = 0;
        end
        % Calculate selection probability of parallel direction jump
        DREAMPar.pparallel = 1 - DREAMPar.psnooker - DREAMPar.pkalman;
    case {'mtdream_zs'}
        % Determine sample indices for proposal generation
        if DREAMPar.delta == 1
            % Define DREAMPar select
            DREAMPar.select =  3 * DREAMPar.mt;
            % Fixed indices randomly chosen samples of Z
            DREAMPar.R = reshape ( 1 : DREAMPar.select , 3 , DREAMPar.mt );
        elseif DREAMPar.delta > 1
            % Define DREAMPar select
            DREAMPar.select =  2 * DREAMPar.delta * DREAMPar.mt;
            % Fixed indices randomly chosen samples of Z
            DREAMPar.R = reshape ( 1 : DREAMPar.select , ...
                2 * DREAMPar.delta , DREAMPar.mt );
        end
        % In DREAM_Suite we use DREAMPar.R = DREAMPar.R' !!!
end

% Check whether we work in normalized parameter space or not
if Par_info.norm == 1
    Par_info.minun = Par_info.min; Par_info.maxun = Par_info.max;
    Par_info.min = zeros(1,DREAMPar.d); Par_info.max = ones(1,DREAMPar.d);
else
    % do nothing
    % Par_info.minun = Par_info.min; Par_info.maxun = Par_info.max;
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    % Now lets define step size
    Par_info.step_size = ( Par_info.max - Par_info.min ) ./ Par_info.steps;
end

% # training data observations
if isfield(Meas_info,'Y')
    Meas_info.n = numel(Meas_info.Y);
else
    Meas_info.n = 0;
end

% Do we use summary metrics?
if isfield(Meas_info,'S')
    Meas_info.n_S = numel(Meas_info.S);
    % Copy epsilon n_S times, unless epsilon is already specified as vector
    options.epsilon = ones(Meas_info.n_S,1) .* options.epsilon(:);
else
    Meas_info.n_S = 0;
end

if isfield(Meas_info,'Sigma')
    if numel(Meas_info.Sigma) == 1
        % Copty Sigma n different times so we yield a nx1 vector
        Meas_info.Sigma = ones(Meas_info.n,1) * Meas_info.Sigma;
    end
else
    Meas_info.Sigma = [];
end

if ( DREAMPar.lik == 52 )
    if isfield(Meas_info,'C')
        % Compute determinant and inverse of covariance matrix, Meas_info.C
        Meas_info.detC = det(Meas_info.C); Meas_info.invC = inv(Meas_info.C);
    else
        % do nothing
    end
end

% Define prior handle (random samples & evaluation of pdf)
if isfield(Par_info,'prior')
    if iscell(Par_info.prior)
        Par_info.u = 'yes';                 % Univariate case
        for ii = 1 : DREAMPar.d             % Anonymous handle each paramtr
            Par_info.prior_rnd{ii} = ...    % Handle draw initial state
                eval(char(strcat('@(x)',{' '}, ...
                char(strrep(Par_info.prior(ii), ...
                'pdf','rnd')))));
            Par_info.prior{ii} = ...        % Handle prior pdf each paramtr
                eval(char(strcat('@(x)',{' '}, ...
                char(strrep(Par_info.prior(ii), ...
                'pdf(','pdf(x,')))));
        end
    else
        Par_info.u = 'no';                  % Multivariate case
        pr_name = char(Par_info.prior);     % Turn handle to string
        [pr_var,idcp] = ...                 % Extract variable names prior
            extract_names(pr_name);         % & index closing parenthesis
        pr_name = pr_name(idcp(1)+1:end);   % Prior handle without @()
        n_var = numel(pr_var);              % # variabls without x (=n_p-1)
        % n_p = nargin(Par_info.prior)
        for z = 1:n_var                     % Add Par_info. to varble names
            pr_name = strrep(pr_name,...
                pr_var{z},strcat('Par_info.',pr_var{z}));
        end
        Par_info.prior_rnd = ...            % Handle draw initial states
            eval(char(strcat('@(x)',{' '}, ...
            char(strrep(pr_name, ...
            'pdf(x,','rnd(')))));
        Par_info.prior = ...                % Handle pdf each paramtr vectr
            eval(char(strcat('@(x)',{' '}, ...
            pr_name)));
    end
    % Now determine whether user returned the pdf or log(pdf)
    M = 100; Par_info.pr = check_prior(Par_info,DREAMPar,M);
    % Par_info.pr = 'pdf' or Par_info.pr = 'logpdf'
    fprintf(char(strcat(['DREAM_Suite: Code analyzed that prior' ...
        ' handle returns a'],{' '},Par_info.pr,{' '},'\n')));
end

% Remove unnecessary fields
switch method
    case {'dream','dream_d'}
        remove_fields = {'k','m0','m'};
    case {'dream_zs','dream_dzs','dream_kzs'}
        remove_fields = {'outlier'};
    case {'mtdream_zs'}
        remove_fields = '';
end
for j = 1 : numel(remove_fields)
    if isfield(DREAMPar,remove_fields(j))
        eval(strcat('DREAMPar = rmfield(DREAMPar,remove_fields(j))'));
    end
end

% Order and print to screen DREAMPar fields
switch method
    case {'dream','dream_d'}
        DREAMPar = orderfields(DREAMPar,{'d','N','T','lik','delta',...
            'nCR','lambda','zeta','p_unit_gamma','thinning','beta0',...
            'GLUE','adapt_pCR','outlier','steps','R','pparallel',...
            'psnooker','pkalman','mt'});
        not_print = 5;
    case {'dream_zs','dream_dzs'}
        DREAMPar = orderfields(DREAMPar,{'d','N','T','lik','delta',...
            'k','pparallel','psnooker','m0','nCR','lambda','zeta',...
            'p_unit_gamma','thinning','beta0','GLUE','adapt_pCR','m',...
            'steps','select','R','pkalman','mt'});
        not_print = 5;
    case {'dream_kzs'}
        DREAMPar = orderfields(DREAMPar,{'d','N','T','lik','delta',...
            'k','pparallel','psnooker','pkalman','M','a_1','a_2','m0',...
            'nCR','lambda','zeta','p_unit_gamma','thinning','beta0',...
            'GLUE','adapt_pCR','m','steps','select','R','mt'});
        not_print = 3;
    case {'mtdream_zs'}
        % Order and print to screen DREAMPar fields
        DREAMPar = orderfields(DREAMPar,{'d','N','T','lik','delta',...
            'k','mt','pparallel','psnooker','m0','nCR','lambda','zeta',...
            'p_unit_gamma','thinning','beta0','GLUE','adapt_pCR','m',...
            'steps','select','R'});
        not_print = 5;
end

% Now print to screen
fprintf('--- Summary of algorithmic settings ---\n');
fprintf('  DREAMPar\n');
F = fieldnames(DREAMPar);
for i = 1 : numel(F) - not_print
    pr = num2str(eval(char(strcat('DREAMPar.',F(i,:)))));
    fprintf(char(strcat('    .',F(i,:),{' '},'=',{' '},pr,'\n')));
end
fprintf('----- End of algorithmic settings -----\n');

% Now add this [at bottom so that field does not print to output]
if strcmp(method,'dream_kzs')
    % Now save original value of DREAMPar.pkalman
    DREAMPar.oldpkalman = DREAMPar.pkalman;
    if DREAMPar.a_1 > 0
        DREAMPar.pkalman = 0;
    else
        % DREAMPar.pkalman is what has been specified by user (or default)
    end
end

% Create new structure Lik_info with information about likelihood function
[Lik_info,DREAMPar,Par_info] = setup_lik(Func_name,DREAMPar,Par_info,...
    Meas_info); disp(Lik_info)

end