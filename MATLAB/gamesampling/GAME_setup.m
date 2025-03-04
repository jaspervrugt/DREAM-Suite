function [DREAMPar,metric,K,N,M,steps,Par_info,Meas_info,Lik_info, ...
    options,slash_dir] = GAME_setup(DREAMPar,Func_name, ...
    user_options,Par_info,Meas_info,options);                       %#ok
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
% GGGG   AA   M   M  EEEE  SSSS   AA   M   M  PPPP  L     I  N    N  GGGG %
% G  G  A  A  M   M  E     S     A  A  M   M  P  P  L     I  NN   N  G  G %
% G  G  A  A  MM MM  EEE   S     A  A  MM MM  P  P  L     I  N N  N  G  G %
% GGGG  AAAA  M M M  EEE - SSSS  A  A  M M M  PPPP  L     I  N  N N  GGGG %
%    G  A  A  M   M  E        S  AAAA  M   M  P     L     I  N   NN     G %
% GGGG  A  A  M   M  EEEE  SSSS  A  A  M   M  P     LLLL  I  N    N  GGGG %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Initializes the main variables used in GAME sampling                    %
%                                                                         %
% SYNOPSIS: [DREAMPar,metric,K,N,M,steps,Par_info,Meas_info,Lik_info, ... %
%               GAMEoptions,slash_dir] = GAME_setup(DREAMPar, ...         %
%               Func_name,user_options,Par_info,Meas_info,options)        %
%                                                                         %
% © Written by Jasper A. Vrugt, Aug. 2015                                 %
% University of California Irvine                                         %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

rng(1+round(100*rand),'twister');   % Random seed

% Define operating system
options.PC = 'no'; options.MAC = 'no'; options.UNIX = 'no';

% Check whether running on PC or MAC/UNIX machine
if ispc, slash_dir = '\'; else, slash_dir = '/'; end

% Now check the operating system
if ispc, options.PC = 'yes'; end
if ismac, options.MAC = 'yes'; end
if isunix, options.UNIX = 'yes'; end

% Set default value to 'No' if not specified
name = {'parallel','IO','modout','save','restart','DB','epsilon'};
% Set default values algorithmic variables DREAM - if not specified
value = {'''no''','''no''','''no''','''no''','''no''','''no''','0.025'};

% Set to "No" those that are not specified
for j = 1 : numel(name)
    if ~isfield(options,name(j))
        % Set variable of DREAMPar to "No"
        evalstr = strcat('options.',char(name(j)),'=',value(j),';'); 
        eval(char(evalstr));
    end
end
% No look at default settings GAME_options
def_GAMEoptions = struct('metric','bic','K',5,'M',1,'N',1e4,'steps',10);
% Extract the various fields from structure def_options
names = fieldnames(def_GAMEoptions);
% Which fields has the use specified?
user_names = fieldnames(user_options);
% Now create structure options as mix of def_GAME_options and user_options
GAMEoptions = def_GAMEoptions;                                      %#ok
% Now check each field name
for i = 1 : numel(names)
    for j = 1 : numel(user_names)
        % Always use lower case for content of each field
        if strcmpi(names(i),user_names(j))
            evalstr = strcat('GAMEoptions.',char(names(i)), ...
                ' = lower(user_options.',char(user_names(j)),');'); 
            eval(char(evalstr));
        end
    end
end
% Now unpack the fields from structure options
for i = 1 : numel(names)
    evalstr = strcat(names(i),' = GAMEoptions.',names(i),';'); 
    eval(char(evalstr));
end

% Initialize DREAMPar.CPU to be one (if 'ris' is used -> for GAME_end)
DREAMPar.CPU = 1;

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

% Check whether we work in normalized parameter space or not
if Par_info.norm == 1
    Par_info.minun = Par_info.min; Par_info.maxun = Par_info.max;
    Par_info.min = zeros(1,DREAMPar.d); Par_info.max = ones(1,DREAMPar.d);
else
    % do nothing
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
    if isscalar(Meas_info.Sigma)
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
    Par_info.pr = check_prior(Par_info,DREAMPar,100);
    % Difference with eDREAM_package_setup as M is variable in GAMEoptions
    % Par_info.pr = 'pdf' or Par_info.pr = 'logpdf'
    fprintf(char(strcat(['GAME_sampling: Code analyzed that prior' ...
        ' handle returns a'],{' '},Par_info.pr,{' '},'\n')));
end

% Create new structure Lik_info with information about likelihood function
[Lik_info,DREAMPar,Par_info] = setup_lik(Func_name,DREAMPar,Par_info,...
    Meas_info); disp(Lik_info)

% Now order in ASCII dictionary order
DREAMPar = orderfields(DREAMPar); Par_info = orderfields(Par_info); 
Meas_info = orderfields(Meas_info); options = orderfields(options);

% Upper case
fprintf('\n');
fprintf(['  ---------------------------------------------------------' ...
    '--------------------------------------------------  \n']);
fprintf('      GGGGGGG     A     M     M  EEEEEE   SSSSSSS    A     M     M  PPPPPP   L       I  N     N  GGGGGGG      \n');
fprintf('      G     G    A A    M     M  E        S         A A    M     M  P     P  L       I  NN    N  G     G      \n');
fprintf('      G     G   A   A   M     M  E        S        A   A   M     M  P     P  L       I  N N   N  G     G      \n');
fprintf('      G     G  A     A  MM   MM  E        S       A     A  MM   MM  P     P  L       I  N N   N  G     G      \n');
fprintf('      GGGGGGD  A     A  M M M M  EEEEE -- SSSSSS  A     A  M M M M  PPPPPP   L       I  N  N  N  GGGGGGG      \n');
fprintf('            G  AAAAAAA  M  M  M  E             S  AAAAAAA  M  M  M  P        L       I  N   N N        G      \n');
fprintf('            G  A     A  M     M  E             S  A     A  M     M  P        L       I  N   N N        G      \n');
fprintf('            G  A     A  M     M  E             S  A     A  M     M  P        L       I  N    NN       GG      \n');
fprintf('      GGGGGGG  A     A  M     M  EEEEEE   SSSSSS  A     A  M     M  p        LLLLLL  I  N     N  GGGGGGG      \n');
fprintf(['  ---------------------------------------------------------' ...
    '--------------------------------------------------  \n']);
fprintf('  © Jasper A. Vrugt, University of California Irvine \n');
fprintf('\n');

end
