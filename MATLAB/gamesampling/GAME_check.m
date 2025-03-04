function [method,DREAMPar,Par_info,options,GAMEoptions,R,d2] = ...
    GAME_check(method,Xp,Func_name,DREAMPar,GAMEoptions,Par_info, ...
    Meas_info,options)
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
% This function verifies the arguments defined by user for GAME sampler   %
%                                                                         %
% SYNOPSIS: [method,DREAMPar,Par_info,options,GAMEoptions,R,d2] = ...     %
%               GAME_check(method,Xp,Func_name,DREAMPar,GAMEoptions, ...  %
%               Par_info,Meas_info,options)                               %
%                                                                         %
% © Written by Jasper A. Vrugt, Jan 2015                                  %
% University of California Irvine                                         %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Determine the number of samples of Xp
[R,d2] = size(Xp);
if (DREAMPar.d+2 ~= d2)
    error(['GAME_sampling ERROR: Matrix Xp should be of ' ...
        'size R x DREAMPar.d+2']);
end

% <><><><><><><><><><><><><><><><> method <><><><><><><><><><><><><><><><><
% Check that method is a string
if ~ischar(method)
    error('GAME_sampling ERROR: method should be a string');
else
    method = lower(method);
end
% Now check if method exists!
if ~sum(strncmp(method,{'ris','is','gb','ob'},inf))
    warning(['GAME_sampling WARNING: Unknown marginal likelihood ' ...
        'estimation method -> Define method = ''ris'' ' ...
        'or ''is'' or ''gb'' or ''ob''!!\n']);
    % Continue without ris
    warning(['GAME_sampling WARNING: Code resorts to default ' ...
        'estimator: reciprocal importance sampling: method = ''ris'' \n']);
    % Now print warning to screen and to file
    method = 'ris';
end

% <><><><><><><><><><><> GAMEoptions & other inputs <><><><><><><><><><><><
if isfield(GAMEoptions,'metric')
    % Check that method is a string
    if ~ischar(GAMEoptions.metric)
        error(['GAME_sampling ERROR: Field metric of GAMEoptions ' ...
            'should be a string with ''bic'' or ''var'' respectively ']);
    else
        GAMEoptions.metric = lower(GAMEoptions.metric);
    end
end

if isfield(GAMEoptions,'metric')
    if ~sum(strncmp(GAMEoptions.metric,{'bic','var'},inf))
        warning(['GAME_sampling WARNING: Unknown metric in field ' ...
            'metric of GAME_options -> Define GAMEoptions.metric = ' ...
            '''bic'' or ''var''!!\n']);
        % Continue without ris
        warning(['GAME_sampling WARNING: Code resorts to default ' ...
            'selection metric: ''bic''\n']);
        % Now print warning to screen and to file
        GAMEoptions.metric = 'bic';
    end
end

if isfield(GAMEoptions,'J')
    if ischar(GAMEoptions.J)
        error(['GAME_sampling ERROR: The field J of GAMEoptions ' ...
            'should be an integer with maximum number of components ' ...
            'of Gaussian mixture distribution']);
    end
end
if isfield(GAMEoptions,'N')
    if ischar(GAMEoptions.N)
        error(['GAME_sampling ERROR: The field N of GAMEoptions ' ...
            'should be an integer with number of importance samples ' ...
            'that is used to compute the marginal likelihood']);
    end
end
if isfield(GAMEoptions,'M')
    if ischar(GAMEoptions.M)
        error(['GAME_sampling ERROR: The field M of GAMEoptions ' ...
            'should be an integer with number of consecutive trials ' ...
            'of each estimator']);
    end
end
if isfield(GAMEoptions,'steps')
    if ischar(GAMEoptions.steps)
        error(['GAME_sampling ERROR: The field steps of GAMEoptions ' ...
            'should be an integer with number of steps used with ' ...
            'geomtric/optimal bridge estimators']);
    end
end

% Check DREAM input data structures
if ~isstruct(DREAMPar)
    error(['GAME_sampling ERROR: input argument DREAMPar should ' ...
        'be a structure with fields']);
end

% Check DREAM input data structures
if ~isstruct(Par_info)
    error(['GAME_sampling ERROR: input argument Par_info should ' ...
        'be a structure with fields']);
end

% Check DREAM input data structures
if ~isstruct(Meas_info)
    error(['GAME_sampling ERROR: input argument Meas_info should ' ...
        'be a structure with fields']);
end

% Check DREAM input data structures
if ~isstruct(options)
    error(['GAME_sampling ERROR: input argument options should ' ...
        'be a structure with fields']);
end

% Which fields has the user specified?
field_names = fieldnames(DREAMPar);

% Now make sure that all strings are lower case
for i = 1 : numel(field_names)
    eval(char(strcat('DREAMPar.',field_names(i),[' = ' ...
        'lower(DREAMPar.'],field_names(i),');')));
end

% Which fields has the user specified?
field_names = fieldnames(Par_info);

% Now make sure that all strings are lower case
for i = 1 : numel(field_names)
    if ~iscell( eval(char(strcat('Par_info.',field_names(i)))))
        eval(char(strcat('Par_info.',field_names(i),[' = ' ...
            'lower(Par_info.'],field_names(i),');')));
    else
        % Cell string with numbers; fine as is
    end
end

% Which fields has the use specified?
field_names = fieldnames(options);

% Now make sure that all strings are lower case
for i = 1 : numel(field_names)
    eval(char(strcat('options.',field_names(i),[' = ' ...
        'lower(options.'],field_names(i),');')));
end

% <><><><><><><><><><><><><><><> Func_name <><><><><><><><><><><><><><><><>
if isempty(Func_name)
    % ERROR -- Func_name is not defined
    error(['GAME_sampling ERROR: The variable Func_name has to be ' ...
        'defined as string (between quotes)']);
end

if isnumeric(Func_name)
    % ERROR -- Func_name is not properly defined
    error(char(strcat(['GAME_sampling ERROR: The variable Func_name ' ...
        'is defined as numerical value'],{' '}, ...
        num2str(Func_name),{' '},[': This should be a string ' ...
        '(between quotes) with name of MATLAB model script (.m file)'])));
end

if ischar(Func_name)
    % Check for Func_name.m and remove blank fields (.m and .mlx)
    file_m = evalc(char(strcat('ls ',{' '},Func_name,'.m')));
    file_m = strtrim(file_m);
    file_mlx = evalc(char(strcat('ls ',{' '},Func_name,'.mlx')));
    file_mlx = strtrim(file_mlx);
    % Now compare (.m and .mlx extensions)
    if strcmp(file_m,strcat(Func_name,'.m')) || ...
            strcmp(file_mlx,strcat(Func_name,'.mlx'))
        if strcmp(file_m,strcat(Func_name,'.m'))
            % Correct - function .m file exists
            char(strcat(['eDREAM_package CHECK: Function specified by ' ...
                'user in variable Func_name exists in directory:'], ...
                {' '},pwd))
        elseif strcmp(file_mlx,strcat(Func_name,'.mlx'))
            % Correct - Live script .mlx exists
            char(strcat(['eDREAM_package CHECK: Live script ' ...
                'specified by ' ...
                'user in variable Func_name exists in directory:'], ...
                {' '},pwd))
        end
    else
        % ERROR -- Func_name is not found in directory
        error(char(strcat(['eDREAM_package ERROR: Function (or Live ' ...
            'Script)' ...
            'specified by user in variable Func_name (CASE-SENSITIVE!) ' ...
            'cannot be found in directory:'],{' '},pwd)));
    end

end

% <><><><><><><><><><><><><><><> DREAMPar <><><><><><><><><><><><><><><><><
if ~isfield(DREAMPar,'d')
    error(['GAME_sampling ERROR: Field ''d'' of structure DREAMPar ' ...
        'undefined']);
end

if ( DREAMPar.d <= 0 )
    % ERROR -- dimensionality should be larger than zero
    error(['GAME_sampling ERROR: Number of parameters should be ' ...
        'integer and larger than zero -> Set DREAMPar.d >= 1']);
end

if ~ismember(DREAMPar.lik,[1 2 11:18 21])
    % Cannot use ML estimation!
    error(['GAME_sampling ERROR: Cannot compute marginal likelihood ' ...
        'for this type of (informal) likelihood/posterior density ' ...
        'function!!']);
end

if ( DREAMPar.lik == 1 || DREAMPar.lik == 2 )
    if ~isempty(fieldnames(Meas_info))
        warning(char(strcat(['GAME_sampling WARNING: No need to ' ...
            'specify fields of structure Meas_info - not used with' ...
            ' likelihood function'],{' '},num2str(DREAMPar.lik), ...
            {' '},'or DREAMPar.lik = ',{' '},num2str(DREAMPar.lik),'\n')));
    end
end

if ( DREAMPar.lik == 12 || DREAMPar.lik == 13 || DREAMPar.lik == 16 )
    if ~isfield(Meas_info,'Sigma')
        % ERROR -- Meas_info.Sigma needs to be specified!!
        error(['GAME_sampling ERROR: Meas_info.Sigma needs to be ' ...
            'defined either an inline function or one (homoscedastic) ' ...
            'or multiple (heteroscedastic) numerical values!!']);
    end
end

if ( DREAMPar.lik == 13 || DREAMPar.lik == 14 || DREAMPar.lik == 17 )
    warning(char(['GAME_sampling WARNING: Likelihood function ' ...
        'selected that contains nuisance variables: Please make ' ...
        'sure that setup is correct! (see DREAM manual)\n']));
end

if ismember(DREAMPar.lik,11:18)  % 31 - 34 not relevant
    if ~isfield(Meas_info,'Y')
        error(['GAME_sampling ERROR: Field ''Y'' of structure ' ...
            'Meas_info has to be defined (stores training ' ...
            'observations)!!']);
    else
        if isempty(Meas_info.Y)
            error(['GAME_sampling ERROR: Field ''Y'' of structure ' ...
                'Meas_info has to contain at least one value ' ...
                '(stores training observations)!!']);
        end
    end
end

if ismember(DREAMPar.lik,11:18) % 31 - 34 not relevant
    if size(Meas_info.Y,2) > 1
        error(['GAME_sampling ERROR: Field ''Y'' of structure ' ...
            'Meas_info has to be a column vector)!!']);
    elseif size(Meas_info.Y,1) == 1
        warning(char(['GAME_sampling WARNING: Only a single ' ...
            'calibration data observation is used/stored in ' ...
            'Meas_info.Y \n']));
    end
end

if ( DREAMPar.lik == 18 )
    if ~isfield(Meas_info,'C')
        error(['DREAM_PACKAGE ERROR: Field ''C'' of structure ' ...
            'Meas_info has to be defined for likelihood 18: ' ...
            'Meas_info.C stores n x n covariance matrix ' ...
            'measurement errors!!']);
    else
        [m,n] = size(Meas_info.C);
        if m ~= n
            error(['DREAM_PACKAGE ERROR: Field ''C'' of structure ' ...
                'Meas_info must be square as it is the covariance ' ...
                'matrix of the measurement data errors!!']);
        end
        if m ~= numel(Meas_info.Y)
            error(['DREAM_PACKAGE ERROR: Field ''C'' of structure ' ...
                'Meas_info must have a similar number of rows/columns ' ...
                'as length of data vector in Meas_info.Y!!']);
        end
    end
end

if ismember(DREAMPar.lik,21) % 22 - 23 not relevant
    % Limits of acceptability
    if ~isfield(Meas_info,'S')
        % ERROR -- Meas_info.S does not exist!!
        error(['GAME_sampling ERROR: Field ''S'' of structure ' ...
            'Meas_info has to be defined (stores simulated summary ' ...
            'metrics / observations)!!']);
    end
    if isfield(Meas_info,'S')
        if isempty(Meas_info.S)
            % ERROR -- Meas_info.S does not exist!!
            error(['GAME_sampling ERROR: Field ''S'' of structure ' ...
                'Meas_info has to be defined (stores simulated ' ...
                'summary metrics / observations)!!']);
        end
    end
end

if ismember(DREAMPar.lik,21) % 22 - 23 cannot be used!
    if size(Meas_info.S,2) > 1
        error(['GAME_sampling ERROR: Field ''S'' of structure ' ...
            'Meas_info has to be a column vector)!!']);
    elseif size(Meas_info.S,1) == 1
        warning(char(['GAME_sampling WARNING: Only a single ' ...
            'summary metric value is used/stored in Meas_info.S \n']));
    end
end

if ismember(DREAMPar.lik,[ 1:2 11:17 ])
    if isfield(options,'epsilon')
        warning(char(['GAME_sampling WARNING: No need to define ' ...
            'field ''epsilon'' of structure options: Not used - ' ...
            'only for ABC/LOA approaches with likelihood function 21 \n']));
    end
    if isfield(options,'rho')
        warning(char(['GAME_sampling WARNING: No need to define ' ...
            'field ''rho'' of structure options: Not used - ' ...
            'only for ABC/LOA approaches with likelihood function 21 \n']));
    end
    if isfield(Meas_info,'S')
        warning(char(['GAME_sampling WARNING: No need to define ' ...
            'field ''S'' of structure Meas_info: Not used - only ' ...
            'for ABC/LOA approaches with likelihood function 21 \n']));
    end
end

if ismember(DREAMPar.lik,21) % 22 cannot be used with ML
    % Warning that multiple epsilon values can be used
    if isfield(options,'epsilon')
        if (numel(Meas_info.S) ~= numel(options.epsilon)) ...
                && ( numel(options.epsilon) > 1 )
            % ERROR -- Meas_info.Sigma incorrect length!!
            error(['GAME_sampling ERROR: Number of elements of ' ...
                '''epsilon'' does not match that of Meas_info.S!!']);
        end
        if numel(options.epsilon) == 1
            warning(char(['GAME_sampling WARNING: If so desired you ' ...
                'can use a different value of ''epsilon'' for each ' ...
                'summary metric stored in Meas_info.S - just define ' ...
                'options.epsilon as vector \n']));
        end
    else
        warning(char(['GAME_sampling WARNING: ABC approach: Default ' ...
            'value of ''epsilon'' will be used (=0.025) and applied ' ...
            'to all summary statistics stored in Meas_info.S \n']));
    end
end

%<><><><><><><><><><><><><><><><> Par_info <><><><><><><><><><><><><><><><>

% If field boundhandling not defined --> define no use
if ~isfield(Par_info,'boundhandling'), Par_info.boundhandling = 'none'; end
% Check: do we sample in normalized [0-1] space, or not?
if isfield(Par_info,'norm')
    if ~isreal(Par_info.norm)
        error(char(strcat(['GAME_sampling ERROR: Par_info.norm' ...
            ' should be a scalar -> Define Par_info.norm = 0 or 1\n'])));
    elseif ~ismember(Par_info.norm,[0 1])
        error(char(strcat(['GAME_sampling ERROR: Par_info.norm' ...
            ' should be zero or one -> Define Par_info.norm = 0 or 1\n'])));
    end
    if Par_info.norm == 1
        if ~isfield(Par_info,'min')
            error(['GAME_sampling ERROR: Parameter normalization is used ' ...
                'but minimum parameter values not defined -> Set Par_info.min!!']);
        end
        if ~isfield(Par_info,'max')
            error(['GAME_sampling ERROR: Parameter normalization is used but' ...
                ' maximum parameter values not defined -> Set Par_info.max!!']);
        end
    end
elseif ~isfield(Par_info,'norm')
    Par_info.norm = 0;
end
if ~any(strcmp(Par_info.boundhandling,{'fold','bound','reflect',...
        'none','reject'}))
    evalstr = char(['GAME_sampling WARNING: Unknown boundary' ...
        ' handling method -> Define Par_info.boundhandling = ' ...
        '''fold''/''bound''/''reflect''/''reject''/''none''!!\n']);
    % Now print warning to screen and to file
    fprintf(evalstr); fprintf(fid,evalstr);
    % Continue without boundary handling
    evalstr = char(['GAME_sampling WARNING: Treating as unbounded ' ...
        'parameter space!!\n']);
    % Now print warning to screen and to file
    fprintf(evalstr); fprintf(fid,evalstr);
end
if any(strcmp(Par_info.boundhandling,{'fold','bound','reflect','reject'}))
    if ~isfield(Par_info,'min')
        error(['GAME_sampling ERROR: Boundary handling is used ' ...
            'but minimum parameter values not defined -> Set Par_info.min!!']);
    end
    if ~isfield(Par_info,'max')
        error(['GAME_sampling ERROR: Boundary handling is used but' ...
            ' maximum parameter values not defined -> Set Par_info.max!!']);
    end
end
% ERROR if Par_info.min or Par_info.max not equal to DREAMPar.d
if ( numel(Par_info.min) ~= DREAMPar.d )
    error(char(strcat(['GAME_sampling ERROR: Number of elements ' ...
        'of field ''min'' of structure Par_info should be ' ...
        'equal to '],{' '},num2str(DREAMPar.d),'!!')));
end
% ERROR if Par_info.min or Par_info.max not equal to DREAMPar.d
if ( numel(Par_info.max) ~= DREAMPar.d )
    error(char(strcat(['GAME_sampling ERROR: Number of elements ' ...
        'of field ''max'' of structure Par_info should be equal' ...
        ' to '],{' '},num2str(DREAMPar.d),'!!')));
end
% % Remove dummy variable to make things easier
% if strcmp(Par_info.boundhandling,'none')
%     Par_info = rmfield(Par_info,'boundhandling');
% end
% JAV: Jan. 2022: Check whether user specified names of parameters
if isfield(Par_info,'names')
    if numel(Par_info.names) ~= DREAMPar.d
        error(char(strcat(['GAME_sampling ERROR: Number of elements' ...
            ' of field ''names'' of structure Par_info should be ' ...
            'equal to '],{' '},num2str(DREAMPar.d),'!!')));
    end
end

% <><><><><><><><><><><><><><><> Meas_info <><><><><><><><><><><><><><><><>
if isfield(Meas_info,'Sigma')
    % Check content of Sigma
    if isreal(Meas_info.Sigma)  % --> true: then not inline function
        if (numel(Meas_info.Sigma) ~= numel(Meas_info.Y)) ...
                && ( numel(Meas_info.Sigma) > 1 )
            % ERROR -- Meas_info.Sigma incorrect length!!
            error(['GAME_sampling ERROR: Length of Meas_info.Sigma ' ...
                'is not equal to that of the observations stored ' ...
                'in Meas_info.Y!!']);
        end
        % Now check this as well
        if any(Meas_info.Sigma < 0)
            fprintf(['GAME_sampling WARNING: One or more entries of ' ...
                'Meas_info.Sigma is negative - we use absolute values\n']);
            Meas_info.Sigma = abs(Meas_info.Sigma);
        end
        if any(Meas_info.Sigma == 0)
            fprintf(['GAME_sampling WARNING: One or more entries of ' ...
                'Meas_info.Sigma are zero - we use smallest positive ' ...
                'value of entered values, otherwise set those entries ' ...
                'of Meas_info.Sigma equal to 1e-3\n']);
            id = find(Meas_info.Sigma == 0);
            i_std = Meas_info.Sigma(Meas_info.Sigma > 0);
            if ~isempty(i_std)
                Meas_info.Sigma(id) = min(i_std);
            else
                Meas_info.Sigma(id) = 1e-3;
            end
        end
    end
end

% <><><><><><><><><><><><><><><>< options ><><><><><><><><><><><><><><><><>
if isfield(options,'epsilon')
    if ( options.epsilon <= 0 )
        error(['GAME_sampling ERROR: Value of ''epsilon'' of ' ...
            'structure options should be larger than zero ' ...
            '-> Set options.epsilon > 0 (default 0.025)']);
    end
end
if isfield(options,'rho')
    if isreal(options.rho)
        error(['GAME_sampling ERROR: Field ''rho'' of structure ' ...
            'options should be defined as inline function ' ...
            '-> Default: options.rho = @(X,Y) abs(X-Y);']);
    end
end
% Check content of each field of structure options
name = fieldnames(options);
% Warning if content not equal to yes or no
for j = 1 : numel(name)
    % Get content of respective field of options
    F = options.(char(name(j)));
    % now check content
    if ~any(strcmp(name(j),{'epsilon','rho','burnin'}))
        if ~any(strcmp(F,{'yes','no'}))
            % ERROR -- content field name(j) of options should be yes or no
            error(char(strcat('GAME_sampling ERROR: Field', ...
                {' '},'''',name(j),'''',{' '},['of structure ' ...
                'options should be set equal to ''yes'' or ''no'''])));
        end
    end
end

% If parallel computing is requested - check whether toolbox is available
if isfield(options,'parallel')
    % If parallel = 'yes' -> check whether user has the distributed
    % computing toolbox
    if strcmp(options.parallel,'yes')
        % Now check whether Distributed Computing toolbox available
        toolbox = license('checkout','Distrib_Computing_Toolbox');
        % If toolbox is zero then list comment
        if ( toolbox == 0 )
            evalstr = char(['GAME_sampling WARNING: Distributed' ...
                ' Computing toolbox not available -> code resorts ' ...
                'to default: options.parallel = ''no'' \n']);
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
            % And update field parallel of options
            options.parallel = 'no';
        end
    end
else
    options.parallel = 'no';
end
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

end
