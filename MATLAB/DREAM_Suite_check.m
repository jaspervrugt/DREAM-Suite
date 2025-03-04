function [DREAMPar,Par_info,options] = DREAM_Suite_check(method, ...
    Func_name,DREAMPar,Par_info,Meas_info,options)
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
% This function verifies the input arguments of DREAM-Suite defined       %
% by the user                                                             %
%                                                                         %
% SYNOPSIS: [DREAMPar,Par_info,options] = DREAM_Suite_check(method,...    %
%               Func_name,DREAMPar,Par_info,Meas_info,options)            %
%                                                                         %
% © Written by Jasper A. Vrugt, Feb 2007                                  %
% Los Alamos National Laboratory 			        	                  %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Derive current time and set deadline
% % deadline = datenum('28-Feb-2025');
% % % Now check whether this is a trial version or not
% % if ( deadline - now ) < 0
% %     % ERROR -- trial version ended
% %     error('DREAM_Suite ERROR: Trial version of DREAM V3.0 has ended');
% % end
% %
% open an output file with warnings

fid = fopen('warning_file.txt','w+');
fprintf(fid,'-------------- DREAM_Suite warning file --------------\n');
% Check DREAM input data structures
if ~isstruct(DREAMPar)
    error(['DREAM_Suite ERROR: input argument DREAMPar ' ...
        'should be a structure with fields']);
end
% Check DREAM input data structures
if ~isstruct(Par_info)
    error(['DREAM_Suite ERROR: input argument Par_info ' ...
        'should be a structure with fields']);
end
% Check DREAM input data structures
if ~isstruct(Meas_info)
    error(['DREAM_Suite ERROR: input argument Meas_info ' ...
        'should be a structure with fields']);
end
% Check DREAM input data structures
if ~isstruct(options)
    error(['DREAM_Suite ERROR: input argument options ' ...
        'should be a structure with fields']);
end
% All strings must be lower case, except parameter names
str = {'DREAMPar','Par_info','options'};
for z = 1:numel(str)
    vrbl = eval(char(str(z))); fd_names = fieldnames(vrbl);
    for i = 1:numel(fd_names)
        if ischar(vrbl.(char(fd_names(i))))
            evalstr = strcat(str(z),['.(char(fd_names(i))) = ' ...
                'lower('],str(z),'.(char(fd_names(i))));');
            eval(char(evalstr))
        end
    end
end

% <><><><><><><><><><><><><><><> Func_name <><><><><><><><><><><><><><><><>
if isempty(Func_name)
    % ERROR -- Func_name is not defined
    error(['DREAM_Suite ERROR: The variable Func_name has ' ...
        'to be defined as string (between quotes)']);
end
if isnumeric(Func_name)
    % ERROR -- Func_name is not properly defined
    error(char(strcat(['DREAM_Suite WARNING: The variable ' ...
        'Func_name is defined as numerical value'],{' '}, ...
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
            char(strcat(['DREAM_Suite CHECK: Function specified by ' ...
                'user in variable Func_name exists in directory:'], ...
                {' '},pwd))
        elseif strcmp(file_mlx,strcat(Func_name,'.mlx'))
            % Correct - Live script .mlx exists
            char(strcat(['DREAM_Suite CHECK: Live script ' ...
                'specified by ' ...
                'user in variable Func_name exists in directory:'], ...
                {' '},pwd))
        end
    else
        % ERROR -- Func_name is not found in directory
        error(char(strcat(['DREAM_Suite ERROR: Function (or Live ' ...
            'Script)' ...
            'specified by user in variable Func_name (CASE-SENSITIVE!) ' ...
            'cannot be found in directory:'],{' '},pwd)));
    end

end

% <><><><><><><><><><><><><><><> DREAMPar <><><><><><><><><><><><><><><><><
if strcmp(method,'dream') || strcmp(method,'dream_d')
    if ~isfield(DREAMPar,'N')
        error(['DREAM/DREAM_{D} ERROR: Field ''N'' of structure ' ...
            'DREAMPar undefined']);
    end
end
if ~isfield(DREAMPar,'d')
    error(['DREAM_Suite ERROR: Field ''d'' of structure ' ...
        'DREAMPar undefined']);
end
if ~isfield(DREAMPar,'T')
    error(['DREAM_Suite ERROR: Field ''T'' of structure ' ...
        'DREAMPar undefined']);
end
if ~isfield(DREAMPar,'lik')
    error(['DREAM_Suite ERROR: Field ''lik'' of structure ' ...
        'DREAMPar undefined']);
end
if ( DREAMPar.d <= 0 )
    % ERROR -- dimensionality should be larger than zero
    error(['DREAM_Suite ERROR: Number of parameters should ' ...
        'be integer and larger than zero -> Set DREAMPar.d >= 1']);
end
if isfield(DREAMPar,'delta')
    if DREAMPar.delta <= 0
        % ERROR -- not enough chains pairs for sampling
        error(['DREAM_Suite ERROR: Number of chains pairs ' ...
            'used for offspring should be integer and larger ' ...
            'than zero -> Use at least DREAMPar.delta \in [1,5] ' ...
            '(default: 1)']);
    elseif DREAMPar.delta > 10
        evalstr = char(strcat(['DREAM_Suite WARNING: Field ' ...
            '''delta'' of structure DREAMPar set rather large ' ...
            '-> recommend to use values of delta in [1,10]\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
    end
end
switch method
    case {'dream','dream_d'}
        if isfield(DREAMPar,'delta')
            delta = DREAMPar.delta;
        else
            delta = 3;
        end
        % Calculate minimum number of chains
        N_min = max( [ 8 ; 2*delta + 1 ; ceil(DREAMPar.d/2) ]);
        if DREAMPar.N < N_min
            % ERROR -- not enough chains for sampling proposals
            error(char(strcat(['DREAM/DREAM_{D} ERROR: ' ...
                'Inufficient number of chains with delta = '], ...
                {' '},num2str(delta),{' '},['-> Use at least ' ...
                'DREAMPar.N = '],{' '},num2str(N_min),{' '},'chains')));
        end
        % Calculate maximum number of chains
        N_max = max( [ 15 ; 2*delta + 1 ; ceil(DREAMPar.d) ]);
        if DREAMPar.N > N_max
            evalstr = char(strcat(['DREAM/DREAM_{D} WARNING: ' ...
                'Field ''N'' of structure DREAMPar set rather large ' ...
                '--> DREAMPar.N = '],{' '},num2str(N_max),{' '}, ...
                'more than sufficient\n'));
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
        end
    case {'dream_zs','dream_dzs','dream_kzs'}
        if isfield(DREAMPar,'delta')
            delta = DREAMPar.delta;
        else
            delta = 1;
        end % default
        if isfield(DREAMPar,'m0')
            m0_min = max ( 10 * DREAMPar.d , max ( 3 * DREAMPar.N , ...
                2 * delta * DREAMPar.N ) );
            if DREAMPar.m0 < m0_min
                error(char(strcat(['' ...
                    'DREAM_{(ZS)}/DREAM_{(DZS)}/DREAM_{(KZS)} ERROR: ' ...
                    'Field ''m0'' of structure DREAMPar should be ' ...
                    'integer and set to at least'],{' '}, ...
                    num2str(m0_min),{' '},'(default: 10 x DREAMPar.d)')));
            elseif DREAMPar.m0 > 50*DREAMPar.d
                evalstr = char(strcat(['' ...
                    'DREAM_{(ZS)}/DREAM_{(DZS)}/DREAM_{(KZS)} WARNING: '...
                    'Field ''m0'' of structure DREAMPar set rather ' ...
                    'large --> recommend to use DREAMPar.m0 = ' ...
                    '10 x DREAMPar.d = '],{' '},num2str(10*DREAMPar.d), ...
                    {' '},'\n'));
                % Now print warning to screen and to file
                fprintf(evalstr); fprintf(fid,evalstr);
            end
        end
        if isfield(DREAMPar,'k')
            if DREAMPar.k < 1
                error(['DREAM_{(ZS)}/DREAM_{(DZS)}/DREAM_{(KZS)} ' ...
                    'ERROR: Field ''k'' of structure DREAMPar should ' ...
                    'be integer and set to at least 1 ' ...
                    '-> DREAMPar.k in [1,20] (default: DREAMPar.k = 10)']);
            elseif DREAMPar.k > 25
                evalstr = char(strcat(['' ...
                    'DREAM_{(ZS)}/DREAM_{(DZS)}/DREAM_{(KZS)} WARNING:' ...
                    ' Field ''k'' of structure DREAMPar set ' ...
                    'rather large ' ...
                    '--> recommend to use DREAMPar.k = 10 (default)\n']));
                % Now print warning to screen and to file
                fprintf(evalstr); fprintf(fid,evalstr);
            end
        end
end
if DREAMPar.T < 2
    % ERROR -- not enough generations
    error(['DREAM_Suite ERROR: Number of generations smaller ' ...
        'than one -> Set at least DREAMPar.T = 2']);
elseif DREAMPar.T > 1e6
    evalstr = char(strcat(['DREAM_Suite WARNING: Field ''T'' ' ...
        'of structure DREAMPar set rather large']));
    % Now print warning to screen and to file
    fprintf(evalstr); fprintf(fid,evalstr);
end
if isfield(DREAMPar,'nCR')
    if DREAMPar.nCR <= 0
        % ERROR -- at least one crossover value
        error(['DREAM_Suite ERROR: Number of crossover values used ' ...
            'for offspring should be integer and larger than zero ' ...
            '-> Use at least DREAMPar.nCR = 1 (default: 3)'])
    elseif ( DREAMPar.nCR < max ( DREAMPar.d / 5 , 1 ) )
        evalstr = char(strcat(['DREAM_Suite SUGGESTION: Number of ' ...
            'crossover values of'],{' '},num2str(DREAMPar.nCR), ...
            {' '},['defined in field ''nCR'' of structure DREAMpar ' ...
            'is rather small for DREAMPar.d = '],{' '}, ...
            num2str(DREAMPar.d),{' '},'(default: 3)\n'));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
        % Now continue
        evalstr = char(strcat(['DREAM_Suite SUGGESTION: Try another ' ...
            'trial later using DREAMPar.nCR = ' ...
            'max(floor(DREAMPar.d/5),1) = '],{' '}, ...
            num2str(max(floor(DREAMPar.d/5),1)),{' '}, ...
            ['(rule of thumb) as this should enhance the sampling ' ...
            'and convergence speed of DREAM\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
    elseif ( DREAMPar.nCR >= 1/4 * DREAMPar.d )
        evalstr = char(strcat(['DREAM_Suite WARNING: Number ' ...
            'of crossover values of'],{' '},num2str(DREAMPar.nCR), ...
            {' '},['defined in field ''nCR'' of structure DREAMpar ' ...
            'is rather large for DREAMPar.d = '],{' '}, ...
            num2str(DREAMPar.d),{' '},'(default: 3)\n'));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
        % Now continue with suggestion
        evalstr = char(strcat(['DREAM_Suite SUGGESTION: Start ' ...
            'another trial using DREAMPar.nCR = ' ...
            'max(floor(DREAMPar.d/5),1) = '],{' '}, ...
            num2str(max(floor(DREAMPar.d/5),1)),{' '}, ...
            ['(rule of thumb) as this should enhance the sampling ' ...
            'and convergence speed of DREAM\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
    end
end
if isfield(DREAMPar,'thinning')
    if ( DREAMPar.thinning < 1 )
        % ERROR -- thinning parameter should be positive
        error(['DREAM_Suite ERROR: Thinning parameter should ' ...
            'be integer and larger than zero ' ...
            '-> Set DREAMPar.thinning >= 1 (default: 1)']);
    elseif ( DREAMPar.thinning > 20 )
        evalstr = char(strcat(['DREAM_Suite WARNING: Field ' ...
            '''thinning'' of structure DREAMPar set rather large ' ...
            '-> recommend to use DREAMPar.thinning in [1,20]\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
    end
end
if isfield(DREAMPar,'beta0')
    if ( DREAMPar.beta0 <= 0 )
        % ERROR -- jump rate multiplier should be larger than zero
        error(['DREAM_Suite ERROR: Multiplier of jump rate ' ...
            'should be larger than zero ' ...
            '-> Set DREAMPar.beta0 > 0 (default: 1)']);
    elseif ( DREAMPar.beta0 > 3 )
        evalstr = char(strcat(['DREAM_Suite WARNING: Field ' ...
            '''beta0'' of structure DREAMPar set rather large ' ...
            '-> recommend to use values of beta0 in [0.1,3]\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
    end
end
if isfield(DREAMPar,'GLUE')
    if ( DREAMPar.GLUE <= 0 )
        % ERROR -- GLUE likelihood variable larger than zero
        error(['DREAM_Suite ERROR: Likelihood variable of ' ...
            '''GLUE'' should be larger than zero ' ...
            '-> Set DREAMPar.GLUE > 0 ']);
    elseif ( DREAMPar.GLUE > 1000 )
        evalstr = char(strcat(['DREAM_Suite WARNING: Field ' ...
            '''GLUE'' of structure DREAMPar set rather large\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
    end
end
if isfield(DREAMPar,'lambda')
    if ( DREAMPar.lambda <= 0 )
        % ERROR -- lambda should be positive
        error(['DREAM_Suite ERROR: Value of ''lambda'' should ' ...
            'be larger than zero -> Set DREAMPar.lambda > 0 ' ...
            '(default: 0.05)']);
    elseif ( DREAMPar.lambda > 1 )
        evalstr = char(strcat(['DREAM_Suite WARNING: Field ' ...
            '''lambda'' of structure DREAMPar set rather large ' ...
            '-> recommend to use DREAMPar.lambda in [0.01,0.25]\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
    end
end
if isfield(DREAMPar,'zeta')
    if ( DREAMPar.zeta <= 0 )
        % ERROR -- zeta should be positive
        error(['DREAM_Suite ERROR: Value of ''zeta'' should be ' ...
            'larger than zero -> Set DREAMPar.zeta > 0 (default: 1e-12)']);
    elseif ( DREAMPar.zeta > 0.1 )
        evalstr = char(strcat(['DREAM_Suite WARNING: Field ''zeta'' ' ...
            'of structure DREAMPar set rather large -> recommend to ' ...
            'use DREAMPar.zeta in [1e-20,1e-3]\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
    end
end
if ~ismember(DREAMPar.lik,[1:2 11:17 21:23 31:34 44 45 61 62 99])
    % Unknown built-in likelihood function
    error(['DREAM_Suite ERROR: Unknown choice of likelihood ' ...
        'function -> Select DREAMPar.lik = ' ...
        '{1,2,11,12,13,14,15,16,17,21,22,23,31,32,33,34,44,45} ' ...
        'or 99 (own likelihood)']);
end
if isfield(DREAMPar,'p_unit_gamma')
    if ( DREAMPar.p_unit_gamma < 0 ) || ( DREAMPar.p_unit_gamma > 1)
        % ERROR -- unit jump rate probability between 0 and 1
        error(['DREAM_Suite ERROR: Probability of unit jump rate ' ...
            'should be between 0 and 1 ' ...
            '-> Set DREAMPar.p_unit_gamma in [0,1] (default: 0.2)']);
    end
end
if strcmp(method,'dream_zs') || strcmp(method,'dream_dzs')
    if isfield(DREAMPar,'psnooker')
        if ( DREAMPar.psnooker < 0 ) || ( DREAMPar.psnooker > 1)
            % ERROR -- snooker jump probability must be between 0 and 1
            error(['DREAM_{(ZS)}/DREAM_{(DZS)}/DREAM_{(KZS)} ERROR: ' ...
                'Probability of snooker jump should be between 0 and 1 ' ...
                '-> Set DREAMPar.psnooker in [0,1] (default: 0.1)']);
        end
    end
end
if strcmp(method,'dream_kzs')
    if ~ismember(DREAMPar.lik,[11:17 31:34 44:45])
        % ERROR -- DREAM_KZS requires simulated model output
        error(['DREAM_{(KZS)} ERROR: Kalman jump cannot be used: ' ...
            'model output is a (log)-density and not a simulation']);
    end
    if ( ~isfield(Meas_info,'R') || isempty(Meas_info.R) ) ...
            && ~isfield(Meas_info,'Sigma')
        % WARNING -- DREAM_KZS will not consider the measurement error 
        % in the Kalman jump
        warning(['DREAM_{(KZS)} WARNING: Field ''R'' and field' ...
            ' ''Sigma'' of structure Meas_info are not specified ' ...
            '-> Kalman jump will assume that data are observed' ...
            ' without a measurement error']);
    end
    if isfield(Meas_info,'R')
        % June 18, 2018: This statement must move to later time as 
        % Meas_info.Y has not been checked whether it exists
        if size(Meas_info.R,1) ~= numel(Meas_info.Y)
            % ERROR -- Measurement error covariance matrix has to 
            % be square and contain Meas_info.n rows
            error(['DREAM_{(KZS)} ERROR: Measurement error ' ...
                'covariance matrix, Meas_info.R, has to have ' ...
                'similar number of rows as length of data vector']);
        end
        if size(Meas_info.R,2) ~= numel(Meas_info.Y)
            % ERROR -- Measurement error covariance matrix has to 
            % be square and contain Meas_info.n columns
            error(['DREAM_{(KZS)} ERROR: Measurement error ' ...
                'covariance matrix, Meas_info.R, has to have ' ...
                'similar number of columns as length of data vector']);
        end
    end
    % Now check selection probabilities of proposal distributions
    if isfield(DREAMPar,'psnooker')
        if ( DREAMPar.psnooker < 0 ) || ( DREAMPar.psnooker > 0.4)
            % ERROR -- snooker jump probability must be between 0 and 0.4
            error(['DREAM_{(KZS)} ERROR: Probability of snooker ' ...
                'jump should be between 0 and 0.4 ' ...
                '-> Set DREAMPar.psnooker in [0,0.4] (default: 0.1)']);
        end
    end
    if isfield(DREAMPar,'pkalman')
        if ( DREAMPar.pkalman < 0 ) || ( DREAMPar.pkalman > 0.5)
            % ERROR -- kalman jump probability must be between 0 and 0.5
            error(['DREAM_{(KZS)} ERROR: Probability of kalman ' ...
                'jump should be between 0 and 0.5 ' ...
                '-> Set DREAMPar.pkalman in [0,1] (default: 0.4)']);
        end
    end
    if isfield(DREAMPar,'psnooker') && isfield(DREAMPar,'pkalman')
        if DREAMPar.psnooker + DREAMPar.pkalman >= 1
            % ERROR -- parallel direction jump probability will be zero or negative
            error(['DREAM_{(KZS)} ERROR: Sum of selection ' ...
                'probabilities of parallel direction and snooker ' ...
                'jump cannot be larger than unity']);
        end
    end
    if isfield(DREAMPar,'a_1')
        if DREAMPar.a_1 < 0
            % ERROR -- a_1 cannot be larger than one
            error(['DREAM_{(KZS)} ERROR: Value of a_1 cannot ' ...
                'be smaller than zero']);
        end
        if DREAMPar.a_1 >= 1
            % ERROR -- a_1 cannot be larger than one
            error(['DREAM_{(KZS)} ERROR: Value of a_1 cannot ' ...
                'exceed value of 0.5']);
        end
    end 
    if isfield(DREAMPar,'a_2')
        if DREAMPar.a_2 < 0
            % ERROR -- a_1 cannot be larger than one
            error(['DREAM_{(KZS)} ERROR: Value of a_2 cannot be ' ...
                'smaller than zero']);
        end
        if DREAMPar.a_2 >= 1
            % ERROR -- a_1 cannot be larger than one
            error(['DREAM_{(KZS)} ERROR: Value of a_2 cannot exceed ' ...
                'value of 0.5']);
        end
    end
    if isfield(DREAMPar,'a_1') && isfield(DREAMPar,'a_2')
        if DREAMPar.a_2 <= DREAMPar.a_1
            % ERROR -- a_1 cannot be larger than one
            error(['DREAM_{(KZS)} ERROR: Value of a_1 cannot exceed ' ...
                'or be equal to value of a_1']);
        end
    end
end
if isfield(DREAMPar,'adapt_pCR')
%    if isempty(strmatch(DREAMPar.adapt_pCR,{'no','yes'}))
    if ~any(strcmp(DREAMPar.adapt_pCR,{'no','yes'}))
        % ERROR -- yes or no adaptation selection prob. crossover values
        error(['DREAM_Suite ERROR: Adaptation of crossover values' ...
            ' should be yes (default) or no ' ...
            '-> Set DREAMPar.adapt_pCR = ''yes'' (default) or ''no''']);
    end
end
if strcmp(method,'dream') || strcmp(method,'dream_d')
    if isfield(DREAMPar,'outlier')
        if ~any(strcmp(DREAMPar.outlier,{'iqr','grubbs',...
                'peirce','chauvenet'}))
            % ERROR -- outlier detection test unknown
            error(['DREAM/DREAM_{D} ERROR: Unknown outlier detect test ' ...
                '-> Set DREAMPar.outlier = ''iqr'' (default) or ' ...
                '''grubbs'' or ''peirce'' or ''chauvenet'' ']);
        end
    end
    if isfield(DREAMPar,'outlier')
        if strcmp(DREAMPar.outlier,'peirce')
            if DREAMPar.N > 60
                % Warning -- Peirce will use r-values of DREAMPar.N = 60
                evalstr = char(strcat(['DREAM/DREAM_{D} WARNING: ' ...
                    'Peirce outlier detect test will use r-values ' ...
                    'assuming 60 different chains\n']));
                % Now print warning to screen and to file
                fprintf(evalstr); fprintf(fid,evalstr);
            end
        end
    end
end
if isfield(DREAMPar,'thinning')
    if ( DREAMPar.T / DREAMPar.thinning ) < 200
        evalstr = char(strcat(['DREAM_Suite WARNING: Cannot compute ' ...
            'within-chain convergence diagnostics due to inufficient ' ...
            'number of chain samples -> Use at least DREAMPar.T = '], ...
            {' '},num2str((200 * DREAMPar.thinning)),'\n'));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
        % Write hint to screen as well
        evalstr = char(['DREAM_Suite HINT: Thinning reduces length ' ...
            'of the sampled chain! -> ratio of DREAMPar.T and ' ...
            'DREAMPar.thinning should be larger than 200 for ' ...
            'within-chain diagnostics\n']);
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
    end
end
if ~isfield(DREAMPar,'thinning')
    if DREAMPar.T < 200
        evalstr = char(strcat(['DREAM_Suite WARNING: Cannot compute' ...
            ' within-chain convergence diagnostics due to inufficient' ...
            ' number of chain samples -> Use at least DREAMPar.T = '], ...
            {' '},num2str(200),'\n'));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
    end
end
if ( DREAMPar.lik == 1 || DREAMPar.lik == 2 )
    if ~isempty(fieldnames(Meas_info))
        evalstr = char(strcat(['DREAM_Suite WARNING: No need to ' ...
            'specify fields of structure Meas_info - not used with ' ...
            'likelihood function'],{' '},num2str(DREAMPar.lik), ...
            {' '},'or DREAMPar.lik = ',{' '},num2str(DREAMPar.lik),'\n'));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
    end
end
if ismember(DREAMPar.lik,12)
    if ~isfield(Meas_info,'Sigma')
        % ERROR -- Meas_info.Sigma needs to be specified!!
        error(['DREAM_Suite ERROR: Meas_info.Sigma needs to be' ...
            ' defined either as one (homoscedastic)' ...
            ' or multiple (heteroscedastic) numerical values!!']);
    end
end
if ismember(DREAMPar.lik,[ 13 14 16 17 44 45]) 
    evalstr = char(['DREAM_Suite WARNING: Likelihood function' ...
        ' selected that contains nuisance variables: Please make' ...
        ' sure that setup is correct! (see DREAM manual)\n']);
    % Now print warning to screen and to file
    fprintf(evalstr); fprintf(fid,evalstr);
end
if ismember(DREAMPar.lik,[11:17 31:34 44 45 99 ])
    if ~isfield(Meas_info,'Y')
        error(['DREAM_Suite ERROR: Field ''Y'' of structure' ...
            ' Meas_info has to be defined (stores calibration' ...
            ' observations)!!']);
    else
        if isempty(Meas_info.Y)
            error(['DREAM_Suite ERROR: Field ''Y'' of structure' ...
                ' Meas_info has to contain at least one value' ...
                ' (stores calibration observations)!!']);
        end
    end
end
if ismember(DREAMPar.lik,[11:17 31:34 44 45 99 ])
    if size(Meas_info.Y,2) > 1
        error(['DREAM_Suite ERROR: Field ''Y'' of structure' ...
            ' Meas_info has to be a column vector)!!']);
    elseif size(Meas_info.Y,1) == 1
        evalstr = char(['DREAM_Suite WARNING: Only a single' ...
            ' calibration data observation is used/stored in ' ...
            'Meas_info.Y \n']);
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
    end
end
if ( DREAMPar.lik == 52 )
    if ~isfield(Meas_info,'C')
        error(['DREAM_Suite ERROR: Field ''C'' of structure ' ...
            'Meas_info has to be defined for likelihood 52: ' ...
            'Meas_info.C stores n x n covariance matrix measurement ' ...
            'errors!!']);
    else
        [m,n] = size(Meas_info.C);
        if m ~= n
            error(['DREAM_Suite ERROR: Field ''C'' of structure' ...
                ' Meas_info must be square as it is the covariance' ...
                ' matrix of the measurement data errors!!']);
        end
        if m ~= numel(Meas_info.Y)
           error(['DREAM_Suite ERROR: Field ''C'' of structure' ...
               ' Meas_info must have a similar number of rows/columns' ...
               ' as length of data vector in Meas_info.Y!!']);
        end 
    end
end
if ismember(DREAMPar.lik,21:23)
    % Limits of acceptability
    if ~isfield(Meas_info,'S')
        % ERROR -- Meas_info.S does not exist!!
        error(['DREAM_Suite ERROR: Field ''S'' of structure' ...
            ' Meas_info has to be defined (stores simulated ' ...
            'summary metrics / observations)!!']);
    end
    if isfield(Meas_info,'S')
        if isempty(Meas_info.S)
            % ERROR -- Meas_info.S does not exist!!
            error(['DREAM_Suite ERROR: Field ''S'' of structure' ...
                ' Meas_info has to be defined (stores simulated' ...
                ' summary metrics / observations)!!']);
        end
    end
end
if ismember(DREAMPar.lik,21:23)
    if size(Meas_info.S,2) > 1
        error(['DREAM_Suite ERROR: Field ''S'' of structure' ...
            ' Meas_info has to be a column vector)!!']);
    elseif size(Meas_info.S,1) == 1
        evalstr = char(['DREAM_Suite WARNING: Only a single summary' ...
            ' metric value is used/stored in Meas_info.S \n']);
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
    end
end
if ismember(DREAMPar.lik,[ 1:2 11:17 31:34 44 45 99 ])
    if isfield(options,'epsilon')
        evalstr = char(['DREAM_Suite WARNING: No need to define' ...
            ' field ''epsilon'' of structure options: Not used - only' ...
            ' for ABC/LOA approaches with likelihood function 21-23 \n']);
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
    end
    if isfield(options,'rho')
        evalstr = char(['DREAM_Suite WARNING: No need to define ' ...
            'field ''rho'' of structure options: Not used - only for' ...
            ' ABC/LOA approaches with likelihood function 21-23 \n']);
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
    end
    if isfield(Meas_info,'S')
        evalstr = char(['DREAM_Suite WARNING: No need to define ' ...
            'field ''S'' of structure Meas_info: Not used - only ' ...
            'for ABC/LOA approaches with likelihood function 21-23\n']);
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
    end
end
if ismember(DREAMPar.lik,21:22)
    % Warning that multiple epsilon values can be used
    if isfield(options,'epsilon')
        if (numel(Meas_info.S) ~= numel(options.epsilon)) ...
                && ( numel(options.epsilon) > 1 )
            % ERROR -- Meas_info.S incorrect length!!
            error(['DREAM_Suite ERROR: Number of elements of ' ...
                '''epsilon'' does not match that of Meas_info.S!!']);
        end
        if numel(options.epsilon) == 1
            evalstr = char(['DREAM_Suite WARNING: If so desired you' ...
                ' can use a different value of ''epsilon'' for each ' ...
                'summary metric stored in Meas_info.S - just define ' ...
                'options.epsilon as vector \n']);
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
        end
    else
        evalstr = char(['DREAM_Suite WARNING: ABC approach: Default ' ...
            'value of ''epsilon'' will be used (=0.025) and applied ' ...
            'to all summary statistics stored in Meas_info.S \n']);
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
    end
end
if ( DREAMPar.lik == 23 )
    % Warning that multiple epsilon values can be used
    if isfield(options,'epsilon')
        if (numel(Meas_info.S) ~= numel(options.epsilon)) ...
                && ( numel(options.epsilon) > 1 )
            % ERROR -- Meas_info.S incorrect length!!
            error(['DREAM_Suite ERROR: Number of elements of ' ...
                '''epsilon'' does not match that of Meas_info.S!!']);
        end
        if numel(options.epsilon) == 1
            evalstr = char(['DREAM_Suite WARNING: Limits of ' ...
                'Acceptability - All observations use same LOA - ' ...
                'You can define ''options.epsilon'' as vector with ' ...
                'a different value for each observation\n']);
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
        end
    else
        evalstr = char(['DREAM_Suite WARNING: Limits of Acceptability' ...
            ' - Default value of ''epsilon'' will be used (=0.025) ' ...
            'and applied to all observations stored in Meas_info.S \n']);
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
    end
end

%<><><><><><><><><><><><><><><><> Par_info <><><><><><><><><><><><><><><><>
if ~isfield(Par_info,'initial')
    error(['DREAM_Suite ERROR: Initial sampling distribution not ' ...
        'defined -> Define Par_info.initial = ''latin'' or ''uniform'' ' ...
        'or ''normal'' or ''prior'' or ''user'' !!']);
end
if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    if ~isfield(Par_info,'steps')
        error(['DREAM_{(D)}/DREAM_{(DZS)} ERROR: The number of ' ...
            'intervals of each parameter (integers) need to be ' ...
            'provided -> Define Par_info.steps!!']);
    end
    if numel(Par_info.steps) ~= DREAMPar.d
        error(char(strcat(['DREAM_{(D)}/DREAM_{(DZS)} ERROR: The ' ...
            'number of elements of Par_info.steps needs to be ' ...
            'equal to '],{' '},num2str(DREAMPar.d),'!!')));
    end
    if any(Par_info.steps < 3)
        error(['DREAM_{(D)}/DREAM_{(DZS)} ERROR: The values of ' ...
            'Par_info.steps need to be integers and larger than 2!!']);
    end
end
if ~any(strcmp(Par_info.initial,{'latin','uniform','normal',...
        'prior','user'}))
    error(['DREAM_Suite ERROR: Initial sampling distribution ' ...
        'unknown -> Set Par_info.initial = ''latin'' or ''uniform''' ...
        ' or ''normal'' or ''prior'' or ''user'' !!']);
end
if strcmp(Par_info.initial,'latin')
    % ERROR -- if lhs is used -> requires explicit parameter ranges
    if ~isfield(Par_info,'min')
        error(['DREAM_Suite ERROR: Latin hypercube sampling ' ...
            'selected but minimum parameter values not defined ' ...
            '-> Set Par_info.min!!']);
    end
    if ~isfield(Par_info,'max')
        error(['DREAM_Suite ERROR: Latin hypercube sampling ' ...
            'selected but maximum parameter values not defined ' ...
            '-> Set Par_info.max!!']);
    end
end
if strcmp(Par_info.initial,'uniform')
    % ERROR -- if uniform initial sampling is used -> requires 
    % explicit parameter ranges
    if ~isfield(Par_info,'min')
        error(['DREAM_Suite ERROR: Uniform initial sampling ' ...
            'selected but minimum parameter values not defined ' ...
            '-> Set Par_info.min!!']);
    end
    if ~isfield(Par_info,'max')
        error(['DREAM_Suite ERROR: Uniform initial sampling ' ...
            'selected but maximum parameter values not defined ' ...
            '-> Set Par_info.max!!']);
    end
end
if strcmp(Par_info.initial,'normal')
    % ERROR -- if normal is used --> mean and covariance of this 
    % distribution need to be defined
    if ~isfield(Par_info,'mu')
        error(['DREAM_Suite ERROR: Normal distribution ' ...
            'selected to sample from but unknown mean ' ...
            '-> Define Par_info.mu!!']);
    end
    if ~isfield(Par_info,'cov')
        error(['DREAM_Suite ERROR: Normal distribution ' ...
            'selected to sample from but unknown covariance ' ...
            '-> Define Par_info.cov!!']);
    end
end
if strcmp(Par_info.initial,'normal')
    % ERROR -- if normal is used --> mean and covariance of this 
    % distribution need to be defined
    if sum(size(Par_info.mu) ~= [1 DREAMPar.d])
        error(char(strcat(['DREAM_Suite ERROR: Mean of normal ' ...
            'distribution (Par_info.mu) should be a row vector with' ...
            ' DREAMPar.d = '],{' '},num2str(DREAMPar.d),{' '},'values')));
    end
    if sum(size(Par_info.cov) ~= [DREAMPar.d DREAMPar.d])
        error(char(strcat(['DREAM_Suite ERROR: Covariance of normal ' ...
            'distribution (''Par_info.cov'') should be a square matrix' ...
            ' of size DREAMPar.d x DREAMPar.d = '],{' '}, ...
            num2str(DREAMPar.d),{' '},'x',{' '},num2str(DREAMPar.d), ...
            {' '},'values')));
    end
end
if strcmp(Par_info.initial,'prior')
    % ERROR -- if explicit prior is used --> marginals need to be defined
    if ~isfield(Par_info,'prior')
        error(['DREAM_Suite ERROR: Prior distribution selected but' ...
            ' unknown field ''prior'' of structure Par_info ' ...
            '-> Define Par_info.prior (see manual)!!']);
    end
end

% % if strcmp(Par_info.initial,'logprior')
% %     % ERROR -- if explicit log prior is used --> marginals need to be defined
% %     if ~isfield(Par_info,'logprior')
% %         error(['DREAM_Suite ERROR: Log prior distribution selected' ...
% %             ' but unknown field ''logprior'' of structure Par_info ' ...
% %             '-> Define Par_info.logprior (see erratum)!!']);
% %     end
% % end

% % if strcmp(Par_info.initial,'prior')
% %     % WARNING -- prior is defined but logprior field also listed
% %     if isfield(Par_info,'logprior')
% %         warning(['DREAM_Suite WARNING: Prior used but logprior ' ...
% %             'defined in field ''logprior'' of structure Par_info ' ...
% %             '-> Will omit and resort to Par_info.prior (see erratum)!!']);
% %     end
% % end

% % if strcmp(Par_info.initial,'logprior')
% %     % WARNING -- logprior is defined but prior field also listed
% %     if isfield(Par_info,'prior')
% %         warning(['DREAM_Suite WARNING: Log prior used but prior' ...
% %             ' defined in field ''prior'' of structure Par_info ' ...
% %             '-> Will omit and resort to Par_info.logprior (see erratum)!!']);
% %     end
% % end
if strcmp(Par_info.initial,'user')
    % ERROR -- if uniform initial sampling is used -> requires 
    % explicit parameter ranges
    if ~isfield(Par_info,'x0')
        error(['DREAM_Suite ERROR: User initial sampling selected ' ...
            'but starting points chain not defined -> Set Par_info.x0!!']);
    else
        if size(Par_info.x0,1) ~= DREAMPar.N
            error(['DREAM_Suite ERROR: Number of rows of matrix' ...
                ' Par_info.x0 does not equal DREAMPar.N']);
        end
        if size(Par_info.x0,2) ~= DREAMPar.d
            error(['DREAM_Suite ERROR: Number of columns of matrix' ...
                ' Par_info.x0 does not equal DREAMPar.d']);
        end
    end
end
% If field boundhandling not defined --> define no use
if ~isfield(Par_info,'boundhandling'), Par_info.boundhandling = 'none'; end
% Check: do we sample in normalized [0-1] space, or not?
if isfield(Par_info,'norm')
    if ~isreal(Par_info.norm)
        error(char(strcat(['DREAM_Suite ERROR: Par_info.norm' ...
            ' should be a scalar -> Define Par_info.norm = 0 or 1\n'])));
    elseif ~ismember(Par_info.norm,[0 1])
        error(char(strcat(['DREAM_Suite ERROR: Par_info.norm' ...
            ' should be zero or one -> Define Par_info.norm = 0 or 1\n'])));
    end
    if Par_info.norm == 1
        if ~isfield(Par_info,'min')
            error(['DREAM_Suite ERROR: Parameter normalization is used ' ...
                'but minimum parameter values not defined -> Set Par_info.min!!']);
        end
        if ~isfield(Par_info,'max')
            error(['DREAM_Suite ERROR: Parameter normalization is used but' ...
                ' maximum parameter values not defined -> Set Par_info.max!!']);
        end
    end    
elseif ~isfield(Par_info,'norm')
    Par_info.norm = 0; 
end
if ~any(strcmp(Par_info.boundhandling,{'fold','bound','reflect',...
        'none','reject'}))
    evalstr = char(['DREAM_Suite WARNING: Unknown boundary' ...
        ' handling method -> Define Par_info.boundhandling = ' ...
        '''fold''/''bound''/''reflect''/''reject''/''none''!!\n']);
    % Now print warning to screen and to file
    fprintf(evalstr); fprintf(fid,evalstr);
    % Continue without boundary handling
    evalstr = char(['DREAM_Suite WARNING: Treating as unbounded ' ...
        'parameter space!!\n']);
    % Now print warning to screen and to file
    fprintf(evalstr); fprintf(fid,evalstr);
end
if any(strcmp(Par_info.boundhandling,{'fold','bound','reflect','reject'}))
    if ~isfield(Par_info,'min')
        error(['DREAM_Suite ERROR: Boundary handling is used ' ...
            'but minimum parameter values not defined -> Set Par_info.min!!']);
    end
    if ~isfield(Par_info,'max')
        error(['DREAM_Suite ERROR: Boundary handling is used but' ...
            ' maximum parameter values not defined -> Set Par_info.max!!']);
    end
end
if any(strcmp(Par_info.initial,{'latin','uniform'})) || ...
        any(strcmp(Par_info.boundhandling,{'fold','bound',...
        'reflect','reject'}))
    % ERROR if Par_info.min or Par_info.max not equal to DREAMPar.d
    if ( numel(Par_info.min) ~= DREAMPar.d )
        error(char(strcat(['DREAM_Suite ERROR: Number of elements ' ...
            'of field ''min'' of structure Par_info should be ' ...
            'equal to '],{' '},num2str(DREAMPar.d),'!!')));
    end
    % ERROR if Par_info.min or Par_info.max not equal to DREAMPar.d
    if ( numel(Par_info.max) ~= DREAMPar.d )
        error(char(strcat(['DREAM_Suite ERROR: Number of elements ' ...
            'of field ''max'' of structure Par_info should be equal' ...
            ' to '],{' '},num2str(DREAMPar.d),'!!')));
    end
end
% Remove dummy variable to make things easier
if strcmp(Par_info.boundhandling,'none')
    Par_info = rmfield(Par_info,'boundhandling');
end
% JAV: Jan. 2022: Check whether user specified names of parameters
if isfield(Par_info,'names')
    if numel(Par_info.names) ~= DREAMPar.d
        error(char(strcat(['DREAM_Suite ERROR: Number of elements' ...
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
            error(['DREAM_Suite ERROR: Length of Meas_info.Sigma ' ...
                'is not equal to that of the observations stored ' ...
                'in Meas_info.Y!!']);
        end
        % Now check this as well
        if any(Meas_info.Sigma < 0)
            fprintf(['DREAM PACKAGE WARNING: One or more entries of ' ...
                'Meas_info.Sigma is negative - we use absolute values\n']);
                 Meas_info.Sigma = abs(Meas_info.Sigma);
        end
        if any(Meas_info.Sigma == 0)
            fprintf(['DREAM PACKAGE WARNING: One or more entries of ' ...
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

if isfield(options,'print')
    if ~ischar(options.print)
        warning(['DREAM_Suite WARNING: Field ''print'' of structure ' ...
            'options should be a string (content equal to ' ...
            '''yes'' or ''no'') ']);
        % Set to no
        options.print = 'yes';
    end
end
if isfield(options,'epsilon')
    if ( options.epsilon <= 0 )
        error(['DREAM_Suite ERROR: Value of ''epsilon'' of ' ...
            'structure options should be larger than zero ' ...
            '-> Set options.epsilon > 0 (default 0.025)']);
    end
end
if isfield(options,'rho')
    if isreal(options.rho)
        error(['DREAM_Suite ERROR: Field ''rho'' of structure ' ...
            'options should be defined as inline function ' ...
            '-> Default: options.rho = inline(''abs(X-Y)'');']);
    end
end
if isfield(options,'burnin')
    if ( options.burnin <= 0 )
        error(['DREAM_Suite ERROR: Value of ''burnin'' of ' ...
            'structure options should be larger than zero ' ...
            '-> Set options.burnin > 0 (default 50 [%])']);
    end
    if ( options.burnin >= 100 )
        error(['DREAM_Suite ERROR: Value of ''burnin'' of ' ...
            'structure options should be smaller than hundred ' ...
            '-> Set options.burnin < 100 (default 50 [%])']);
    end
    warning(['DREAM_Suite WARNING: Burn in percentage defined ' ...
        'in field ''burnin'' of structure options -> Stay ' ...
        'between 50 - 80% - unless you are an expert']);
end
% Check content of each field of structure options
name = fieldnames(options);
% Warning if content not equal to yes or no
for j = 1 : numel(name)
    % Get content of respective field of options
%    F = getfield(options,char(name(j)));
    F = options.(char(name(j)));
    % now check content
    if ~any(strcmp(name(j),{'epsilon','rho','burnin'}))
        if ~any(strcmp(F,{'yes','no'}))
            % ERROR -- content field name(j) of options should be yes or no
            error(char(strcat('DREAM_Suite ERROR: Field', ...
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
            evalstr = char(['DREAM_Suite WARNING: Distributed' ...
                ' Computing toolbox not available -> code resorts ' ...
                'to default: options.parallel = ''no'' \n']);
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
            % And update field parallel of options
            options.parallel = 'no';
        end
    end
end
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
fclose(fid);        % Now close warning_file.txt file

end