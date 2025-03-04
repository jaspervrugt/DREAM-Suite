function [DREAMPar,Par_info,Meas_info,Lik_info,options,f_handle,...
    MAP_info,chain,output,X,FX,S,Table_gamma,CR,pCR,cCR,TdCR,...
    sdX_Xp,loglik,EX,std_mX,id_r,id_c,iloc,it,g,LLX,T_start] = ...
    DREAM_Suite_restart(fname)                                   %#ok
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
% This function reloads variables for a restart run                       %
%                                                                         %
% SYNOPSIS: [DREAMPar,Par_info,Meas_info,Lik_info,options,f_handle,...    %
%               MAP_info,chain,output,X,FX,S,Table_gamma,CR,pCR,cCR,...   %
%               TdCR,sdX_Xp,loglik,EX,std_mX,iloc,it,g,LLX,...            %
%               T_start] = DREAM_Suite_restart(fname)                     %
%                                                                         %
% © Written by Jasper A. Vrugt, Feb 2007                                  %
% Los Alamos National Laboratory                                          %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Try ... catch statement in case restart .mat file does not exist
try
    % String with data to read from fname
    evalstr = strcat('load',{' '},fname,{' '},['DREAMPar Par_info ' ...
        ['Meas_info Lik_info options MAP_info output X FX S chain ' ...
        'Table_gamma CR pCR cCR TdCR sdX_Xp loglik EX std_mX iloc ' ...
        'id_r id_c it g LLX t Func_name plugin']]);
    % Load the restart data
    eval(char(evalstr));
catch
    % Create error warning to screen
    evalstr = char(strcat(['DREAM_Suite ERROR: Cannot restart --> ' ...
        'File '],{' '},fname,{' '},['does not exist. Next run, to ' ...
        'avoid this problem, set field ''save'' of structure options' ...
        ' to ''yes'' ']));
    % Now print error
    error(evalstr);
end

% Sample randomly the crossover values [ 1 : DREAMPar.nCR ]
CR = reshape(randsample(1:DREAMPar.nCR,DREAMPar.N * DREAMPar.steps, ...
    'true',pCR),DREAMPar.N,DREAMPar.steps); %#ok
% Open warning file
fid = fopen('warning_file.txt','a+');

% Check whether previous run was aborted early or not
switch ( t < DREAMPar.T )
    case 1      % Finish previous budget
        evalstr = char(strcat(['DREAM_Suite RESTART: Starting with' ...
            ' t = '],{' '},num2str(t),{' '},['but still using ' ...
            'old budget of'],{' '},'DREAMPar.T = ',{' '}, ...
            num2str(DREAMPar.T),'\n'));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);

    otherwise   % Assign new budget and run

        fprintf(fid,['-------------- DREAM_Suite warning file -------' ...
            '------ \n']);
        fname = 'T.txt';    % Now load file
        % Check whether file T.txt exists
        if exist(fname,'file') == 2 % Yes
            evalstr = char(strcat(['DREAM_Suite RESTART: Located ' ...
                'file ''T.txt'' in respective example directory - now ' ...
                'checking its content\n']));
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
            % Now check content of T.txt to avoid crash of code
            T_new = checkfile_T(fname);
            % If all are satisfied then
            evalstr = char(strcat(['DREAM_Suite RESTART: User has ' ...
                'requested/listed'],{' '},num2str(T_new),{' '}, ...
                'additional generations in file ''T.txt''\n'));
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
            % Now print to file/screen the final message
            evalstr = char(strcat(['DREAM_Suite RESTART: Initial ' ...
                't = '],{' '},num2str(t),{' '},'and completing',{' '}, ...
                num2str(T_new),{' '},'additional generations so that', ...
                {' '},'DREAMPar.T = ',{' '}, ...
                num2str(DREAMPar.T  + T_new),'\n'));
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
        
        else % File does not exist
        
            evalstr = char(strcat(['DREAM_Suite RESTART: Did not' ...
                'locate file ''T.txt'' in respective example directory' ...
                'so per default we double the "previous" budget of' ...
                'generations\n']));
            % Now print warning to screen and to file   
            fprintf(evalstr); fprintf(fid,evalstr);
            % T_new set equal to DREAMPar.T (= double budget)
            T_new = DREAMPar.T;        
            % Now print to file/screen the final message
            evalstr = char(strcat('DREAM_Suite RESTART: Initial t = ', ...
                {' '},num2str(t),{' '},'and completing',{' '}, ...
                num2str(T_new),{' '},...
                'additional generations so that',{' '},'DREAMPar.T = ', ...
                {' '},num2str(DREAMPar.T  + T_new),'\n'));
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);

        end
        [T,d2,N] = size(chain);         %#ok                     
        % Add nan to chains
        chain(T+1:T+T_new,1:d2,1:N) = nan(T_new,d2,N);
        % Update DREAMPar.T with T_new
        DREAMPar.T = DREAMPar.T + T_new;   
end

% close warning_file.txt
fclose(fid);
% Define starting value of T
T_start = t + 1;
% Now initialize computational framework
[DREAMPar,f_handle] = DREAM_Suite_calc_setup(DREAMPar,Func_name,...
    options,plugin);

end
