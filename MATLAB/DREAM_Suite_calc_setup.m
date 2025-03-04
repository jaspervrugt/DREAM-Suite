function [DREAMPar,f_handle,T] = DREAM_Suite_calc_setup(DREAMPar,...
    fname,options,plugin)
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
% This function initializes computational environment of DREAM-Suite      %
%                                                                         %
% SYNOPSIS: [DREAMPar,f_handle,T] = DREAM_Suite_calc_setup(...            %
%               DREAMPar,fname,options,plugin)                            %
%                                                                         %
% © Written by Jasper A. Vrugt, Feb 2007                                  %
% Los Alamos National Laboratory                                          %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

if isempty(plugin)
    f_handle = eval(['@(x)',char(fname),'(x)']);        % Function handle
else
    f_handle = eval(['@(x)',char(fname),'(x,plugin)']);
end
T = 2;                                                  % First generation
switch options.parallel                                 % Parallel chains?
    case 'no'
        DREAMPar.CPU = 1;                               % 1 CPU (processor)
    case 'yes'
        if verLessThan('matlab','8.2')                  % MATLAB version
            vMATLAB = 'old';
        else
            vMATLAB = 'new';
        end       
        switch vMATLAB                                  % Open cluster, MATLAB version
            case 'old'  % Use matlabpool
                isOpen = matlabpool('size') > 0;        %#ok<DPOOL> Cluster open?
                if ~isOpen, matlabpool open, end        %#ok<DPOOL> If not, new pool created 
                workers = matlabpool('size');           %#ok<DPOOL> # processors
                if ( workers > DREAMPar.N )             % Close the cluster
                    matlabpool close force local        %#ok<DPOOL> 
                    evalstr = strcat(['matlabpool ' ... % Reopen DREAMPar.N cores
                        'open'],{' '}, ...
                        num2str(DREAMPar.N)'); 
                    eval(char(evalstr));
                    workers = DREAMPar.N;               % # workers is DREAMPar.N
                end
            case 'new'  % Use parpool
                pool = gcp;                             % Cluster open?
                if isempty(pool)                        % If not, new pool created
                    pool = parpool('local'); 
                end
                workers = pool.NumWorkers;              % # processors
                if ( workers > DREAMPar.N )
                    delete(gcp)                         % Close the cluster
                    parpool('local',DREAMPar.N);        % Reopen cluster DREAMPar.N cores
                    workers = DREAMPar.N;               % # workers is DREAMPar.N
                end
        end
        DREAMPar.CPU = workers;                         % DREAMPar.CPU is "workers"
        fid = fopen('warning_file.txt','a+');           % Open warning_file.txt
        evalstr = char(strcat(['DREAM_Suite ' ...       % Print to screen
            'PARALLEL: MATLAB pool ' ...                
            'opened with'],{' '},num2str(DREAMPar.CPU), ...
            {' '},'workers for',{' '}, ...
            num2str(DREAMPar.N),{' '},'chains \n'));
        fprintf(evalstr); fprintf(fid,evalstr);         % Print warning to screen & file
        fclose(fid);                                    % Close warning.txt
    
        if strcmp(options.IO,'yes')                     % If input/output writing, 
                                                        % we need directories each worker
            if (ispc || ismac)                          % Copy files to work directories
                copyfile('*.*',[pwd,'/temp_DREAM'])     % copy first to temporary directory
                for ii = 1:DREAMPar.CPU                 % now move to this directory
                    copyfile([pwd,'/temp_DREAM/', ...
                        '*.*'],[pwd,'/',num2str(ii)]);
                end
                rmdir([pwd,'/temp_DREAM'],'s');         % remove temporary directory
            end
            % LINUX/UNIX ENVIRONMENT BUT WORKS FOR PC AS WELL
        end
    
end

end
