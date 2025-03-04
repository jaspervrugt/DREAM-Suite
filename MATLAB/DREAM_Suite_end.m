function DREAM_Suite_end(DREAMPar,options)
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
% This function closes pool of workers and removes IO files               %
%                                                                         %
% SYNOPSIS: DREAM_Suite_end(DREAMPar,options)                             %
%                                                                         %
% © Written by Jasper A. Vrugt, Feb 2007                                  %
% Los Alamos National Laboratory                                          %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

if DREAMPar.CPU > 1
    if verLessThan('matlab','8.2')      
        parpool('close');               % Close MATLAB pool
    else
        delete(gcp('nocreate'));        % Close in latest version
    end
    if strcmp(options.IO,'yes')         % If IO writing, remove directories
        try                             % Remove each worker directory
            for ii = 1:DREAMPar.CPU
                rmdir([pwd,'/', ...
                    num2str(ii)],'s');
            end
        catch                           % Otherwise, return to header
            return 
        end
    end
end

end
