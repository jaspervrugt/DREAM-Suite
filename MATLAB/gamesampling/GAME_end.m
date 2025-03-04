function GAME_end(DREAMPar,options)
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
% This function close the MATLAB pool (if CPU > 1) and/or removes files   %
%                                                                         %
% SYNOPSIS: GAME_end(DREAMPar,options)                                    %
%                                                                         %
% © Written by Jasper A. Vrugt, Jan 2015                                  %
% University of California Irvine                                         %
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

