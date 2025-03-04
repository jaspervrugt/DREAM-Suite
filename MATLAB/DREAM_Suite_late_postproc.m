function varargout = DREAM_Suite_late_postproc
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
% This function generates tables and figures of DREAM-Suite results       %
% outside main program (DREAM_Suite.m)                                    %
%                                                                         %
% SYNOPSIS: [chain,output,FX,Z] = DREAM_Suite_late_postproc               %
%                                                                         %
% Step 1: load DREAM_Suite.mat                                            %
% Step 2: DREAM_Suite_late_postproc                                       %
%                                                                         %
% © Written by Jasper A. Vrugt, Feb 2007                                  %
% Los Alamos National Laboratory 			        	                  %
% University of California Irvine                                         %
%                                                                         %
% Release Version Jan. 2018     Publication of DREAM manual in EMS        %
% Update          Nov. 2021     Chngd visual outpt & DREAM_Package_end    %
% Update          June 2023     Included new likelihood functions         %
% Update          July 2023     Scoring rules and other additions         %
% Update          Aug. 2024     Major overhaul code & integration MTDREAM %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Load the memory of DREAM Package
load DREAM_Suite.mat method DREAMPar f_handle Par_info ...
    Meas_info Lik_info options MAP_info chain output it iloc

% Run postprocessor (tables and figures)
[chain,output,FX,Z] = DREAM_Suite_postproc(method,DREAMPar,f_handle,...
    Par_info,Meas_info,Lik_info,options,chain,output,it,iloc,MAP_info);
% Close the cluster - if used
DREAM_Suite_end(DREAMPar,options);
% Now define return arguments
varargout(1:4) = {chain,output,FX,Z};

end
