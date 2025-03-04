function [accept,id,ii] = mtMetropolis_rule(log_wp,log_wr,DREAMPar,id_pop)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%  MM   MM TTTTTTTT    DDDDD   RRRRR  EEEEEE   AA   MM   MM               %
%  MM   MM TTTTTTTT    DDDDDD  RRRRRR EEEEEE  AAAA  MM   MM               %
%  MMM MMM TT TT TT    DD   DD RR  RR EE     AA  AA MMM MMM               %
%  MMMMMMM    TT    -- DD   DD RR  RR EEEE   AA  AA MMMMMMM ZZZZZZ SSSSSS %
%  MMM MMM    TT    -- DD   DD RRRRRR EEEE   AAAAAA MMM MMM    ZZZ SSS    %
%  MM   MM    TT       DD   DD RR  RR EE     AAAAAA MM   MM   ZZZ  SSSSSS %
%  MM   MM    TT       DDDDDD  RR  RR EEEEEE AA  AA MM   MM ZZZ       SSS %
%  MM   MM    TT       DDDDD   RR  RR EEEEEE AA  AA MM   MM ZZZZZZ SSSSSS %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
% Multi-try Metropolis rule for acceptance or rejection of proposals      %
%                                                                         %
% SYNOPSIS:                                                               %
%  [accept,id,ii] = mtMetropolis_rule(log_wp,log_wr,DREAMPar,id_pop)      %
% where                                                                   %
%  log_wp      [input] mtxN matrix Σ logpr, loglik & logsnooker proposals %
%  log_wr      [input] mtxN matrix Σ logpr, loglik & logsnooker reference %
%  DREAMPar    [input] Structure with algorithmic variables               %
%  id_pop      [input] Nx1 vector of selected proposals in MT population  %
%  accept      [outpt] Nx1 vector of 1 (accept) or 0 (reject) proposals   %
%  id          [outpt] vector of element # accepted proposals population  %
%  ii          [outpt] vector of element # accepted proposals             %
%                                                                         %
%  © Written by Jasper A. Vrugt, Feb 2007                                 %
%  Los Alamos National Laboratory                                         %
%                                                                         %
% FURTHER CHECKING                                                        %
%  Website:  http://faculty.sites.uci.edu/jasper                          %
%  Papers: http://faculty.sites.uci.edu/jasper/publications/              %
%  Google Scholar: https://scholar.google.com/citations?user=...          %
%                       zkNXecUAAAAJ&hl=nl                                %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

accept = zeros(DREAMPar.N,1);       % Initialize accept
log_wr = max(log_wr,-realmax);      % -Inf or NaN --> -realmax 
log_wr = min(log_wr,realmax);       % Inf --> realmax
% log_wr(isnan(log_wr)) = -realmax;   % Remove nan numbers and replace with 
%                                     % bad log-likelihood value

mx = max(max(log_wp),max(log_wr));  % Maximum of log(w) each chain
wp = max(exp(log_wp - mx),realmin); % Importance weights candidate points
wr = max(exp(log_wr - mx),realmin); % Importance weights reference points
alfa = sum(wp) ./ sum(wr);          % Metropolis acceptance probability
Z = rand(1,DREAMPar.N);             % Draw uniform random numbers

ii = alfa > Z;                      % alfa's greater than Z, accepted
accept(ii,1) = 1;                   % These chains have been accepted
id = id_pop(ii);                    % Relate accept. samples to indx id_pop

end

% Old
% alfa = nan(DREAMPar.N,1);           % Initlze Metropolis accept probability
% % Scale the weights for under or overflow and compute alfa(X --> X_p)
% for zz = 1:DREAMPar.N
%     mx = max(max(log_wp( ...        % Scale log density with maximum to
%         1:DREAMPar.mt,zz)), ...     % avoid problems with zero
%         max(log_wr(1: ...
%         DREAMPar.mt,zz)));
%     w1 = max(exp(log_wp( ...        % Candidate importance weights
%         1:DREAMPar.mt,zz)-mx), ...
%         realmin);
%     w2 = max(exp(log_wr( ...        % Reference importance weights
%         1:DREAMPar.mt,zz)-mx), ...
%         realmin);
%     alfa(zz) = sum(w1)/sum(w2);     % Metropolis acceptance probability
% end
% Z = rand(DREAMPar.N,1);             % Draw uniform random numbers

