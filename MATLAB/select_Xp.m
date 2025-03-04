function [log_w,id_pop] = select_Xp(DREAMPar,log_PRL_Xp,log_sn_Xp)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%  MM   MM TTTTTTTT    DDDDD   RRRRR   EEEEEE  AA   MM   MM               %
%  MM   MM TTTTTTTT    DDDDDD  RRRRRR  EEEEEE AAAA  MM   MM               %
%  MMM MMM TT TT TT    DD   DD RR  RR  EE    AA  AA MMM MMM               %
%  MMMMMMM    TT    -- DD   DD RR  RR  EEEE  AA  AA MMMMMMM ZZZZZZ SSSSSS %
%  MMM MMM    TT    -- DD   DD RRRRRR  EEEE  AAAAAA MMM MMM    ZZZ SSS    %
%  MM   MM    TT       DD   DD RR  RR  EE    AAAAAA MM   MM   ZZZ  SSSSSS %
%  MM   MM    TT       DDDDDD  RR  RR EEEEEE AA  AA MM   MM ZZZ       SSS %
%  MM   MM    TT       DDDDD   RR  RR EEEEEE AA  AA MM   MM ZZZZZZ SSSSSS %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
% Computes the importance weights of the multiple try candidate points of %
% each chain. This weights are normalized and used to sample a preferred  %
% candidate point for each chain.                                         %
%                                                                         %
% SYNOPSIS:                                                               %
%  [log_w,id_pop] = select_Xp(DREAMPar,log_PRL_Xp,log_sn_Xp)              %
% where                                                                   %
%  DREAMPar    [input] Structure with algorithmic variables               %
%  log_PRL_Xp  [input] mtx1 vector of log-densities (logpr + loglik)      %
%  log_sn_Xp   [input] mtx1 vector of log snooker transition densities    %
%  log_w       [input] mtxN matrix Σ logpr, loglik & logsnooker proposals %
%  id_pop      [input] Nx1 vector indices preferred candidate points      %
%                                                                         %
%  © Written by Jasper A. Vrugt, Feb 2012                                 %
%  Los Alamos National Laboratory                                         %
%                                                                         %
% FURTHER CHECKING                                                        %
%  Website:  http://faculty.sites.uci.edu/jasper                          %
%  Papers: http://faculty.sites.uci.edu/jasper/publications/              %
%  Google Scholar: https://scholar.google.com/citations?user=...          %
%                       zkNXecUAAAAJ&hl=nl                                %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Check: Eq. 10 of https://arxiv.org/pdf/1801.09065 
log_w = sum(log_PRL_Xp,2) - log_sn_Xp;  % Sum log_prior + loglik - log g(x)
log_w = max(log_w,-realmax);            % -Inf or NaN --> -realmax 
log_w = min(log_w,realmax);             % Inf --> realmax
%log_w(isnan(log_w)) = -realmax;         % Replace nan with minimum value 
                                        % (before -inf is reached)
log_w = reshape(log_w,DREAMPar.mt, ...  % Reshape log_w 
    DREAMPar.N);                        % Separate column for each chain

log_ws = bsxfun(@minus,log_w, ...       % Rescale loglik (numrcl underflow)    
    max(log_w));   
ws = exp(log_ws);                       % Compute rescaled weights 
w = bsxfun(@rdivide,ws,sum(ws));        % Normalize weights each chain
id_pop = nan(DREAMPar.N,1);             % Initialize column vector id_sel

% Now select the best among the DREAMPar.mt proposals in each chain
for i = 1:DREAMPar.N
    % Now select the "best" proposal from the normalized weights
    id = randsample(1:DREAMPar.mt, ...
        1,'true',w(1:DREAMPar.mt,i));
    % Now determine which element of y_mt this is for each chain
    id_pop(i) = id + (i-1) * DREAMPar.mt; 
end

end
