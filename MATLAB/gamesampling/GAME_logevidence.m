function logZ = GAME_logevidence(t,logq00,logq10,logq01,logq11,n1)
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
% Bridge sampling estimation of Bayesian model evidence. This function    %
% will compute the log normalizing constant of q1, which is the           %
% logarithmic value of the Bayesian model evidence                        %
%                                                                         %
% MAIN IDEA: Take q1 as unnormalized posterior, and q0 as "importance"    %
% distribution, i.e. a normalized density that ideally approximates the   %
% shape of q1                                                             %
%                                                                         %
% SYNOPSIS: logZ = GAME_logevidence(t,logq00,logq10,logq01,logq11,n1)     %
%                                                                         %
% WHERE                                                                   %
%  t           [input] specifies the type of estimator                    %
%   = 0            reciprocal imprtnce smplng [= harm est. if q0 = prior] %
%   = 1            importance sampling                                    %
%   = [0,1]        geometric bridge sampling                              %
%   = otherwise    optimal bridge sampling                                %
%  logq00      [input] n0x1 vector log(q0) n0 samples q0 [not used t=0]   %
%  logq10      [input] n0x1 vector log(q1) n0 samples q0 [not used t=0]   %
%  logq01      [input] n1x1 vector log(q0) n1 samples q1 [not used t=1]   %
%  logq11      [input] n1x1 vector log(q1) n1 samples q1 [not used t=1]   %
%  n1          [input] OPT: # samples of q1                               %
%  logZ        [outpt] logarithmic value of Bayesian model evidence       %
%                                                                         %
% © Written by Gerrit Schoups, July 2016                                  %
% Delft Technical University, The Netherlands                             %
% Modified by Jasper A. Vrugt, Aug. 2016                                  %
% University of California Irvine                                         %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

if nargin < 6, n1 = numel(logq01); end

if (t == 0)
    % qhalf = q0 (RIS)
    logZ = -logmeanexp(logq01-logq11);
elseif (t == 1)
    % qhalf = q1 (IS)
    logZ = logmeanexp(logq10-logq00);
elseif (t > 0 && t < 1)
    % qhalf = geometric bridge = q0^(1-t) * q1^t
    logZ = logmeanexp(t*(logq10-logq00)) - ...
        logmeanexp((1-t)*(logq01-logq11));
else
    % qhalf = optimal bridge = 1 / (Z*s0/q1 + s1/q0)
    n0 = numel(logq00);
    % n1 = numel(logq01);
    logn  = log(n0 + n1);
    logs0 = log(n0) - logn;
    logs1 = log(n1) - logn;
    logl0 = logq10 - logq00;
    logl1 = logq11 - logq01;
    % initial guess
    if (t == 10)
        % use RIS as initial guess
        logZ = -logmeanexp(-logl1);
    elseif (t > 10 && t < 11)
        % use geometric bridge as initial guess
        t = t - 10;
        logZ = logmeanexp(t*logl0) - logmeanexp(-(1-t)*logl1);
    else
        % use IS as initial guess
        logZ = logmeanexp(logl0);
    end
    % iterate
    for i = 1:10
        logs0Z0 = (logs0+logZ)*ones(n0,1);
        logs0Z1 = (logs0+logZ)*ones(numel(logq01),1);   % ones(n1,1)
        logdenom0 = logsumexp([logs0Z0 logs1+logl0],2); % log(s0*Z + s1*l0)
        logdenom1 = logsumexp([logs0Z1 logs1+logl1],2); % log(s0*Z + s1*l1)
        logZ = logmeanexp(logl0-logdenom0) - logmeanexp(-logdenom1);
    end
end

end
