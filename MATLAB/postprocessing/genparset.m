function P = genparset(chain)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% This function generates a matrix P from sampled chain trajectories DREAM Package   %%
%%                                                                                    %%
%% SYNOPSIS: P = genparset(chain)                                                     %%
%%                                                                                    %%
%% Example: If N = 3 (DREAM_ZS) then sample number (= row) 1, 2 and 3 of P equal the  %%
%%          initial chain states of chains 1,...,N. Then, P(4:6,:) correspond to the  %%
%%          second chain samples of 1,...,N and sample # (row) 18 of P equals the 6th %%
%%          sample of chain 3.                                                        %%
%%          If thinning (DREAMPar.thinning > 1) then, this applies to P as well.      %%
%%          For example, DREAMPar.thinning = 2 (only every 2nd sample is stored)      %%
%%          then P(4:6,:) would store the 3rd sample of the chains and P(18,:) is the %%
%%          12th sample of chain 3.                                                   %%
%%                                                                                    %%
%% (c) Written by Jasper A. Vrugt, Feb 2007                                           %%
%% Los Alamos National Laboratory                                                     %%
%% University of California Irvine                                                    %%
%%                                                                                    %%
%% Update          July 2024     Initialization of P matrix and switch command        %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

[T,d,N] = size(chain);  % Determine # samples, # pars - 2, # chains
switch T
    case 0
        P = [];             % ParSet is empty
    otherwise
        id = (1:T)';        % ID for each chain sample
        P = nan(N*T,d+1);   % Initialize return argument P: N*T x d matrix
        for z = 1:N         % Copy each chain to P and sample IDs
            P((z-1)*T+1:z*T,1:d+1) = [chain(1:T,1:d,z) , id];
        end
        P = sortrows(P,d+1); P = P(1:N*T,1:d);  % Sort on ID and remove ID
end

end