function R = LH_sampling(mn,mx,N)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Latin hypercube sampling of initial chain states DREAM Package                     %%
%%                                                                                    %%
%% SYNOPSIS: R = LH_sampling(mn,mx,N)                                                 %%
%%  where                                                                             %%
%%   mn        [input] REQUIRED: 1 x d vector of lower bound values                   %%
%%   mx        [input] REQUIRED: 1 x d vector of upper bound values                   %%
%%   N         [input] REQUIRED: # of Latin hypercube samples                         %%
%%   r         [outpt] Nxd matrix of Latin hypercube samples                          %%
%%                                                                                    %%
%% Implementation based on the following paper                                        %%
%%  Minasny, B., and A.B. McBratney (2006), A conditioned Latin hypercube method for  %%
%%      sampling in the presence of ancillary information, Computers & Geosciences,   %%
%%      Vol. 32 (9), pp. 1378-138                                                     %%
%%                                                                                    %%
%% © Written by Jasper A. Vrugt, Feb 2007                                             %%
%% Los Alamos National Laboratory                                                     %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

d = size(mn,2);                 % # parameters 
rng = mx - mn;                  % 1xd vector with parameter ranges
y = rand(N,d);                  % Nxd matrix with uniform random labels
[~,id] = sort(rand(N,d));       % Draw at random 1:N without replacement
M = (id - y)/N;                 % Multiplier matrix (y introduces randomness)
R = bsxfun(@plus, ...
    bsxfun(@times,M,rng),mn);   % Nxd matrix of stratified LH samples
%r = mn + M .* rng;

end

% % % Draw LH samples
% % for j = 1:d
% %     P = (randperm(N)' - y(:,j))/N;
% %     r(:,j) = mn(j) + P * rng(j);    % Sample jth parameter
% % end