function [Xr,v] = Boundary_handling(X,Par_info)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Function checks whether parameter values are within prior (min/max) bounds         %%
%% If outside of bound, then parameter values are corrected to satisfy bounds         %%
%%                                                                                    %%
%% SYNOPSIS: [Xr,v] = Boundary_handling(X,Par_info)                                   %%
%%  where                                                                             %%
%%   X         [input]  REQUIRED: N x d matrix of candidate points                    %%
%%   Par_info  [input]  REQUIRED: Parameter structure with min/max and bndry treatmnt %%
%%   Xr        [outpt]  Nxd matrix with revised candidate points                      %%
%%   v         [outpt]  Nx1 vector with 0 for in bound and 1 for out of bound         %%
%%                                                                                    %%
%% © Written by Jasper A. Vrugt, Feb 2007                                             %%
%% Los Alamos National Laboratory                                                     %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

Xr = X; [N,d] = size(X); v = zeros(N,1) > 1; % logical array

mn = repmat(Par_info.min,N,1); mx = repmat(Par_info.max,N,1);
%     otherwise
%         mn = repmat(Par_info.minun,N,1); mx = repmat(Par_info.maxun,N,1);
% end

id_l = find(X < mn);            % Smaller than lower bound
id_u = find(X > mx);            % Larger than upper bound

switch Par_info.boundhandling
    case 'reflect'  % reflection
        Xr(id_l) = 2 * mn(id_l) - X(id_l);          % reflect in min
        Xr(id_u) = 2 * mx(id_u) - X(id_u);          % reflect in max
    case 'bound'    % set to bound
        Xr(id_l) = mn(id_l); Xr(id_u) = mx(id_u);
    case 'fold'     % folding
        Xr(id_l) = mx(id_l) - (mn(id_l)-X(id_l));   % Fold lower values
        Xr(id_u) = mn(id_u) + (X(id_u)-mx(id_u));   % Fold upper values
    case 'reject'   % rejection
        o = zeros(N,d); o(id_l) = 1; o(id_u) = 1;   % out of bound
        v = sum(o,2) > 0;                           % Nx1 0/1 vector
        % --> likelihood = -inf (assigned later to v==1)
end
% Reflection/folding: check if all elements are within bounds
% Both methods can go out of bound if violation exceeds mx - mn
% Then, take a random point in space: maybe better to take a random point
% between Xold and bound [FUTURE RELEASE?]
if any(strcmp(Par_info.boundhandling,{'reflect','fold'}))
    id_l = find(Xr < mn);                   % Smaller than lower bound
    id_u = find(Xr > mx);                   % Larger than upper bound
    Xr(id_l) = unifrnd(mn(id_l),mx(id_l));  % Random draw in [mn , mx]
    Xr(id_u) = unifrnd(mn(id_u),mx(id_u));  % Random draw in [mn , mx]
end

end