function X_dis = Discrete_space(X,Par_info)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% This function transforms continuous space of x to discrete values                  %%
%%                                                                                    %%
%% SYNOPSIS: X_d = Discrete_space(X,Par_info)                                         %%
%%  where                                                                             %%
%%   X         [input]  REQUIRED: N x d matrix of candidate points                    %%
%%   Par_info  [input]  REQUIRED: Parameter structure                                 %%
%%   X_dis     [outpt]  N x d matrix of discrete parameter values                     %%
%%                                                                                    %%
%% (c) Written by Jasper A. Vrugt, Feb 2007                                           %%
%% Los Alamos National Laboratory 			        	                              %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

method = 2;                                     % latest MATLAB release
N = size(X,1);                                  % # candidate vectors? 

% Step 1: Transform continuous x to integer between 0 and # steps
% Step 2: Back transform to discrete space
switch method
    case 1 % Proper method - for all releases
        X_min = repmat(Par_info.min,N,1);              
        X_int = round(repmat(Par_info.steps,N,1) .* ... 
            ((X-X_min) ./ repmat(Par_info.max - Par_info.min,N,1))); 
        X_dis = X_min + X_int .* repmat(Par_info.step_size,N,1); 
    case 2 % New MATLAB - later releases
        X_int = round(Par_info.steps .* ((X - Par_info.min) ./ ...
            (Par_info.max - Par_info.min)));
        X_dis = Par_info.min + X_int .* Par_info.step_size; 
end

end