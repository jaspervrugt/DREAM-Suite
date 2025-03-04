function [i,I,flag] = Haverkamp_I_der(eta,plugin,rtol,kmax,itol)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% This function solves for the I(t) relationship of Haverkamp using Newton's method  %%
%%  SYNOPSIS: I = Haverkamp_I_patch(eta,t)                                            %%
%%            I = Haverkamp_I_patch(eta,t,rtol)                                       %%
%%            I = Haverkamp_I_patch(eta,t,rtol,kmax)                                  %%
%%            I = Haverkamp_I_patch(eta,t,rtol,kmax,itol)                             %%
%% where                                                                              %%
%%  eta  [input]       4x1 vector with S [cm/h^0.5], Ks [cm/h], beta [-], Ki [cm/h]   %%
%%  t    [input]       nx1 vector with time, t, in hours (h)                          %%
%%  rtol [opt. input]  tolerance on function value at root        (default: 1e-12)    %%
%%  kmax [opt. input]  maximum number of Newton iterations        (default: 20)       %%
%%  itol [opt. input]  tolerance on constant rate assumption      (default: 1e-10)    %%
%%  I    [output]      cumulative infiltration, I (cm), as function of time, t (h)    %%
%%  i    [output]      infiltration rate, i (cm/h), as function of time, t (h)        %%
%%  flag [output]      exit flag: [1] exact [2] approximate                           %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
if nargin < 5, itol = 1e-10; end
if nargin < 4, kmax = 20; end
if nargin < 3, rtol = 1e-12; end
%% Initialization
S = eta(1); Ks = eta(2); B = eta(3); Ki = 0;        % Unpack parameter values
dK = Ks - Ki; xi = dK/S^2; t = plugin.t(:);         % Compute dK, xi and t column vector
n = numel(t); [I,i] = deal(nan(n,1)); flag = 1;     % Initialize cum. inf., rate/flag
if t(1) == 0, I(1)=0; i(1)=inf; j=2; else, j=1; end % Set I(1) zero if t(1) is zero
dt = [ t(1) ; diff(t) ];                            % Compute delta time
r = @(I,t) I-1/(2*xi)*log(exp(2*xi*B*(I-Ki*t))/B + ...
    (B-1)/B) - (dK*(1-B)+Ki)*t;                     % Residual function
dr = @(I,t) (1-(B*exp(2*xi*B*(I-Ki*t)))./ ...
    (exp(2*xi*B*(I-Ki*t))+B-1));                    % Derivative function
%% Dynamic part
I_up = S*sqrt(t) + t*Ks; I_low = max(I_up(1)/2,.1); % Upper/lower limits root bracket
gamma = 2*xi*B*(I_up-Ki*t);                         % Compute gamma
for m = j:n                                         % For loop
    if gamma(m) < 709.783                               % Check overflow condition
        y(1) = I_low; k = 1;                            % Root guess and iter.
        while (abs(r(y(k),t(m)))>rtol) && (k<kmax)          % While loop
            y(k+1) = y(k) - r(y(k),t(m))/dr(y(k),t(m));         % Next iterate
            k = k + 1;                                          % Increment iteration
        end                                                 % End while loop
        [I(m),I_low] = deal(y(k));                          % I(m) is equal to root
        A = exp(2*B*xi*(I(m)-Ki*t(m)));                     % Compute alpha
        i(m) = (((1-B)*Ks+B*Ki)*(A+B-1)-A*B*Ki)/(A+B-A*B-1);% Compute inf. rate
    elseif (m > 1)                                      % Patch for overflow
        if abs(i(m-1)-Ks) < itol                            % Constant inf rate?
            i(m) = Ks; I(m) = I(m-1) + dt(m)*i(m);              % Constant rate
        else                                                % Else
            flag = 2; break;                                    % flag = 2; break loop
        end                                                 % End if statement    
    else                                                % Else
        flag = 2; break;                                    % flag = 2; break loop
    end                                                 % End if statement
end                                                 % End for loop    
i = diff(I); 
% % function [i,I,flag] = Haverkamp_I_der(eta,plugin,rtol,kmax)
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% % %% This function solves for the I(t) relationship of Haverkamp using Newton's method  %%
% % %%  SYNOPSIS: I = Haverkamp_I(eta,t)                                                  %%
% % %%            I = Haverkamp_I(eta,t,rtol)                                             %%
% % %%            I = Haverkamp_I(eta,t,rtol,kmax)                                        %%
% % %% where                                                                              %%
% % %%  eta  [input]       4x1 vector with S [cm/h^0.5], Ks [cm/h], beta [-], Ki [cm/h]   %%
% % %%  t    [input]       nx1 vector with time, t, in hours (h)                          %%
% % %%  rtol [opt. input]  tolerance on function value at root        (default: 1e-12)    %%
% % %%  kmax [opt. input]  maximum number of Newton iterations        (default: 20)       %%
% % %%  I    [output]      cumulative infiltration, I (cm), as function of time, t (h)    %%
% % %%  i    [output]      infiltration rate, i (cm/h), as function of time, t (h)        %%
% % %%  flag [output]      exit flag: [1] exact [2] approximate                           %%
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% % if nargin < 4, kmax = 20; end 
% % if nargin < 3, rtol = 1e-12; end
% % if numel(eta) == 3, eta(4) = 0; end
% % %% Initialization
% % S = eta(1); Ks = eta(2); B = eta(3); Ki = eta(4);   % Unpack parameter values
% % dK = Ks-Ki; xi = dK/S^2; t = plugin.t(:);           % Compute dK and xi, # elements t
% % n = numel(t); I = nan(n,1); flag = 1;               % Initial cum. inf.
% % if t(1) == 0, I(1) = 0; j = 2; else, j = 1; end
% % r = @(I,t) I-1/(2*xi)*log(exp(2*xi*B*(I-Ki*t))/B + ...
% %     (B-1)/B) - (dK*(1-B)+Ki)*t;                     % Residual function
% % dr = @(I,t) (1-(B*exp(2*xi*B*(I-Ki*t)))./ ...
% %     (exp(2*xi*B*(I-Ki*t))+B-1));                    % Derivative function
% % %% Dynamic part
% % I_up = S*sqrt(t) + t*Ks; I_low = max(I_up(1)/2,.1); % Upper/lower limits root bracket
% % gamma = 2*xi*B*(I_up-Ki*t);                         % Compute gamma 
% % for m = j:n                                         % For loop
% %      if gamma(m) < 709.783                              % Check overflow condition   
% %         y(1) = I_low; k = 1;                                % Yes: root guess and iter.    
% %         while (abs(r(y(k),t(m)))>rtol) && (k<kmax)          % While loop
% %             y(k+1) = y(k) - r(y(k),t(m))/dr(y(k),t(m));         % Next iterate
% %             k = k + 1;                                          % Increment iteration 
% %         end                                                 % End while loop    
% %         [I(m),I_low] = deal(y(k));                          % I(m) is equal to root
% %     else                                                % Otherwise 
% %         flag = 2; break                                     % No: flag=2, break loop 
% %     end                                                 % End if statement
% % end                                                 % End for loop
% % A = exp(2*B*xi*(I-Ki*t));                           % Compute alpha
% % i = (((1-B)*Ks+B*Ki)*(A+B-1)-A*B*Ki)./(A+B-A*B-1);  % Compute infiltration rate