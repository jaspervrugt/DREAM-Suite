function [dtdI,flag] = Haverkamp_t_der(eta,plugin,itol)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% This function solves for the t(I) relationship of Haverkamp                        %%
%%  SYNOPSIS: t = Haverkamp_t_patch(eta,I)                                            %%
%%  SYNOPSIS: t = Haverkamp_t_patch(eta,I,itol)                                       %%
%% where                                                                              %%
%%  eta  [input]       3x1 vector with S [cm/h^0.5], Ks [cm/h], and beta [-]          %%
%%  I    [input]       nx1 vector with cumulative infiltration, I, in cm              %%
%%  itol [opt. input]  tolerance on constant rate assumption      (default: 1e-10)    %%
%%  t    [output]      time, t, in hours corresponding to I (cm)                      %%
%%  flag [output]      exit flag: [1] exact [2] approximate                           %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
if nargin < 3, itol = 1e-10; end
%% Initialization part
S = eta(1); Ks = eta(2); B = eta(3); Ki = 0;        % Unpack parameter values
dK = Ks-Ki; xi = dK/S^2; I = plugin.I(:);           % Compute dK, xi, I is column vector
n = numel(I); t = nan(n,1); flag = 1;               % Initialize time
%% Dynamic part: Evaluate time expression
gamma = 2*B*xi*I;                                   % Compute gamma
if I(1) == 0; t(1)=0; i=inf; j=2; else, j=1; end    % t(1) equals zero if I(1) is zero
for m = j:n                                         % For loop
    if gamma(m) < 709.783                               % Numerical overflow?
        t(m) = 1./(dK*(B-1)*xi).*(1/2* ...                  % No: compute inf. time
            log(exp(gamma(m))/B + (B-1)/B) - xi*I(m));      
        A = exp(2*B*xi*(I(m)-Ki*t(m)));                     % Compute alpha
        i = (1-B)*Ks*(A+B-1)/(A+B-A*B-1);                   % Compute inf. rate        
    elseif (m > 1)                                      % Patch for overflow    
        if abs(i-Ks) < itol                                 % Constant inf rate?
            t(m) = t(m-1) + (I(m)-I(m-1))/Ks; flag = 2;         % Constant rate
        else                                                % Else
            flag = 2; break                                     % flag = 2; break loop
        end                                                 % End if statement    
    else                                                % Else
        flag = 2; break                                     % flag = 2; break loop
    end                                                 % End if statement
end                                                 % End for loop
%dtdI = (exp((2*B*I*Ks)/S^2)) ./ ((B - 1)*(exp((2*B*I(m)*Ks)/S^2)/B + (B - 1)/B)) - 1/(Ks*(B - 1));
%dtdI = diff(t);
% A = exp(2*B*xi*I(m));


% % function [dtdI,flag] = Haverkamp_t_der(eta,plugin)
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% % %% This function solves for the t(I) relationship of Haverkamp                        %%
% % %%  SYNOPSIS: t = Haverkamp_t(eta,I)                                                  %%
% % %% where                                                                              %%
% % %%  eta  [input]  3x1 vector with S [cm/h^0.5], Ks [cm/h], and beta [-]               %%
% % %%  I    [input]  nx1 vector with cumulative infiltration, I, in cm                   %%
% % %%  dtdi [output] derivative of time with respect to I (cm)                           %%
% % %%  flag [output] exit flag: [1] finite [2] infinite time values                      %%
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% % 
% % %% Initialization part
% % S = eta(1); Ks = eta(2); B = eta(3); Ki = 0;        % Unpack parameter values
% % dK = Ks-Ki; xi = dK/S^2; I = plugin.I(:);           % Compute dK and xi
% % n = numel(I); [t,dtdI] = deal(nan(n,1)); flag = 1;  % Initialize time
% % %% Dynamic part: Evaluate time expression
% % gamma = 2*B*xi*I;                                   % Compute gamma
% % if I(1) == 0; t(1) = 0; j = 2; else, j = 1; end     % t(1) equals zero if I(1) is zero
% % for m = j:n                                         % For loop
% %     if gamma(m) < 709.783                               % Provision numerical overflow
% %         t(m) = 1./(dK*(B-1)*xi).*(1/2* ...
% %             log(exp(gamma(m))/B + (B-1)/B) - xi*I(m));      % Compute inf. time
% %         dtdI(m) = (exp((2*B*I(m)*Ks)/S^2))/((B - 1)*(exp((2*B*I(m)*Ks)/S^2)/B + (B - 1)/B)) - 1/(Ks*(B - 1));
% %         % A = exp(2*B*xi*I(m));
% %         % dtdI(m) = (A*B*Ks-A-B+1)./((B-1)*Ks*(A+B-1));
% % % %         A/((B - 1)*(A/B + (B - 1)/B)) - 1/(Ks*(B - 1));
% % % %         A/((B - 1)*( ( A + B - 1)/B)) - 1/(Ks*(B - 1));
% % % %        (A*B*Ks-A-B+1)./((B-1)*Ks*(A+B-1)); %(A-1)*B/((B - 1) * ( A + B + B*Ks - 1 ));
% % %          A = exp(2*B*xi*(I-Ki*t));                           % Compute alpha
% % %          i = (((1-B)*Ks+B*Ki)*(A+B-1)-A*B*Ki)./(A+B-A*B-1);  % Compute infiltration rate
% %     else
% %         flag = 2; break                                 % Break for loop; flag = 2
% %     end                                                 % End if statement
% % end                                                 % End for loop