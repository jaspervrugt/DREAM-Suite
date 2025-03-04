function J = jac_Haverkamp_anal(x,plugin,app)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Analytic estimate of the nxd Jacobian matrix of the Haverkamp         %%
%% infiltration equation in time and infiltration form                   %%
%%                                                                       %%
%%  SYNOPSIS: J = jac_Haverkamp_anal(x,plugin,app)                       %%
%% where                                                                 %%
%%  x           [input] 1xd vector of Haverkamp parameters               %%
%%                  3-par form, S [cm/h^.5], Ks [cm/h], ß [-]            %%
%%                  4-par form, S [cm/h^.5], Ks [cm/h], ß [-], Ki [cm/h] %%
%%  plugin      [input] structure time, cum. inf. experimental info      %%
%%  app         [input] Residuals to minimimze                           %%
%%   = 1            Minimize cumulative infiltration residuals           %%
%%   = 2            Minimize time residuals [valid if Ki = 0 -> d = 3]   %%
%%  J           [outpt] nxd Jacobian matrix of time/inf form Haverkamp   %%
%%                                                                       %%
%%  LITERATURE                                                           %%
%%  Vrugt, J.A., J.W. Hopmans, Y. Gao, M. Rahmati, J. Vanderborght, and  %%
%%      H. Vereecken (2023), The time validity of Philip's two-term      %%
%%      infiltration equation: An elusive theoretical quantity? Vadose   %%
%%      Zone Journal, e20309, pp. 1-25,                                  %%
%%          https://doi.org/10.1002/vzj2.20309                           %%
%%  Vrugt, J.A. and y. Gao (2022), On the three-parameter infiltration   %%
%%      equation of Parlange et al. (1982): Numerical solution,          %%
%%      experimental design, and parameter estimation, Vadose Zone       %%
%%      Journal, 21:e20167, pp. 1-25,                                    %%
%%          https://doi.org/10.1002/vzj2.20167                           %%
%%                                                                       %%
%%  © Written by Jasper A. Vrugt, Jan. 2019                              %%
%%  University of California Irvine                                      %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Unpack parameter values
S = x(1); Ks = x(2); 
% Check wkether we estimate Ki or not
d = numel(x);
% Check input arguments
if d == 4 && app == 2
    error(['jac_Haverkamp_anal:Time form of Haverkamp implies Ki = 0, ' ...
        'and, consequently, for d = 4, app ≠ 2. For app = 2 to work ' ...
        '(time minimization), one must use the 3-parameter variant of ' ...
        'Haverkamp''s infiltration equation, then Ki = 0']);
end
% Assign parameters
switch d
    case 2 % ß fixed, Ki is zero
        b = plugin.b; Ki = 0;
    case 3 % ß estimated, Ki is zero
        b = x(3); Ki = 0;
    case 4 % ß and Ki estimated
        b = x(3); Ki = x(4);
end
% Infiltration or time form of Haverkamp's equation?
switch app
    case 1 % Infiltration form
        % We know t - and simulate cum. inf., I (cm)
        I = Haverkamp_I(x,plugin); t = plugin.t; t = t(:);
        % Vrugt and Gao, VZJ, 2022
        % Define α: nx1 vector as it depends on time
        a = exp(2*b*(Ks-Ki)*(I-Ki*t)/S^2);
        % Define common terms numerator/denominator, nx1 vectors
        c1 = a + b - 1; c2 = a + b - a*b - 1;
        % ∂I/∂S, Equation 19a
        dIdS = (2*(I-Ki*t).*c2 + 2*t*(b-1)*(Ks-Ki) .* c1) ./ ...
            + (S * c2);
        % Compute ∂I/∂Ks, Equation 19b 
        dIdKs = (a*b.*(I-Ki*t) - (I - Ki*t + 2*t*(b-1)*(Ks-Ki)).*c1)./ ...
            ((Ks-Ki)*c2);   % [mistake VZJ, [I + Kit ... -> [I - Kit ...
        % Compute ∂I/∂ß, Equation 19c
        dIdb = ( S^2*(1-a)/(2*b) + a*(Ks-Ki).*(I-Ki*t) ...
            - t .* c1 * (Ks-Ki)^2) ./ ((Ks-Ki)*c2);
        switch d
            case 2
                J = [dIdS, dIdKs]; 
            case 3
                J = [dIdS, dIdKs, dIdb];
            case 4
                % Compute ∂I/∂Ki, Equation 19d
                dIdKi = ( (I - Ki*t + (2*b - 1)*t*(Ks-Ki)) .* c1 - ...
                    a*b.*t*(Ks-Ki) - a*b.*(I - Ki*t) ) ./ ((Ks-Ki)*c2);
                J = [dIdS, dIdKs, dIdb, dIdKi];
        end
        if t(1) == 0 
            J(1,1:d) = zeros(1,d);
        end
        % CORRECT
        % drds = (2*(Ki-Ks)*(I-Ki*t) - 2*t*(B-1)*(Ki-Ks)^2)/S^3 - ...
        %     2*B/S^3*(Ki-Ks)*(I-Ki*t).*exp(-(2*B*(Ki-Ks)*(I-Ki*t))/ ...
        %     S^2) ./ ( exp(-(2*B*(Ki-Ks)*(I-Ki*t))/S^2) + B - 1);
        % CORRECT
        % drdks =  ( (I-Ki*t)-2*t*(B-1)*(Ki-Ks) )/S^2 - ...
        %     B/S^2 * (I-Ki*t) .* exp(-(2*B*(Ki-Ks)*(I-Ki*t))/S^2)./ ...
        %     (exp(-(2*B*(Ki-Ks)*(I-Ki*t))/S^2) + B - 1);
        % CORRECT
        % drdb = -1/(2*B) * ( 1 - exp(-(2*B*(Ki-Ks)*(I-Ki*t))/S^2) ) ./ ...
        %     ( exp(-(2*B*(Ki-Ks)*(I-Ki*t))/S^2) + B - 1) + (Ki-Ks) * ...
        %     (I-Ki*t)/S^2 .* exp(-(2*B*(Ki-Ks)*(I-Ki*t))/S^2) ./ ...
        %     ( exp(-(2*B*(Ki-Ks)*(I-Ki*t))/S^2) + B - 1) + ...
        %     t*(Ki-Ks)^2/S^2;
        % CORRECT
        % drdki = (t*(Ki-Ks) - (I-Ki*t))/S^2 + (2*B*(I-Ki*t) - ...
        %     2*B*t*(Ki-Ks)).*exp(-(2*B*(Ki-Ks)*(I-Ki*t))/S^2) ./ ...
        %     ( 2*S^2 * ( exp(-(2*B*(Ki-Ks)*(I-Ki*t))/S^2) + B - 1)) + ...
        %     2*t*(B-1)*(Ki-Ks)/S^2;
        % CORRECT
        % drdi =  B/S^2 * (Ki-Ks)*exp(-(2*B*(Ki-Ks)*(I-Ki*t))/S^2) ./ ...
        %     (exp(-(2*B*(Ki-Ks)*(I-Ki*t))/S^2) + B - 1) - (Ki - Ks)/S^2;   
        % Compute Jacobian
        % J = [ -drds./drdi , -drdks./drdi , -drdb./drdi ];
        % if d == 4
        %     J(:,4) = -drdki./drdi;
        % end
        % Check if first entry of t/I equals zero 
        % if t(1) == 0 
        %     J(1,1:d) = zeros(1,d);
        % end

        % Enter limit of B -> 0
        % limB0 = ( t*(Ki - Ks)^2.*((2*(Ki - Ks)*(I - Ki*t))/S^2 - 1) + ...
        %     ((Ki - Ks)^2*(I - Ki*t).^2)/S^2 ) ./ ( (2*(Ki - Ks)^2 * ...
        %     (I - Ki*t))/S^2 )
        % limB0 = ( -t .* [2*(Ks - Ki)*(I - Ki*t) + S^2] + ...
        %     (I - Ki*t).^2 ) ./ ( (2*(I - Ki*t)) )
        % limB0 = -t*(Ks-Ki) -1/2*t*S^2./(I - Ki*t) + 1/2*(I-Ki*t);
        % limB2 = ( -(S^2*(exp(-(4*(Ki - Ks)*(I - Ki*t))/S^2)/2 - ...
        %     1/2))/2 - exp(-(4*(Ki - Ks)*(I - Ki*t))/S^2)*(Ki - Ks) .* ...
        %     (I - Ki*t) - t*(Ki - Ks)^2.*(exp(-(4*(Ki - Ks) * ...
        %     (I - Ki*t))/S^2) + 1) ) ./ ( (Ki - Ks)*(exp(-(4* ...
        %     (Ki - Ks)*(I - Ki*t))/S^2) - 1) );
        % [limB0 limB2]

    case 2 % Time form
        % We know I - and simulate time (hour)
        t = Haverkamp_t(x,plugin); I = plugin.I(:);
        % Vrugt and Gao, VZJ, 2022      
        % Define α: nx1 vector as it depends on time (or I)
        a = exp(2*b*(Ks-Ki)*(I-Ki*t)/S^2); 
        % --> reduces to a = exp(2*b*Ks*I/S^2) as Ki = 0
        % Define common terms numerator/denominator, nx1 vectors
        c1 = a + b - 1; c2 = a + b + a*b - 1;
        % ∂I/∂S, Equation 29a
        dtdS = (S^2*c1.*log(c1/b) - 2*I.*a*b*Ks) ./(S*Ks^2*(b-1)*c1);
         % Compute ∂I/∂Ks, Equation 29b
        dtdKs = (I*Ks.*c2 - S^2*c1.*log(c1/b)) ./ (Ks^3*(b-1)*c1); 
        % Compute ∂I/∂ß, Equation 19c
        dtdb = ( I * Ks.*(a * b + b - 1) -.5*S^2*c1.*log(c1/b) ...
            -.5*(a-1)*(b-1)*S^2/b ) ./ (Ks^2*(b-1)^2*c1); 
        switch d
            case 2
                J = [dtdS, dtdKs]; 
            case 3
                J = [dtdS, dtdKs, dtdb];
        end
        % dtdS = @(I,S,Ks,B) ( 2*S*log(1/B * exp( (2*B*I*Ks)/S^2) + ...
        %     (B-1)/B ) ) / ( Ks^2 * (2*B-2)) - ...
        %     4*I.*exp(2*B*I*Ks/S^2) ./ (Ks*S*(2*B-2) * ...
        %     (1/B * exp(2*B*I*Ks/S^2) + (B-1)/B) );
        % dtdKs = @(I,S,Ks,B) I/(Ks^2*(B-1)) - ( 2*S^2*log(1/B * ...
        %     exp( (2*B*I*Ks)/S^2) + (B-1)/B ) ) ./ ...
        %     ( Ks^3 * (2*B-2)) + 2*I.*exp(2*B*I*Ks/S^2 ) ./ ...
        %     (Ks^2*(2*B-2)*(1/B * exp(2*B*I*Ks/S^2) + (B-1)/B) );
        % dtdB = @(I,S,Ks,B) I/(Ks*(B-1)^2) - ( 2*S^2*log(1/B * ...
        %     exp( (2*B*I*Ks)/S^2) + (B-1)/B ) ) / ...
        %     ( Ks^2 * (2*B-2)^2) - S^2 * (1/B^2 * ...
        %     exp(2*B*I*Ks/S^2 ) - 1/B + (B-1)/B^2 - ...
        %     2*I*Ks/(B*S^2).*exp(2*B*I*Ks/S^2)) ./ ...
        %     (Ks^2*(2*B-2)*(1/B * exp( (2*B*I*Ks)/S^2) + (B-1)/B ));
end

end
