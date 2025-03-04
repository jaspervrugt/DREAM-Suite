function H = DTB_model(x)
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> %
% <><><><><>  DTB_model - Simulation of the depth to bedrock  <><><><><> %
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> %
%                                                 Version:   29.Jan.2016 %
%                                                    (c) Guilherme Gomes %
%                                              guilhermejcg@yahoo.com.br %
%                                            Modified by Jasper A. Vrugt %
% This program simulates the bedrock depth                               %
%                                                                        %
% SYNOPSIS                                                               %
%                                                                        %
% function H = DTB_model(x)                                              %
%                                                                        %
% Input:  x       vector of parameters (Phi, lambda_1, Sc, lambda_2)     %
% Output: H       vector of simulated bedrock depths                     %
%                                                                        %
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> %

% Retain variables in local memory
persistent Ld Zx Zy Ld_n nabla_z id

if isempty(Ld)
    data = load('Input_data.txt');                      % load input data
    Zx = reshape(data(2:end,1),data(1,1),data(1,2));    % slope gradient x-direction
    Zy = reshape(data(2:end,2),data(1,1),data(1,2));    % slope gradient y-direction
    Ld = reshape(data(2:end,3),data(1,1),data(1,2));    % drainage distance
    Ld_n = (Ld - min(Ld(:)))./(max(Ld(:))-min(Ld(:)));  % normalized drainage distance
    max(Ld(:)), min(Ld(:))
    Ld(1,:)
    size(Ld)
    nabla_z = Zx.^2 + Zy.^2;                            % gradient norm
    DTBdata = load('DTB_data.txt'); id = DTBdata(:,2);  % define calibration points
    % NOTE: Can use 6th variable "plugin" in example_25 
    % instead to port all this stuff as second variable 
    % to this DTB_model function
end

Lambda = exp( - x(2) * (1 - Ld_n).^x(4) );              % bedrock valley shape term
Psi = 1 - ( nabla_z ./ x(3)^2 );                        % threshold angle of mass movement
H = ( Psi ./ Lambda ) .* sqrt ( x(1) * Ld.^2 );         % Modeled depth to bedrock
H = H(id);                                              % return simulated data

end