function [log_L] = DCT_GPR(x)
% DCT inversion of geophysical travel time data measured with GPR

% Store local variables in memory
persistent func dummy

% Now define the local variables - only once
if isempty(dummy)
    
    % Standard deviation of Gaussian traveltime error
    func.error = 0.5;
    % Inversion parameters in x-direction (DCT order)
    func.parx = 8;
    % Inversion parameters in z-direction (DCT order)
    func.parz = 8;
    % Now setup forward problem
    [ func ] = setup_GPR ( func );
    % Update dummy
    dummy = 1;
    
end

% -------------------------------------------------------------------------
%                               Model script
% -------------------------------------------------------------------------

% Initialize DCT model
model_DCT = zeros(func.dimver,func.dimhor);

% Assign proposed model
for k = 1 : func.parz
    for i = 1 : func.parx
        model_DCT(k,i) = x ( ( k - 1 ) * func.parx + i ); 
    end
end
% Do inverse DCT
model = idct2(model_DCT);           
% Transform from logarithmic scale;
model = 10.^model;
% Transpose
model = model';
% Find minimum and maximum value
minvalue = min(min(model)); maxvalue = max(max(model));
% Calculate residual
data = func.J*model(:)-func.datasim;  
% Calculate weighted residual
data = data./func.error;              
% Calculate log-likelihood
log_L = sum(data.^2);                 
% Penalty of models outside range
if minvalue<(1/0.17)
    log_L = log_L*(1+(1/0.17)/minvalue); 
end
% Penalty of models outside range
if maxvalue>(1/0.05)
    log_L = log_L*(1+maxvalue/(1/0.05));
end
% Return log-likelihood
log_L = -1/2*log_L;