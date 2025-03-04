function [E,q,p_0] = efdc(ID,type,col)
% Empirical flow duration curve

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Written by Jasper A. Vrugt                          %
%                   University of California Irvine                       %
%                           December 2014                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isnumeric(ID)    % Do not do anything - just rename ID
    Q = ID;
elseif ischar(ID)   % Load the data
    Q = load_data_dly(ID);
end

% Take only zero or postive streamflow values ( should be all values)
ii = Q(:,col) >= 0; Q = Q(ii,col);
% size of all data
n = size(Q,1);

% Now return FDC
switch type
    case 'daily'    % Daily FDC
        p_0 = sum(Q==0)/n;                      % probability of zero flows
        Q = Q(Q(1:n,1) > 0,1:end);              % only positive flows
        q = Q(1:end,1);                         % define Y
    case 'weekly'   % Weekly FDC
        ct = 1; p = 0;                          % set initial counter
        for i = 1:7:n                           % average weekly data
            y = Q(i:min(i+6,n),1);              % extract data
            y = y(y>0);                         % remove negative values
            if ~isempty(y)
                q(ct,1) = mean(y);              % average streamflow values
                ct = ct + 1;                    % update counter
            else
                p = p + 1;
            end
            p_0 = p/(ct-1);                     % p_0 = probability of zero weekly flows
        end
    case 'monthly'  % Monthly FDC
        ct = 1; p = 0;                          % Set initial counter
        id_mth = [1;find(diff(Q(:,2))~=0)+1];   % Find months
        for i = 2:size(id_mth,1)                % Average weekly data
            y = Q(id_mth(i-1):id_mth(i)-1,1);   % Extract data
            y = y(y>0);                         % Remove negative values
            if ~isempty(y)                      % Average streamflow values
                q(ct,1) = mean(y);
                ct = ct + 1;                    % Update counter
            else
                p = p + 1;
            end
            p_0 = p/(ct-1);                     % p_0 = probability of zero monthly flows
        end
    case 'yearly'   % Yearly FDC
        ct = 1;                                 % Set counter
        id_yr = [1;find(diff(Q(:,1))~= 0)+1];   % Find months
        for i = 2 : size(id_yr,1)               % Average weekly data
            y = Q(id_yr(i-1):id_yr(i)-1,1);     % Extract data
            y = y(y>0);                         % Remove negative values
            if ~isempty(y)
                q(ct,1) = mean(y);              % Average streamflow values
                ct = ct + 1;                    % Update counter
            else
                p = p + 1;
            end
            p_0 = p/(ct-1);                     % p_0 = probability of zero yearly flows
        end
    otherwise
        error('wrong type --> should be daily/weekly/monthly/yearly');
end

% How many values of q?
n = size(q,1);
% Now sort the data
q_s = sort(q);
% Calculate the exceedance probabilities
Ep = 1 - ((1:n)'-0.5) ./ n;
% Now return streamflow values and corresponding exceedance probabilities
FDC = [ Ep , q_s ];
% Exceedance probabilities and streamflow values
E = FDC(:,1); q = FDC(:,2);