function [plugin,mult,mult_names,fpar_mult] = check_rainfall(plugin)
% Locates individual storm events in rainfall-discharge record 

% Analyze original rainfall and discharge record
plugin.id = analyze_rainfall(plugin.P);
% How many multiplers?
plugin.n_mult = size(plugin.id,1);
% Now assign the ranges of the multiplier
mult.min = 0.05 * ones(1,plugin.n_mult); mult.max = 5 * ones(1,plugin.n_mult);
% Now determine names of multipliers
mult_names = cell(1,plugin.n_mult);
for ii = 1:plugin.n_mult
    mult_names{ii} = strcat('\beta_{',num2str(ii),'}');
end
fpar_mult = ones(1,plugin.n_mult);

end

% Secondary function
function [id,P_new] = analyze_rainfall(P)
% Identifies storm events from precipitation data record only
% Option (commented out) also allows to jointly use discharge as well
% Indeed: discharge can go up while rainfall is zero --> suspicious

% Length of P?
n = size(P,1);
% Take the difference of the discharge values
%delta = [Y(2) - Y(1), diff(Y)];
% Set counter and flag
ct = 1; flag = 0;

% Loop over each day of precipitation data set
for j = 1:n
    % Check whether rainfall is larger than zero
    if ( P(j,1) > 0 )
        if ( flag == 0 )
            % set flag to one as new storm just started
            flag = 1; 
            % Set idx_start
            id_s(ct,1) = j;
        end
    elseif ( P(j,1) == 0 ) && ( flag == 1 )
%    if (( P(j,1) == 0 ) & ( delta(j,1) <= 0 )) & (flag == 1),
%    if ( P(j,1) == 0 ) && ( flag == 1 )
        % end of storm event
        id_e(ct,1) = j - 1;
        % update counter
        ct = ct + 1;
        % set flag to zero
        flag = 0;
    end
end

% Make sure that final idx_end is correct
if P(n,1) > 0
    id_e(ct,1) = n;
end

% Return both idx_start and idx_end in one matrix
id = [id_s id_e]; P_new = P;

end

% % Adjust rainfall
% P_new = P;
% 
% % Adjust the precipitation values if zero --> 
% for j = 1:size(idx,1),
%     % Determine the index
%     ii = [idx(j,1):idx(j,2)];
%     % Make sure that all events have positive rainfall!
%     R = P(ii); 
%     % Check whether any of them is zero
%     uu = find(R > 0); 
%     % If uu is empty (discharge increasing but rain is zero!)
%     if isempty(uu), 
%         % Then set minR to mean rainfall
%         minR = mean(P); 
%     else
%         % Set minR to minimum of the daily rainfall amounts within the storm
%         minR = min(R(uu)); 
%     end
%     % Check whether there is zero rainfall
%     uu = find(R == 0);
%     % If so then adjust rainfall
%     Pnew(ii(uu),1) = minR;
% end