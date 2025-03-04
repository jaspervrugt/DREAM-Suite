function fix_ticklabels(ax,axis_lab,n_add)
% Fix the ticklabels on x and/or y-axis in plotmatrix

if nargin < 2
    axis_lab = 'x'; % x-axis
end
if nargin < 3
    n_add = 0; % adds an extra digit
end

% Hold current axes
hold(ax,'on');
% Extract from current axes the xticklabels or yticklabels
str_lab = eval(char(strcat(axis_lab,'ticklabels(ax)')));
% How many labels
n = size(str_lab,1);
% Initialize number of entries before and after the dot
[n_bdot,n_adot] = deal(nan(1,n));
% How many values before and after the dot?
for z = 1:n
    lab = char(str_lab(z,:));
    % It is possible that MATLAB has tick entries that are empty
    switch isempty(lab)
        case 0
            idx_dot = strfind(lab,'.');
            idx_pnlty = 0; if strcmp(lab(1),'-'), idx_pnlty = 1; end
            % check how many empty spaces
            idx_empty = strfind(lab,' ');
            if isempty(idx_dot)
                n_bempty = numel(idx_empty);
                n_bdot(z) = numel(lab) - idx_pnlty - n_bempty; n_adot(z) = 0;
            else
                n_bempty = sum(idx_empty < idx_dot); 
                n_aempty = sum(idx_empty > idx_dot);
                n_bdot(z) = idx_dot - 1 - idx_pnlty - n_bempty; 
                n_adot(z) = numel(lab) - n_bdot(z) - idx_pnlty - 1 - n_aempty;
            end
        case 1
            % do nothing
    end
end
n_tot = max(n_bdot + n_adot) + n_add; n_adot = max(n_adot) + n_add;

% Now adjust format on xticklabels or yticklabels
eval(char(strcat(axis_lab,'tickformat(ax,''%',num2str(n_tot),'.',...
    num2str(n_adot),'f''',')')));