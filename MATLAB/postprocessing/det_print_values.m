function [X_ticks,x_min,x_max] = det_print_values(X,nx)
% This function determines which values will be used for tick labels

% INPUT: X = N by p matrix with N entries of p variables
%        nx = scalar that determines how many ticklabels are useds
% OUTPUT: X_ticks = p by n vector with nx values for each p variables
%
p = size(X,2);              % How many variables
nbins = 2*nx;               % number of bins
x_min = nan(p,1); x_max = nan(p,1);
X_ticks = nan(p,nx);
idx = 2:2:nbins;
for i = 1:p
    [~,x] = histcounts(X(:,i),nbins);
    X_ticks(i,1:nx) = x(idx);
    x_min(i,1) = x(1); x_max(i,1) = x(nbins+1);
end
% % % % How many x-values should be printed on each axis?
% % % mnmx_f = deal(nan(p,2));    % minmax values of p variables
% % % num_x = 2*nx+1;
% % % switch nx
% % %     case 2 % num_x = 5;
% % %         idx_col = [1 3 5]; 
% % %     otherwise
% % %         idx_col = 1:2:num_x;
% % % end
% % % n_col = numel(idx_col);
% % % % Initialize X_ticks with nan
% % % X_ticks = nan(p,numel(idx_col));
% % % % Determine ranges for all axes
% % % mnmx_t = [ min(X)' max(X)' ];
% % % % Determine how many 10 multiples for each entry
% % % ord_p = order_no(mnmx_t);
% % % for i = 1:p
% % %     % Note you can try +1 instead of +0 --> changes x/y values at ticks
% % %     n1 = 10^(-ord_p(i,1)+1); mnmx_f(i,1) = floor(mnmx_t(i,1)*n1)/n1;
% % %     % Other option is to use "floor"
% % %     % Note you can try +1 instead of +0 --> changes x/y values at ticks
% % %     n2 = 10^(-ord_p(i,2)+1); mnmx_f(i,2) = ceil(mnmx_t(i,2)*n2)/n2;
% % % end
% % % % Now check that two columns of mm are different
% % % idx_zero = find(diff(mnmx_f,[],2) == 0);
% % % % Otherwise, reset mnmx values to original values
% % % mnmx_f(idx_zero,1:2) = mnmx_t(idx_zero,1:2);
% % % % Now extract min and max of each variable
% % % x_min = mnmx_f(:,1); x_max = mnmx_f(:,2);
% % % % Set small values of x_min to zero;
% % % for i = 1:p
% % %     if x_min(i) < 1e-3
% % %         if x_min(i) < (x_max(i)-x_min(i))
% % %             x_min(i) = 0;
% % %         end
% % %     end
% % % end
% % % % Get difference between x_max and x_min
% % % dx_mnmx = (x_max - x_min)/(2*nx+1);
% % % ord_dx = order_no(dx_mnmx);
% % % n_dx = 10.^(-ord_dx+1); 
% % % dx_new = floor(n_dx .* dx_mnmx)./n_dx; 
% % % % Now use dx_new to create xticks
% % % for i = 1:p
% % %     x_vals = x_min(i) : dx_new(i) : x_max(i);
% % %     X_ticks(i,1:n_col) = x_vals(idx_col);
% % % end