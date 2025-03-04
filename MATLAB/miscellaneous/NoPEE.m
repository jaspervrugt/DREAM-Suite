function [y_sig,tab,a,b] = NoPEE(y,u,m,verbose)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Estimate measurement error standard deviation of a time series of data             %%
%% Requirements [1] underlying data-generating process, h(t), is sufficiently smooth  %%
%%              [2] sampling interval is high compared to the time-scale of h(t)      %%
%%              [3] measurement errors have a constant or heteroscedastic variance    %%
%%                                                                                    %%
%% SYNOPSIS:    [y_sig,tab,a,b] = NoPEE(y);                                           %%
%%              [y_sig,tab,a,b] = NoPEE(y,u);                                         %%
%%              [y_sig,tab,a,b] = NoPEE(y,u,m);                                       %%
%%              [y_sig,tab,a,b] = NoPEE(y,u,m,verbose);                               %%
%%  where                                                                             %%
%%   y          [input] n x 1 vector with values of measured signal                   %%
%%   u          [input] difference operator applied u times    (def: u = 3)           %%
%%   m          [input] tab window size to left and right      (def: m = 10)          %%
%%   verbose    [input] plotting to screen                     (def: verbose = 1)     %%
%%   y_sig      [outpt] (n-u)x2 matrix with data pairs of y, measurement error sigma  %%
%%   tab        [outpt] (n-u)x2 matrix with sliding average of sigma (window: [m m])  %%
%%   a          [outpt] slope of the (Y,sig) function                                 %%
%%   b          [outpt] intercept of the (Y,sig) function                             %%
%%                                                                                    %%
%% LITERATURE                                                                         %%
%%   Vrugt, J.A., C.G.H. Diks, W. Bouten, H.V. Gupta, and J.M. Verstraten (2005),     %%
%%       Improved treatment of uncertainty in hydrologic modeling: Combining the      %%
%%       strengths of global optimization and data assimilation, Water Resources      %%
%%       Research, 41(1), W01017, doi:10.1029/2004WR003059                            %%
%%                                                                                    %%
%% EXAMPLE I: HOMOSCEDASTIC ERROR                                                     %%
%%   n = 5000; x = linspace(0,2*pi,n)';                                               %%
%%   sig_err = 0.05;                                                                  %%
%%   y = 2 + sin(x); err = sig_err * randn(n,1); y_err = y + err;                     %%
%%   plot(x,y_err,'r'); hold on; plot(x,y,'k');                                       %%
%%   [y_sig,tab,a,b] = NoPEE(y_err,3);                                                %%
%%   -> if scatter does not have a slope -> measurement error has a constant variance %%
%%      with sigma ~ sqrt(mean(y_sig(:,2).^2));                                       %%
%%   NOTE: value of a ~ 0 (as it should be) and b ~ 0.05 (= sig_err)                  %%
%%                                                                                    %%
%% EXAMPLE II: HETEROSCEDASTIC ERROR                                                  %%
%%   n = 5000; X = linspace(0,2*pi,n)';                                               %%
%%   a_err = 0.1;                                                                     %%
%%   y = 2 + sin(x); err = normrnd(0,abs(a_err * y)); y_err = y + err;                %%
%%   plot(x,y_err,'r'); hold on; plot(x,y,'k');                                       %%
%%   [y_sig,tab,a,b] = NoPEE(y_err,3);                                                %%
%%   -> if scatter does have a slope -> measurement error has a non-constant variance %%
%%      the least squares regression fit determines the measurement error model       %%
%%   NOTE: value of a ~ 0.10 (= a_err)                                                %%
%%                                                                                    %%
%% © Written by Jasper A. Vrugt, Dec. 2004 & modified thereafter                      %%
%% University of Amsterdam                                                            %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Check input arguments for their validity - report if default values used
if nargin < 1
    error(['ERROR: user should input n-vector with values of measured',...
        'signal']);
end
if nargin == 1
    if numel(y) < 100
        fprintf(['WARNING: results may be unreliable as data set, y',...
            ', has relatively few elements\n']);
    end
end
if nargin < 2
    fprintf(['WARNING: The difference operator is not defined - ',...
        'we resort to default of u = 3\n']);
    u = 3;
end
if nargin >= 2
    if isempty(u)
        error(['ERROR: u should not be empty but an integer larger ',...
            'than zero']);
    end
    if ~ismember(u,1:50)
        error('ERROR: u should be an integer and larger than zero');
    end
end
if nargin < 3
    fprintf(['WARNING: Window size, m, not specified - we resort to',...
        'default of m = 10\n']);
    m = 10;
end
if nargin >= 3
    if ~ismember(m,1:numel(y))
        error(['ERROR: m (window to left and right) should be an ',...
            'integer and larger than zero']);
    end
    if m > round(numel(y)/10)
        fprintf(['WARNING: m (window to left and right) set rather ',...
            'large (default: m = 10)']);
    end
end
if nargin < 4
    fprintf(['WARNING: Verbose not specified - we resort to default ',...
        'of verbose = 1 (printing to screen)\n']);
    verbose = 1;
end
if nargin >= 4
    if ~ismember(verbose,[0 1])
        error(['ERROR: verbose should be set to 0 (no figure printing)',...
            'or 1 (figure printing)']);
    end
end

% Apply differencing operator u successive times
dy = diff(y,u);
% What is the mean of the y-values corresponding to △y?
% u = 1: [1 2], [2 3], [3 4], [4 5], [5 6]        w = [ 1 1 ]; (equal weights)
% u = 2: [1 2 2 3], [2 3 3 4], [3 4 4 5]          w = [ 1 2 1 ];
% u = 3: [1 2 2 2 3 3 3 4], [ 2 3 3 3 4 4 4 5]    w = [ 1 3 3 1 ];
% etc. (thus for u = 3: y_1 and y_4 are used once, y_2 and y_3 are used 3 times
% Initialize w for u = 1; add zeros to match final size of w
w = [ 1 1 zeros(1,u-1) ];
for j = 2:u
    w_add = [ 0 w(1:end-1) ];
    w = w + w_add;
end
% Now normalize so that weights add up to one
weights = w./sum(w);
% Now compute the mean y values that correspond to △Y
y_m = nan(numel(dy),1);
for j = 1:numel(dy)
    % Inner product of w and corresponding data points
    y_m(j,1) = weights * y(j:j+u,1);
end

% MAYBE SIMPLER TO UNDERSTAND
% % dy = diff(y,u);                                 % Apply differencing operator u times
% % y_m = y;                                        % Initialize mean Y values
% % for delta = 1:u                                 % Now loop u times
% %     y2 = [y_m(1:end-1,1) y_m(2:end,1)];         % Combine Y_m values
% %     y_m = mean(y2,2);                           % Compute mean Y adjacent values
% % end

% Now calculate the std of the measurement error
sig = sqrt( 1/nchoosek(2*u,u) * dy.^2);

% Remove small measurement errors - subjective - can be removed
% id = sig > 1e-4; y_m = y_m(id); sig = sig(id); 
% Prepare return argument - sort first column for moving average in tab
y_sig = sortrows([ y_m sig ],1);

% Now make sure that Y_m elements are unique (for spline interpolation)
% [~,id] = unique(y_m); y_sig = [ y_m(id) sig(id) ];
% NOTE: Non-unique rows are removed - but better that we average sigma values
% of duplicate y_m values instead - this is important if length of y is small

% Now use moving average with window to get average estimates of sigma
tab = [ y_sig(:,1) sqrt(movmean(y_sig(:,2).^2,[m m])) ];
% Now fit line to tabular (averaged) data
P = polyfit(tab(:,1),tab(:,2),1); a = P(1); b = P(2);

% Print figures to screen?
if verbose == 1
    % Open figure
    figure('name',['Standard deviation of measurement error versus ',...
        'measured value'],'units','normalized','outerposition',[0 0 1 1])
    subplot(1,2,1)
    loglog(y_sig(:,1),y_sig(:,2),'r.','markersize',15); hold on
    set(gca,'fontsize',18);
    axis square
    xlabel('$y$','interpreter','latex','fontsize',21,...
        'Units','normalized','Position',[0.5, -0.1, 0],...
        'fontweight','normal');
    ylabel('$\widehat{\sigma}_{y} \; {\rm in \; units \; of} \; y$',...
        'interpreter','latex','fontsize',21,'Units','normalized',...
        'Position',[-0.13, 0.5, 0],'fontweight','normal');
    % Now add tab and regression line to figure
    plot(tab(:,1),tab(:,2),'k--','linewidth',3);
    % Add legend
    evalstr_u = strcat('$\; {\rm data:}','\; u\;=\;',num2str(u),'$');
    evalstr_m = strcat('$\; {\rm tab:}','\; m\;=\;',num2str(m),'$');
    legend({char(evalstr_u),char(evalstr_m)},'interpreter','latex',...
        'fontsize',21,'location','northwest','box','off');
    title(['${\rm MEASUREMENT \; ERROR \; USING \; A \; LOG-LOG \; ',...
        'SCALE}$'],'interpreter','latex','Units', 'normalized',...
        'Position',[0.5, 1.03, 0],'fontsize',18);   
    set(gca,'TickDir','out');
    set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 1]);
    % get handle to current axes 
    h1 = gca;
    % set box property to off and remove background color
    set(h1,'box','off','color','none')
    % create new, empty axes with box but without ticks
    h2 = axes('Position',get(h1,'Position'),'box','on','xtick',[],...
        'ytick',[]);
    % set original axes as active
    axes(h1)
    % link axes in case of zooming and make rounding box square
    linkaxes([h1 h2]); axis(h2,'square');
%     axis([min(0.5*min(Y_dY(:,1)),1.5*min(Y_dY(:,1))) 1.5*max(Y_dY(:,1)) ...
%         0.5*min(Y_dY(:,2)) 1.5*max(Y_dY(:,2))]); 
    % Plot new figure on linear scale
    subplot(1,2,2)
    plot(y_sig(:,1),y_sig(:,2),'r.','markersize',15); hold on
    plot(y_sig(:,1),a*y_sig(:,1) + b,'b-','linewidth',3);
    set(gca,'fontsize',18);
    axis square
    xlabel('$y$','interpreter','latex','fontsize',21,...
        'Units','normalized','Position',[0.5, -0.1, 0],...
        'fontweight','normal');
    ylabel('$\widehat{\sigma}_{y} \; {\rm in \; units \; of} \; y$',...
        'interpreter','latex','fontsize',21,'Units','normalized',...
        'Position',[-0.12, 0.5, 0],'fontweight','normal');
    % Add function
    a_str = num2str(a); b_str = num2str(b);
    if a < 0
        num_values_a = 6;
    else
        num_values_a = 5;
    end
    if b < 0
        num_values_b = 6;
    else
        num_values_b = 5;
    end
    a_str = a_str(1:min(numel(a_str),num_values_a));
    b_str = b_str(1:min(numel(b_str),num_values_b));
    if b > 0
        evalstr = strcat('$\;\; \widehat{\sigma}_{y}',{' '},...
            '=',{' '},a_str,'y',{' '},'+',{' '},b_str,'$');
    else
        evalstr = strcat('$\;\; \widehat{\sigma}_{y}',{' '},...
            '=',{' '},a_str,'y',{' '},'-',{' '},b_str(2:end),'$');
    end
    % Add legend
    legend({char(evalstr_u),char(evalstr)},'interpreter','latex',...
        'fontsize',21,'location','northwest','box','off');
    title(['${\rm MEASUREMENT \; ERROR \; USING \; A \; LINEAR \; ',...
        'SCALE}$'],'interpreter','latex','Units',...
        'normalized', 'Position',[0.5, 1.03, 0],'fontsize',18);    
%     axis([min(0.9*min(Y_dY(:,1)),1.1*min(Y_dY(:,1))) 1.1*max(Y_dY(:,1)) ...
%         0.9*min(Y_dY(:,2)) 1.1*max(Y_dY(:,2))]); 
    % Now determine minor ticks on y-axis
    set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 1]); %
    set(gca,'TickDir','out');
    % get handle to current axes
    h1 = gca;
    % set box property to off and remove background color
    set(h1,'box','off','color','none')
    % create new, empty axes with box but without ticks
    h2 = axes('Position',get(h1,'Position'),'box','on','xtick',[],...
        'ytick',[]);
    % set original axes as active
    axes(h1)
    % link axes in case of zooming and make rounding box square
    linkaxes([h1 h2]); axis(h2,'square');
else
    fprintf('WARNING: No printing to screen as verbose = 0\n');
end