% Runs all samples from SWIG data base
close all hidden; clc; clear all

save_fig = 0;
fontsize_axes = 18;
fontsize_labels = 18;
fig_5 = 0;

load SWIG_1D
% Rename content of file
data_SWIG = SWIG_1D;
% Unpack the data
n_soil = size(data_SWIG,1);
% Summarize data
for soil_type = 1:n_soil
    dat = data_SWIG{soil_type};
    Data(soil_type,1:4) = [ soil_type size(dat,1) max(dat(:,1)) max(dat(:,2)) ];
end
%ii_data = Data(:,3)>= 1; idx_data = Data(ii_data,1);
% Define bounds on parameters
%eta_min = [ 1e-5 1e-5 1e-5 ]'; eta_max = [ inf inf 2 ]';
%eta_min = [ 1e-3 1e-3 0.3 ]'; eta_max = [ inf inf 2 ]';
eta_min = [ 1e-3 1e-3 ]'; eta_max = [ inf inf ]';
% Number of parameters
d = 2;
% How many trials with LM from different starting points?
Ntrials = 10; nu = 2; %1.3; % nu = 10;
% Define beta values to use
Bfix = [0.001 0.1:0.1:0.9 0.95 1.05 1.1:0.1:1.9 1.999]; NB = numel(Bfix);
% Initialize return matrices
[opt_LM_SWIG,min95,max95,min90,max90] = deal(nan(n_soil,d,NB,2)); [n_LM_SWIG,SSR_LM_SWIG] = deal(nan(n_soil,NB,2));
LM_res_SWIG = nan(n_soil,4,NB,2); warning off

%% CHECK PAR_SWIG_BETA --> USES PARFLOW - MUCH FASTER
% Now loop over each soil and call DREAM Package
for approach = 1:2
    approach
    for uu = 1:numel(Bfix)
        Bfixed = Bfix(uu); uu
        for soil_type = 1:n_soil
            % Unpack data
            dat = data_SWIG{soil_type};
            % Define plugin structure
            n_data = size(dat,1);
            % Store time and infiltration
            t_meas = dat(1:n_data,1); I_meas = dat(1:n_data,2);
            % Execute Levenberg Marquardt method
            [eta_opt,SSR_opt,K] = least_squares_haverkamp(t_meas,I_meas,approach,Bfixed,eta_min,eta_max,nu,Ntrials);
            % How many trials converge to same optimum?
            [nr,Km,Kstd,Kopt] = analyze_opt(eta_opt,SSR_opt,K);
            % Store
            LM_res_SWIG(soil_type,1:4,uu,approach) = [ nr Km Kstd Kopt ];
            %% Store output of least_squares
            zz = find(SSR_opt == min(SSR_opt)); zz = zz(1);
            % Least squares values of LM method
            opt_LM_SWIG(soil_type,1:d,uu,approach) = eta_opt(zz,1:d);
            % SSR of least squares values of LM method
            SSR_LM_SWIG(soil_type,uu,approach) = SSR_opt(zz);
            % Store number of data points
            n_LM_SWIG(soil_type,uu) = n_data;
            % Compute uncertainty at optimum
            [std_opt,r90,r95] = det_uncertainty(I_meas,t_meas,Bfixed,eta_opt(zz,1:d),SSR_opt(zz),eta_min,eta_max,n_data,approach);
            % Now store 95%
            min95(soil_type,1:d,uu,approach) = r95(1:d,1)';
            max95(soil_type,1:d,uu,approach) = r95(1:d,2)';
            min90(soil_type,1:d,uu,approach) = r90(1:d,1)';
            max90(soil_type,1:d,uu,approach) = r90(1:d,2)';
        end
        sum(SSR_LM_SWIG(1:n_soil,1,uu,approach))
    end
end
% Save optimal values
save SWIG_646_beta opt_LM_SWIG SSR_LM_SWIG data_SWIG n_soil Data min95 max95 min90 max90 LM_res_SWIG

idx_notsel = [ 5 27 30 48 57 60 61 72 74 76 77 78 80 95 100 109 110 111 113 ...
    116 120 121 126 127 128 131 144 145 147 175 176 177 187 192 193 196  ...
    202 203 204 221 222 225 226 228 275 276 269 302 ...
    303 305 306 371 363 377 385 386 390 393 394 400 410 411 421 424 427 428 435 ...
    445 452 464 465 477 485 513 515 516 520 521 522 540 542 543 544 545 546 547 563 564 566 567 581 592 600 ...
    601:637 639 641 ];
idx_sel = [1:n_soil]; idx_sel(idx_notsel) = [];
idx_new = [ 2 7 30 48 60 69 72 76 95 100 109 110 120 145 146 177 200 221 222 225 289 303 391 445 496 ...
    520 530 546 547 600 601:605];
idx = [1:n_soil]; idx(idx_new) = [];

% Plot results of LM
for par = 1:3
    figure(par+200)
    % Plot all data
    plot(opt_LM_SWIG(1:n_soil,par,1),opt_LM_SWIG(1:n_soil,par,2),'ro','markersize',5,'markerfacecolor','w','linewidth',2);
    % Now only low RMSE
    RMSE(1:n_soil,1) = sqrt(SSR_LM_SWIG(1:n_soil,1,1) ./ Data(1:n_soil,2));
    RMSE(1:n_soil,2) = sqrt(SSR_LM_SWIG(1:n_soil,1,2) ./ Data(1:n_soil,2));
    %     %idx_RMSE = find(SSR_LM_SWIG(1:n_soil,1,1) < 0.1);
    idx_RMSE = find(RMSE(:,2) < 0.05);
    % idx_RMSE = idx_sel;
    idx_RMSE = idx;
    hold on
    plot(opt_LM_SWIG(idx_RMSE,par,1),opt_LM_SWIG(idx_RMSE,par,2),'ro','markersize',5,'markerfacecolor','r','linewidth',2);
    hold on
    set(gca,'fontsize',16);
    switch par
        case 1
            xlabel('${\rm Soil\;sorptivity}, S {\rm (cm/h^{1/2})}$','interpreter','latex','fontsize',18);
            ylabel('${\rm Soil\;sorptivity}, S {\rm (cm/h^{1/2})}$','interpreter','latex','fontsize',18);
        case 2
            xlabel('${\rm Saturated\;hydraulic\;conductivity}, K_{\rm s} {\rm (cm/h)}$','interpreter','latex','fontsize',18);
            ylabel('${\rm Saturated\;hydraulic\;conductivity}, K_{\rm s} {\rm (cm/h)}$','interpreter','latex','fontsize',18);
        case 3
            xlabel('${\rm Coefficient\;}, \beta (-)$','interpreter','latex','fontsize',18);
            ylabel('${\rm Coefficient\;}, \beta (-)$','interpreter','latex','fontsize',18);
    end
    a = axis;
    min_data = min(a(1:2));
    max_data = max(a(1:2));
    % 1:1 relationship
    plot([min_data max_data],[min_data max_data],'k','linewidth',1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 5: PLOT MEASURED VERSUS FITTED CUMULATIVE INFILTRATION DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx_not = [ 2 7 30 48 60 69 72 76 95 100 109 110 120 145 146 177 200 221 222 225 289 303 391 445 496 ...
    520 530 546 547 600 601:605];
idx_sel = [1:n_soil]; idx_sel(idx_not) = [];
idx_sel = idx_not;
if fig_5 == 1
    % Create figure name
    fig_name = 'Plot of observed and fitted data';
    % Define maximum of
    %    max_I =   [ 56 75 280 4000 8000 35 350 1200 75 120 8   20 ];
    %max_x = 240;
    %    delta_I = [ 14 15  70 1000 2000 7  70  300  15 30  2   5 ];
    %    minortick_I = [ 2 3 10 200 400  1  10  50   3  5  0.4  1 ];
    idx = [ 0 1 2 3 0 1 2 3 0 1 2 3 ];
    idy = [ 0 0 0 0 1 1 1 1 2 2 2 2 ];
    ct = 1;
    n_graphs = 12;
    fig_num = 1; rem_fig = 1:n_graphs:n_soil;
    for soil_type = 1:numel(idx_sel); n_soil
        % Create figure name
        if rem(soil_type,rem_fig(fig_num)) == 0
            fig_name = char(strcat('Plot of observed and fitted cumulative infiltration:',{' '},num2str(fig_num)));
            % Create figure
            figure('unit','inches','name',fig_name,'PaperOrientation','portrait','position',[0.1 0.4 19.7 14]);
            % Set cf to 1;
            cf = 1; fig_num = fig_num + 1; fig_num = min(fig_num,ceil(n_soil/n_graphs));
        end
        dat = data_SWIG{idx_sel(soil_type)};
        fx_LM = []; plugin.t = dat(:,1);
        % Approach 1 - predict I
        fx_LM(:,1) = Haverkamp_I(opt_LM_SWIG(idx_sel(soil_type),1:3,1),plugin);
        % Approach 2 - predict time
        plugin.I = dat(:,2);
        fx_LM(:,2) = Haverkamp_t(opt_LM_SWIG(idx_sel(soil_type),1:3,2),plugin);
        % define max_x
        x_max = max(dat(:,1)); y_max = max(dat(:,2));
        % define axis position
        ax_str = strcat('ax',num2str(soil_type));
        evalstr = strcat('ax',num2str(soil_type),{' '},'= axes(''units'',''inches'');'); eval(char(evalstr));
        evalstr = strcat('axpos',num2str(soil_type),{' '},'= [ 1.25 + idx(cf)*4.75 , 7.7-idy(cf)*3.4 3.8 2.65 ];'); eval(char(evalstr));
        evalstr = strcat('set(','ax',num2str(soil_type),',','''position''',',axpos',num2str(soil_type),');'); eval(char(evalstr));
        % plot
        plot(eval(ax_str),dat(:,1),dat(:,2),'ro','markersize',5,'linewidth',1.5,'markerfacecolor',[1 1 1]); hold on
        set(gca,'fontsize',fontsize_axes);
        % Now plot optimum of DREAM and LM
        for u = 1:2
            switch u
                case 1 % I optimized LM
                    plot(eval(ax_str),dat(:,1),fx_LM(:,1),'b--','linewidth',1.5);
                case 2 % t optimized LM
                    plot(eval(ax_str),fx_LM(:,2),dat(:,2),'k--','linewidth',1.5);
            end
        end
        % Add labels on x and y-axis
        if ismember(cf,9:12)
            tx = xlabel('${\rm Time}, t \; {\rm (\,h\,)}$','interpreter','latex',...
                'fontsize',fontsize_labels,'Units','normalized','Position',...
                [0.5, -0.3, 0],'fontweight','normal');
            pos_tx = get(tx,'position'); set(tx,'position',[ pos_tx(1) pos_tx(2)+0.08 pos_tx(3)]);
        end
        if ismember(cf,[1 5 9])
            ty = ylabel('$\widetilde{I} \; {\rm (cm)}$','interpreter',...
                'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                [-0.23, 0.5, 0],'fontweight','normal');
            pos_ty = get(ty,'position');
            if ismember(cf,[1 9])
                set(ty,'position',[ pos_ty(1)+0.08 pos_ty(2:3)]);
            else
                set(ty,'position',[ pos_ty(1)+0.005 pos_ty(2:3)]);
            end
        end
        % set axis
        axis(eval(ax_str),[0 x_max 0 y_max]);
        % Add label with name of soil
        str = data_SWIG{idx_sel(soil_type),3};
        % find upper case letter
        idx_str = isstrprop(str,'upper'); ii = find(idx_str==1); n_str = numel(ii);
        % change to string of individual letters
        str_char = num2str(str);
        % check if SandyClayloam
        if strcmp(str_char,'SandyClayloam'), str_f = 'Sandy clay loam'; n_str = 4; end
        % reorganize string
        switch n_str
            case 1 % Single capital letter
                str_f = str_char(1:end);
            case 2 % Two capital letters
                str_f = strcat(str_char(1:ii(2)-1),{' '},lower(str_char(ii(2):end)));
            case 3 % Three capital letters
                str_f = strcat(str_char(1:ii(2)-1),{' '},lower(str_char(ii(2):ii(3)-1)),{' '},lower(str_char(ii(3):end)));
        end
        axis tight; a = axis;
        %        axis([0 x_max 0 y_max]);
        %        a = axis;
        str_f = char(str_f);
        text(eval(ax_str),a(1)+0.5*(a(2)-a(1)),a(3)+1.01*(a(4)-a(3)),str_f,'interpreter','latex','fontweight',...
            'bold','fontsize',fontsize_labels,'HorizontalAlignment', 'center');
        text(eval(ax_str),a(1)+0.02*(a(2)-a(1)),a(3)+1.01*(a(4)-a(3)),strcat(num2str(soil_type),':',num2str(cf),')'),'fontsize',fontsize_labels,'interpreter','latex');
        str_s = data_SWIG{idx_sel(soil_type),2};
        str_sample = char(strcat('${\rm Code}:',{' '},num2str(str_s),'$'));
        text(eval(ax_str),a(1)+0.96*(a(2)-a(1)),a(3)+0.05*(a(4)-a(3)),str_sample,'interpreter','latex','fontweight',...
            'bold','fontsize',fontsize_labels,'HorizontalAlignment', 'right');
        % Now determine xticks and xticklabels and minorxticks
        xticklabel = get(eval(ax_str),'xticklabel');
        xtick = get(eval(ax_str),'xtick');
        % remove current labels - and replace with own values closer to axis
        set(eval(ax_str),'xticklabel',[]);
        % now plot manually - for all except 12
        dy = a(3) - 0.10*(a(4)-a(3));
        for rr = 1:numel(xtick)
            text(eval(ax_str),xtick(rr),dy,xticklabel(rr),'fontsize',fontsize_axes,'HorizontalAlignment', 'center');
        end
        %        xtickformat('%.0f')
        ytickformat('%.1f');
        
        %       text(eval(ax_str),0.02*x_max,0.95*y_max,char(index(soil_type)),'interpreter','latex','fontweight',...
        %           'normal','fontsize',fontsize_labels); %max_x = 240;
        % Add t_char - 2.5 percentile
        set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 1]); %
        set(gca,'TickDir','out');
        %xtickformat('%.0f')
        ytickformat('%.0f')
        % set axis
        axis(eval(ax_str),[0 x_max 0 y_max]);
        % get handle to current axes (next lines remove xticks on top line of figure
        h1 = gca;
        % set box property to off and remove background color
        set(h1,'box','off','color','none')
        % create new, empty axes with box but without ticks
        h2 = axes('Position',get(h1,'Position'),'box','on','xtick',[],'ytick',[]);
        % set original axes as active
        axes(h1);
        % link axes in case of zooming and make rounding box square
        linkaxes([h1 h2]);
        axis(eval(ax_str),[0 x_max 0 y_max]);
        % Update cf
        cf = cf + 1;
        % check before we move to new figure
        if rem(soil_type+1,rem_fig(fig_num)) == 0
            set(gcf,'color','w');
            % set(legb(1),'Color','r'); set(legb(4),'linewidth',1.5,'markersize',8); set(legb(5),'XData',[0.0603 0.1838],'linewidth',3); set(legb(2),'Color','b');
            evalstr = strcat('export_fig figure_SWIG_',num2str(fig_num-1),{' '},'-','pdf');
            % now save file
            if save_fig == 1, eval(char(evalstr)); end
        end
        
    end
end
