%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Plot Jacobian entries of Haverkamp infiltration equation              %%
%% Figures 6 and 10 of the following paper                               %%
%%  Vrugt, J.A. and Y. Gao (2022), On the three-parameter infiltration   %%
%%      equation of Parlange et al. (1982): Numerical solution,          %%
%%      experimental design, and parameter estimation, Vadose Zone       %%
%%      Journal, 21:e20167, pp. 1-25, https://doi.org/10.1002/vzj2.20167 %%
%%                                                                       %%
%% © Written by Jasper A. Vrugt, Feb 2019                                %%
%% University of California Irvine                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

close all hidden; clc; clear

% FIGURES
fig_6 = 1;      % Jacobian matrices of time solution Haverkamp
fig_9 = 1;      % Example of residuals of infiltration and time forms
fig_10 = 1;     % Jacobian matrices of infiltration solution Haverkamp
save_fig = 0;   % save figures or not

color_conf_line = [0.6 0.6 0.6];
color_scatter = [0 102 255]/255;
color_LMtime = [ 96 168 48]/255; % [113, 123, 40]/255;
fontsize_axes = 21; fontsize_labels = 21;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                            FIGURES                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 6: JACOBIAN OF INFILTRATION FORM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if fig_6 == 1
    % Define default parameter values
    S = 2; Ks = 1; B = 1.5; Ki = 1e-1; eta = [S Ks B Ki];
    I1 = [0.0:0.1:20]'; plugin.I = I1;
    % Do forward run to determine Infiltration corresponding to
    [t,flag] = Haverkamp_t(eta,plugin);
    def1 = t; plugin.t = t;
    % Compute Jacobian with kappa (delta) = 0.01 - central differences
    J_anal = jac_Haverkamp_anal(eta,plugin,1);
    I2 = [0:1:20]'; plugin.I = I2;
    [t,flag] = Haverkamp_t(eta,plugin); plugin.t = t;
    J_num = jac_Haverkamp_num(eta,plugin,1);
    %    J = jac_Hav(eta,plugin);
    x_max = max(I1);
    %    y_max = [ 0.2 1.5 0.6 ];
    %    y_min = [ -1.8 -18 -0.03 ];
    y_max = [  2.0   18    0.08  2.20 ];
    y_min = [ -0.1 -1.0  -0.70 -0.08 ];
    % idx and idy
    id_x = [ 0 1  0 1 ];
    id_y = [ 0 0  1 1 ];
    % Create figure name
    fig_name = 'Figure 6 of VZJ paper: Jacobian of infiltration equation - infiltration form';
    % Create figure
    figure('unit','inches','name',fig_name,'PaperOrientation', ...
        'portrait','position',[0.3 0.5 12 14],'numbertitle','off');
    % Define index
    index = {'(A)','(B)','(C)','(D)','(E)','(F)','(G)','(H)','(I)','(J)','(K)','(L)','(M)'};
    % Plot 3 figures
    for zz = 1:4
        ax_str = strcat('ax',num2str(zz));
        evalstr = strcat('ax',num2str(zz),{' '},'= axes(''units'',''inches'')'); eval(char(evalstr));
        evalstr = strcat('axpos',num2str(zz),{' '},'= [ 1.35 + id_x(zz)*5.3 , 5.8-id_y(zz)*5 3.8 3.8 ]'); eval(char(evalstr));
        evalstr = strcat('set(','ax',num2str(zz),',','''position''',',axpos',num2str(zz),')'); eval(char(evalstr));
        switch zz
            case 1
                plot1 = plot(eval(ax_str),I1,J_anal(:,1),'k','linewidth',2);
                ylabel(eval(ax_str),'$\partial I/\partial S\;\; {\rm (h^{1/2}/cm)}$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [-0.18, 0.5, 0],'fontweight','normal');
            case 2
                plot(eval(ax_str),I1,J_anal(:,2),'k','linewidth',2);
                ylabel(eval(ax_str),'$\partial I/\partial K_{\rm s}\;\; {\rm (h)}$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [-0.16, 0.5, 0],'fontweight','normal');
            case 3
                plot(eval(ax_str),I1,J_anal(:,3),'k','linewidth',2);
                ylabel(eval(ax_str),'$\partial I/\partial \beta\;\; {\rm (cm)}$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [-0.2, 0.5, 0],'fontweight','normal');
            case 4
                plot(eval(ax_str),I1,J_anal(:,4),'k','linewidth',2);
                ylabel(eval(ax_str),'$\partial I/\partial K_{\rm i}\;\; {\rm (h)}$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [-0.18, 0.5, 0],'fontweight','normal');
        end
        hold on
        % Now plot the numerical Jacobian
        plot2 = plot(eval(ax_str),I2,J_num(:,zz),'ro','color',color_scatter,'markerfacecolor','w','linewidth',2); hold on
        axis([0 x_max y_min(zz) y_max(zz) ]);
        % Add (A), (B) and (C)
        set(gca,'fontsize',fontsize_axes);
        % Add labels on x and y-axis
        xlabel(eval(ax_str),'${\rm Cumulative\;infiltration}, I \; {\rm (\,{\rm cm}\,)}$','interpreter','latex',...
            'fontsize',fontsize_labels,'Units','normalized','Position',[0.5, -0.14, 0],'fontweight','normal');
        a = axis;
        set(eval(ax_str),'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 1]); %
        set(eval(ax_str),'TickDir','out','clipping','off');
        % Now determine xticks and xticklabels and minorxticks
        xticklabel = get(eval(ax_str),'xticklabel');
        xtick = get(eval(ax_str),'xtick');
        % remove current labels - and replace with own values closer to axis
        set(eval(ax_str),'xticklabel',[]);
        % now plot manually - for all except 12
        dy = a(3) - 0.07*(a(4)-a(3));
        for rr = 1:numel(xtick)
            text(eval(ax_str),xtick(rr),dy,xticklabel(rr),'fontsize',fontsize_axes,'HorizontalAlignment', 'center');
        end
        %        xtickformat('%.0f')
        switch zz
            case {1,3,4}
                ytickformat('%.1f')
            case {2}
                ytickformat('%.0f')
        end
        % Now create minor ticks - y axis
        %        eval(strcat(ax_str,'.XAxis.MinorTickValues = xminortick'));
        axis(eval(ax_str),[0 x_max a(3) y_max(zz) ]);
        % set box property to off and remove background color
        set(eval(ax_str),'box','off','color','none')
        % create new, empty axes with box but without ticks
        h2 = axes('Position',get(eval(ax_str),'Position'),'box','on','xtick',[],'ytick',[]);
        % set original axes as active
        % axes(eval(ax_str))
        % link axes in case of zooming and make rounding box square
        linkaxes([eval(ax_str) h2]);
        % add legend
        if zz == 1
            [lega,legb,legc,legd] = legend([plot1 plot2],{'${\rm \;\;Analytic}$','${\rm \;\;Numerical}$'},'interpreter','latex','box','off','location',...
                'southeast','fontsize',fontsize_labels);
            set(legb(1),'Color','k'); set(legb(5),'linewidth',1.5,'markersize',8); set(legb(2),'Color',color_scatter);
            set(legb(3),'linewidth',2.5);
            % set(eval(ax_str),'PaperPositionMode','auto');
            set(gcf, 'Renderer', 'painters')
        end
        axis([0 x_max y_min(zz) y_max(zz) ]);
        a = axis;
        % Add label (A), (B), or (C)
        text(eval(ax_str),a(1)+0.02*(a(2)-a(1)),a(3)+0.955*(a(4)-a(3)),char(index(zz)),...
            'fontsize',fontsize_labels+2,'interpreter','latex');
        set(gcf,'color','w');
       
        if zz < 3
            % add 2nd x-axis with time
            yloc = a(3)+1.02*(a(4)-a(3));
            line(eval(ax_str),[0 x_max],[yloc yloc],'color','k','linewidth',0.5);
            % now set vertical lines at right places
            mult = x_max / max(def1); ttick = [ 0 4 8 12 16 ]; tticklabel = {'0','4','8','12','16'};
            dy2 = yloc + 0.072*(a(4)-a(3));
            for rr = 1:numel(ttick)
                line(eval(ax_str),[ttick(rr) ttick(rr)]*mult,[yloc yloc+0.03*(a(4)-a(3))],'color','k','linewidth',0.5);
                % add label
                text(eval(ax_str),mult*ttick(rr),dy2,tticklabel(rr),'fontsize',fontsize_axes,'HorizontalAlignment', 'center');
            end
            % minor ticks
            minorttick = [0:18]; minorttick = minorttick(~ismember(minorttick,ttick)); %[ 1 2 3 5 6 7 9 10 11 13 14 15 17 18 ];
            for rr = 1:numel(minorttick)
                line(eval(ax_str),[minorttick(rr) minorttick(rr)]*mult,[yloc yloc+0.015*(a(4)-a(3))],'color','k','linewidth',0.5);
            end
            % add axes label
            dy3 = yloc + 0.15*(a(4)-a(3));
            text(eval(ax_str),10,dy3,'${\rm Time}, t{\rm \;(h)}$','fontsize',fontsize_labels,'HorizontalAlignment', 'center','interpreter','latex');
        end
    end
    %  set(legb(1),'Color','k'); set(legb(3),'linewidth',2.5); set(legb(4),'markersize',8); set(legb(2),'Color','r');
    % add a line at bottom with infiltration times
    %  set(gca,'clipping','off');
    set(legb(1),'Color','k'); set(legb(5),'linewidth',1.5,'markersize',8); set(legb(2),'Color',color_scatter);
    set(legb(3),'linewidth',2.5);
    % addpath to get nicest pdf from printing to pdf within matlab
    addpath('C:\Jasper Vrugt\Papers\Book\Chapter 1\MATLAB\Figure_resolution');
    % now use export_fig function from Figure_resolution directory
    evalstr = strcat('export_fig figure_jac_infil',{' '},'-','pdf');
    % now save file
    if save_fig == 1, eval(char(evalstr)); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 9: RESIDUAL EXAMPLE OF INFILTRATION AND TIME FORM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Selection from data comparison after the fact (LM_SWIG figures)
if fig_9 == 1
    plugin.t = [ 0 : 0.1 : 1 ]';
    x_min = 0; y_min = 0;
    x_max = max(plugin.t); y_max = 2.6;
    % Create figure name
    fig_name = 'Figure 9 of VZJ paper: Residuals of infiltration form of Haverkamp';
    eta = [ 2 1 1.5 0.1 ]'; 
    color_leg = {'color_scatter','color_LMtime'};
    % Create figure (change zz = 1:1 and zz = 2:2 for plotting)
    for zz = 2:2
        figure('unit','inches','name',fig_name,'PaperOrientation', ...
            'portrait','position',[0.3 0.5 8.5 6],'numbertitle','off');
        ax_str = strcat('ax',num2str(zz));
        evalstr = strcat('ax',num2str(zz),{' '},'= axes(''units'',''inches'');'); eval(char(evalstr));
        evalstr = strcat('axpos',num2str(zz),{' '},'= [ 1.35 1.2 6.5 4.5 ];'); eval(char(evalstr));
        evalstr = strcat('set(','ax',num2str(zz),',','''position''',',axpos',num2str(zz),');'); eval(char(evalstr));
        % Create data and imulated output
        I_dat = Haverkamp_I(eta,plugin); % Data
        plugin.I = I_dat;
        I_dat = I_dat + normrnd(0,0.25,numel(plugin.t),1);
        I_mod = Haverkamp_I([2 1 1.3],plugin); % Model        
        t_dat = Haverkamp_t(eta,plugin); 
        t_dat = t_dat + normrnd(0,0.075,numel(plugin.t),1);
        t_mod = Haverkamp_t([2 1 1.3],plugin); % Model
        % Plot residuals
        for i = 2:numel(plugin.t)
            switch zz
                case 1
                    plot3 = line(eval(ax_str),[plugin.t(i) plugin.t(i)],[I_dat(i) I_mod(i)],...
                        'color',color_conf_line,'linewidth',0.1); 
                case 2
                    plot3 = line(eval(ax_str),[t_dat(i) t_mod(i)],[plugin.I(i) plugin.I(i)],...
                        'color',color_conf_line,'linewidth',0.1);
            end
            if i == 2; hold on; end
        end        
        % plot simulated and measured data
        switch zz
            case 1
                % Plot data
                plot1 = plot(eval(ax_str),plugin.t(2:end),I_dat(2:end),'ro',...
                    'linewidth',1.5,'markersize',6,'markerfacecolor','w');
                % Rerun model to get nice curve
                plugin.t = [0:0.001:1]'; I_nice = Haverkamp_I([2 1 1.3],plugin);
                plot2 = plot(eval(ax_str),plugin.t,I_nice,'color',color_scatter,...
                    'linewidth',1.5);
            case 2
                % Plot data
                plot1 = plot(eval(ax_str),t_dat(2:end),plugin.I(2:end),'ro','color','r',...
                    'linewidth',1.5,'markersize',6,'markerfacecolor','w');
                % Rerun model to get nice curve
                plugin.t = [0:0.001:1]'; I_nice = Haverkamp_I([2 1 1.3],plugin);
                plugin.I = I_nice; t_nice = Haverkamp_t([2 1 1.3],plugin);
                plot2 = plot(eval(ax_str),t_nice,I_nice,'color',color_LMtime,...
                    'linewidth',1.5);
        end
        % Add labels on x and y-axis
        xlabel(eval(ax_str),'${\rm Time},t\;({\rm h})$','interpreter',...
            'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
            [ 0.5 , -0.14, 0],'fontweight','normal');
        ylabel(eval(ax_str),'${\rm Cumulative\;infiltration},I\;({\rm cm})$','interpreter',...
            'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
            [ -0.12, 0.50, 0],'fontweight','normal');
        axis(eval(ax_str),[x_min x_max y_min y_max ]);
        % Add (A), (B) and (C)
        set(gca,'fontsize',fontsize_axes);
        % Add labels on x and y-axis
        a = axis;
        set(eval(ax_str),'XMinorTick','on','YMinorTick','on','TickLength',[0.02, 1]); %
        set(eval(ax_str),'TickDir','out','clipping','off');
        % Now determine xticks and xticklabels and minorxticks
        xtickformat('%.1f');
        xticklabel = get(eval(ax_str),'xticklabel');
        xtick = get(eval(ax_str),'xtick');
        % remove current labels - and replace with own values closer to axis
        set(eval(ax_str),'xticklabel',[]);
        % now plot manually - for all except 12
        dy = a(3) - 0.075*(a(4)-a(3));
        for rr = 1:numel(xtick)
            text(eval(ax_str),xtick(rr),dy,xticklabel(rr),'fontsize',fontsize_axes,'HorizontalAlignment', 'center');
        end
        % Now use the same labels on y-axis
        set(eval(ax_str),'ytick',[0:0.5:2.5],'yticklabel',[0:0.5:2.5]); 
        ytickformat('%.1f');
        % set box property to off and remove background color
        set(eval(ax_str),'box','off','color','none')
        % create new, empty axes with box but without ticks
        h2 = axes('Position',get(eval(ax_str),'Position'),'box','on','xtick',[],'ytick',[]);
        % set original axes as active
        linkaxes([eval(ax_str) h2]);
        % set axis
       % axis(eval(ax_str),[x_min x_max y_min y_max ]);
        % add legend
        if zz == 1
            [lega,legb,legc,legd] = legend([plot1 plot2 plot3],{'${\rm \;\;Measured\;data}$',...
                '${\rm \;\;Infiltration\;form}$','${\rm \;\;residuals}$'},'interpreter','latex','box','off','location',...
                'northwest','fontsize',fontsize_labels);
        else
            [lega,legb,legc,legd] = legend([plot1 plot2 plot3],{'${\rm \;\;Measured\;data}$',...
                '${\rm \;\;Time\;form}$','${\rm \;\;residuals}$'},'interpreter','latex','box','off','location',...
                'northwest','fontsize',fontsize_labels);
        end
        set(legb(1),'Color','r'); set(legb(2),'Color',eval(char(color_leg(zz)))); set(legb(3),'Color',color_conf_line);
        set(legb(5),'linewidth',1.5,'markersize',8);
        set(legb(6),'linewidth',2.5);
        set(legb(7),'linewidth',1);
        set(gcf, 'Renderer', 'painters');
        %leg_pos = get(lega,'Position');
        %set(lega,'Position',[leg_pos(1) + 0.22 leg_pos(2:4)]);
        axis(eval(ax_str),[x_min x_max y_min y_max ]);
        set(gcf,'color','w');
        set(legb(1),'Color','r'); set(legb(2),'Color',eval(char(color_leg(zz)))); set(legb(3),'Color',color_conf_line);
        set(legb(5),'linewidth',1.5,'markersize',8);
        set(legb(6),'linewidth',2.5);
        set(legb(7),'linewidth',1);
        axis([0 1 0 2.6]);
        % now use export_fig function from Figure_resolution directory
        switch zz
            case 1
                evalstr = strcat('export_fig figure_residuals_I',{' '},'-','pdf');
            case 2
                evalstr = strcat('export_fig figure_residuals_t',{' '},'-','pdf');
        end
        % now save file
        if save_fig == 0, eval(char(evalstr)); end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 10: JACOBIAN ANALYTIC VERSUS NUMERICAL: TIME FORM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fig_10 == 1
    % Define default parameter values
    S = 2; Ks = 1; B = 1.5;
    eta = [ S Ks B ];
    I1 = [0.01:0.01:20]'; plugin.I = I1;
    % What is time length?
    def1 = Haverkamp_t(eta,plugin);
    % Compute Jacobian
    J_anal = jac_Haverkamp_anal(eta,plugin,2);
    I2 = [0:1:20]'; plugin.I = I2;
    J_num = jac_Haverkamp_num(eta,plugin,2);
    % Plot
    x_max = max(I1);
    %  y_max = [ 0.2 0.3 0.5 ];
    y_max = [ 0.2 1.5 0.6 ];
    y_min = [ -1.8 -18 -0.03 ];
    % Create figure name
    fig_name = 'Figure 10 of VZJ paper: Jacobian of infiltration equation: time - form';
    % Create figure
    figure('unit','inches','name',fig_name,'PaperOrientation', ...
        'portrait','position',[0.3 0.5 18.5 6.5],'numbertitle','off');
    % Define index
    index = {'(A)','(B)','(C)','(D)','(E)','(F)','(G)','(H)','(I)','(J)','(K)','(L)','(M)'};
    % Plot 3 figures
    for zz = 1:3
        ax_str = strcat('ax',num2str(zz));
        evalstr = strcat('ax',num2str(zz),{' '},'= axes(''units'',''inches'')'); eval(char(evalstr));
        evalstr = strcat('axpos',num2str(zz),{' '},'= [ 1.35 + (zz-1)*5.9 , 1.0 4.5 4.5 ]'); eval(char(evalstr));
        evalstr = strcat('set(','ax',num2str(zz),',','''position''',',axpos',num2str(zz),')'); eval(char(evalstr));
        switch zz
            case 1
                plot1 = plot(eval(ax_str),I1,J_anal(:,1),'k','linewidth',2);
                ylabel(eval(ax_str),'$\partial t/\partial S\;\; {\rm (h^{3/2}/cm)}$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [-0.18, 0.5, 0],'fontweight','normal');
            case 2
                plot(eval(ax_str),I1,J_anal(:,2),'k','linewidth',2);
                ylabel(eval(ax_str),'$\partial t/\partial K_{\rm s}\;\; {\rm (h^{2}/cm)}$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [-0.2, 0.5, 0],'fontweight','normal');
            case 3
                plot(eval(ax_str),I1,J_anal(:,3),'k','linewidth',2);
                ylabel(eval(ax_str),'$\partial t/\partial \beta\;\; {\rm (h)}$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [-0.16, 0.5, 0],'fontweight','normal');
        end
        hold on
        % Now plot the numerical Jacobian
        plot2 = plot(eval(ax_str),I2,J_num(:,zz),'ro','color',color_LMtime,'markerfacecolor','w','linewidth',2);
        axis([0 x_max y_min(zz) y_max(zz) ]);
        % Add (A), (B) and (C)
        set(gca,'fontsize',fontsize_axes);
        % Add labels on x and y-axis
        xlabel(eval(ax_str),'${\rm Cumulative\;infiltration}, I \; {\rm (\,{\rm cm}\,)}$','interpreter','latex',...
            'fontsize',fontsize_labels,'Units','normalized','Position',[0.5, -0.14, 0],'fontweight','normal');
        a = axis;
        set(eval(ax_str),'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 1]); %
        set(eval(ax_str),'TickDir','out','clipping','off');
        % Now determine xticks and xticklabels and minorxticks
        xticklabel = get(eval(ax_str),'xticklabel');
        xtick = get(eval(ax_str),'xtick');
        % remove current labels - and replace with own values closer to axis
        set(eval(ax_str),'xticklabel',[]);
        % now plot manually - for all except 12
        dy = a(3) - 0.07*(a(4)-a(3));
        for rr = 1:numel(xtick)
            text(eval(ax_str),xtick(rr),dy,xticklabel(rr),'fontsize',fontsize_axes,'HorizontalAlignment', 'center');
        end
        %        xtickformat('%.0f')
        ytickformat('%.1f')
        % Now create minor ticks - y axis
        %        eval(strcat(ax_str,'.XAxis.MinorTickValues = xminortick'));
        axis(eval(ax_str),[0 5 a(3) y_max(zz) ]);
        % set box property to off and remove background color
        set(eval(ax_str),'box','off','color','none')
        % create new, empty axes with box but without ticks
        h2 = axes('Position',get(eval(ax_str),'Position'),'box','on','xtick',[],'ytick',[]);
        % set original axes as active
        % axes(eval(ax_str))
        % link axes in case of zooming and make rounding box square
        linkaxes([eval(ax_str) h2]);
        % add legend
        if zz == 1
            [lega,legb,legc,legd] = legend([plot1 plot2],{'${\rm \;\;Analytic}$','${\rm \;\;Numerical}$'},'interpreter','latex','box','off','location',...
                'northeast','fontsize',fontsize_labels);
            set(legb(1),'Color','k'); set(legb(5),'linewidth',1.5,'markersize',8); set(legb(2),'Color',color_LMtime);
            set(legb(3),'linewidth',2.5);
            % set(eval(ax_str),'PaperPositionMode','auto');
            set(gcf, 'Renderer', 'painters')
        end
        axis([0 x_max y_min(zz) y_max(zz) ]);
        a = axis;
        % Add label (A), (B), or (C)
        text(eval(ax_str),a(1)+0.02*(a(2)-a(1)),a(3)+0.96*(a(4)-a(3)),char(index(zz)),...
            'fontsize',fontsize_labels+2,'interpreter','latex');
        set(gcf,'color','w');
        % add 2nd x-axis with time
        yloc = a(3)+1.02*(a(4)-a(3));
        line(eval(ax_str),[0 x_max],[yloc yloc],'color','k','linewidth',0.5);
        % now set vertical lines at right places
        mult = x_max / max(def1); ttick = [ 0 4 8 12 16 ]; tticklabel = {'0','4','8','12','16'};
        dy2 = yloc + 0.072*(a(4)-a(3));
        for rr = 1:numel(ttick)
            line(eval(ax_str),[ttick(rr) ttick(rr)]*mult,[yloc yloc+0.03*(a(4)-a(3))],'color','k','linewidth',0.5);
            % add label
            text(eval(ax_str),mult*ttick(rr),dy2,tticklabel(rr),'fontsize',fontsize_axes,'HorizontalAlignment', 'center');
        end
        % minor ticks
        minorttick = [0:18]; minorttick = minorttick(~ismember(minorttick,ttick)); %[ 1 2 3 5 6 7 9 10 11 13 14 15 17 18 ];
        for rr = 1:numel(minorttick)
            line(eval(ax_str),[minorttick(rr) minorttick(rr)]*mult,[yloc yloc+0.015*(a(4)-a(3))],'color','k','linewidth',0.5);
        end
        % add axes label
        dy3 = yloc + 0.15*(a(4)-a(3));
        text(eval(ax_str),10,dy3,'${\rm Time}, t{\rm \;(h)}$','fontsize',fontsize_labels,'HorizontalAlignment', 'center','interpreter','latex');
        
    end
    set(legb(1),'Color','k'); set(legb(5),'linewidth',1.5, ...
        'markersize',8); set(legb(2),'Color',color_LMtime);
    set(legb(3),'linewidth',2.5);
    % add a line at bottom with infiltration times
    %  set(gca,'clipping','off');
    
    % addpath to get nicest pdf from printing to pdf within matlab
    addpath('C:\Jasper Vrugt\Papers\Book\Chapter 1\MATLAB\Figure_resolution');
    % now use export_fig function from Figure_resolution directory
    evalstr = strcat('export_fig figure_jac_time',{' '},'-','pdf');
    % now save file
    if save_fig == 1, eval(char(evalstr)); end
end
