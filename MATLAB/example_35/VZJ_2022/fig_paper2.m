% Creates figures from paper
close all hidden; clc; clear all

file_LM_SWIG = 'SWIG_646_LM_test4_Jac';
load_file = strcat('load',{' '},file_LM_SWIG);
eval(char(load_file));

% TABLES
tab_3 = 0;  % DONE: Posterior moments of S, Ks and beta
tab_4 = 0;  % DONE: Correlation among bivariate posterior samples
tab_5 = 0;  % DONE: Least squares values of S, Ks and beta
tab_6 = 0;  % DONE: least squares value of c, time validity and t_char
% FIGURES
fig_0 = 0;  % DONE: Example of residuals of infiltration and time forms
fig_1a = 0;  % DONE: Jacobian matrices of time solution Haverkamp
fig_1b = 0;  % DONE: Jacobian matrices of infiltration solution Haverkamp
fig_2 = 0;  % DONE: RESIDUAL FUNCTION OF IMPLICIT SOLUTION
fig_3a = 0;  % NOT USED: RESIDUAL FUNCTION OF IMPLICIT SOLUTION 1 SOIL - MANY TIMES
fig_3b = 0;  % DONE: RESIDUAL FUNCTION OF IMPLICIT SOLUTION 12 ORIGINAL SOILS - MANY TIMES
fig_4 = 1;  % DONE: NEWTON'S METHOD APPLIED TO RESIDUAL FUNCTION (= IMPLICIT FORM)
fig_5 = 0;  % DONE: OBSERVED VERSUS SIMULATED DATA OF SWIG - all soils
fig_6 = 0;  % DONE: OBSERVED VERSUS SIMULATED DATA OF SWIG - all soils (SELECTION)
fig_7 = 0;  % DONE: SWIG LM VALUES OF S, Ks AND BETA FROM INFILTRATION AND TIME
fig_8 = 0;  % DONE: SWIG LM BETA VALUES WITH CONFIDENCE INTERVAL
fig_9 = 0;  % DONE: HISTOGRAM OF NUMERICAL ERROR IMPLICIT SOLUTION
fig_9b = 0; % DONE: NUMERICAL ERROR INVESTIGATED
fig_10 = 0; % DONE: HISTOGRAM OF MAXIMUM MEASUREMENT TIME OF SWIG SAMPLES
fig_11 = 0;  % DONE: NOT USED SWIG LM VALUES OF S, Ks AND BETA FROM INFILTRATION AND TIME
fig_12 = 0;  % DONE: NOT USED SWIG DREAM BETA VALUES WITH CONFIDENCE INTERVAL
fig_13 = 0;  % DONE: posterior distribution of DREAM: S, Ks and beta: Haverkamp (SELECTION)
fig_14 = 0;  % TBD: Bivariate scatter plots of DREAM: (S,Ks); (S,beta) + Ks,beta (SELECTION)
fig_15 = 0;  % DONE: Number of iterations Newton's method different I_0's
fig_16 = 0;  % DONE: NUMBER OF ITERATIONS NEWTON'S METHOD AS FUNCTION OF DELTA T
fig_17 = 0;  % DONE: COMPARISON OF UNCERTAINTIES DERIVED FROM DREAM AND LM METHOD

save_fig = 0; % save figures or not

% % addpath([pwd,'\Data'],[pwd,'\Colormaps'],[pwd,'\Supporting Codes'],...
% %     [pwd,'\Export Figures']);
addpath('C:\Jasper Vrugt\Papers\Book\Chapter 1\MATLAB\Figure_resolution');
addpath('C:\Jasper Vrugt\Papers\Book\Chapter 1\MATLAB\response_surface');
addpath('C:\Jasper Vrugt\Papers\Book\Chapter 1\MATLAB\Export figures');

load SWIG_1D
% Rename content of file
data_SWIG = SWIG_1D;
% Unpack the data
n_soil = size(data_SWIG,1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                            FIGURES                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 0: RESIDUAL EXAMPLE OF INFILTRATION AND TIME FORM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Selection from data comparison after the fact (LM_SWIG figures)
color_conf_line = [0.6 0.6 0.6];
color_scatter = [0 102 255]/255;
color_LMtime = [ 96 168 48]/255; % [113, 123, 40]/255;
fontsize_axes = 21; fontsize_labels = 21;
color_leg = {'color_scatter','color_LMtime'};
if fig_0 == 1
    plugin.t = [ 0 : 0.1 : 1 ]';
    x_min = 0; y_min = 0;
    x_max = max(plugin.t); y_max = 2.6;
    % Create figure name
    fig_name = 'Residuals of infiltration form of Haverkamp';
    eta = [ 2 1 1.5 0.1 ]'; 
    % Create figure (change zz = 1:1 and zz = 2:2 for plotting)
    for zz = 2:2
        figure('unit','inches','name',fig_name,'PaperOrientation','portrait','position',[0.3 0.5 8.5 6]);
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
        % addpath to get nicest pdf from printing to pdf within matlab
    %    addpath('C:\Jasper Vrugt\Papers\Book\Chapter 1\MATLAB\Figure_resolution');
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
% FIGURE 1a: JACOBIAN ANALYTIC VERSUS NUMERICAL: TIME FORM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

color_scatter = [0 102 255]/255;
color_LMtime = [ 96 168 48]/255; % [113, 123, 40]/255;
if fig_1a == 1
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
    fig_name = 'Jacobian of infiltration equation: time - form';
    % Create figure
    figure('unit','inches','name',fig_name,'PaperOrientation','portrait','position',[0.3 0.5 18.5 6.5]);
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
        % %         % add 2nd x-axis with time
        % %         yloc = a(3)-0.25*(a(4)-a(3));
        % %         line(eval(ax_str),[0 x_max],[yloc yloc],'color','k','linewidth',0.5);
        % %         % now set vertical lines at right places
        % %         mult = x_max / max(def1); ttick = [ 0 4 8 12 16 ]; tticklabel = {'0','4','8','12','16'};
        % %         dy2 = yloc - 0.06*(a(4)-a(3));
        % %         for rr = 1:numel(ttick)
        % %             line(eval(ax_str),[ttick(rr) ttick(rr)]*mult,[yloc yloc-0.03*(a(4)-a(3))],'color','k','linewidth',0.5);
        % %             % add label
        % %             text(eval(ax_str),mult*ttick(rr),dy2,tticklabel(rr),'fontsize',fontsize_axes,'HorizontalAlignment', 'center');
        % %         end
        % %         % minor ticks
        % %         minorttick = [0:18]; minorttick = minorttick(~ismember(minorttick,ttick)); %[ 1 2 3 5 6 7 9 10 11 13 14 15 17 18 ];
        % %         for rr = 1:numel(minorttick)
        % %             line(eval(ax_str),[minorttick(rr) minorttick(rr)]*mult,[yloc yloc-0.015*(a(4)-a(3))],'color','k','linewidth',0.5);
        % %         end
        % %         % add axes label
        % %         dy3 = yloc - 0.14*(a(4)-a(3));
        % %         text(eval(ax_str),10,dy3,'${\rm Time},\;t{\rm \;\;(\,h\,)}$','fontsize',fontsize_labels,'HorizontalAlignment', 'center','interpreter','latex');
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
    set(legb(1),'Color','k'); set(legb(5),'linewidth',1.5,'markersize',8); set(legb(2),'Color',color_LMtime);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 1b: JACOBIAN OF INFILTRATION FORM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if fig_1b == 1
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
    fig_name = 'Jacobian of infiltration equation - infiltration form';
    % Create figure
    figure('unit','inches','name',fig_name,'PaperOrientation','portrait','position',[0.3 0.5 12 14]);
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
        % %         switch zz
        % %             case 1
        % %                 ylabel(eval(ax_str),'$\partial I/\partial S\;\; {\rm (h^{1/2})}$','interpreter',...
        % %                     'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
        % %                     [-0.18, 0.5, 0],'fontweight','normal');
        % %             case 2
        % %                 ylabel(eval(ax_str),'$\partial I/\partial K_{\rm s}\;\; {\rm (h)}$','interpreter',...
        % %                     'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
        % %                     [-0.16, 0.5, 0],'fontweight','normal');
        % %             case 3
        % %                 ylabel(eval(ax_str),'$\partial I/\partial \beta\;\; {\rm (cm)}$','interpreter',...
        % %                     'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
        % %                     [-0.20, 0.5, 0],'fontweight','normal');
        % %             case 4
        % %                 ylabel(eval(ax_str),'$\partial I/\partial K_{\rm i}\;\; {\rm (h)}$','interpreter',...
        % %                     'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
        % %                     [-0.18, 0.5, 0],'fontweight','normal');
        % %         end
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
        % %         % add 2nd x-axis with time
        % %         yloc = a(3)-0.25*(a(4)-a(3));
        % %         line(eval(ax_str),[0 x_max],[yloc yloc],'color','k','linewidth',0.5);
        % %         % now set vertical lines at right places
        % %         mult = x_max / max(def1); ttick = [ 0 4 8 12 16 ]; tticklabel = {'0','4','8','12','16'};
        % %         dy2 = yloc - 0.06*(a(4)-a(3));
        % %         for rr = 1:numel(ttick)
        % %             line(eval(ax_str),[ttick(rr) ttick(rr)]*mult,[yloc yloc-0.03*(a(4)-a(3))],'color','k','linewidth',0.5);
        % %             % add label
        % %             text(eval(ax_str),mult*ttick(rr),dy2,tticklabel(rr),'fontsize',fontsize_axes,'HorizontalAlignment', 'center');
        % %         end
        % %         % minor ticks
        % %         minorttick = [0:18]; minorttick = minorttick(~ismember(minorttick,ttick)); %[ 1 2 3 5 6 7 9 10 11 13 14 15 17 18 ];
        % %         for rr = 1:numel(minorttick)
        % %             line(eval(ax_str),[minorttick(rr) minorttick(rr)]*mult,[yloc yloc-0.015*(a(4)-a(3))],'color','k','linewidth',0.5);
        % %         end
        % %         % add axes label
        % %         dy3 = yloc - 0.14*(a(4)-a(3));
        % %         text(eval(ax_str),10,dy3,'${\rm Time},\;t{\rm \;\;(\,h\,)}$','fontsize',fontsize_labels,'HorizontalAlignment', 'center','interpreter','latex');
        
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
% FIGURE 2: RESIDUAL FUNCTION OF IMPLICIT SOLUTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if fig_2 == 1
    % Define default parameter values
    S = 2; Ks = 1; B = 1.5; Ki = 0; dK = Ks - Ki; xi = dK/S^2;
    % Define residual function
    r = @(I,t) I-1/(2*xi)*log(exp(2*xi*B*(I-Ki*t))/B + ...
        (B-1)/B) - (dK*(1-B)+Ki)*t;                     % Residual function
    dr = @(I,t) 1-(B*exp(2*xi*B*(I-Ki*t))) / ...
        (exp(2*xi*B*(I-Ki*t))+B-1);                     % Derivative function
    % Define a time
    t = 3; % hours
    % Determine root
    fun = @(I) r(I,t); Iroot = fzero(fun,2);
    % Now plot
    I = [0.01:0.01:20]';
    y_max = 2;
    y_min = -8;
    x_max = max(I);
    % Create figure name
    fig_name = 'Residual function of implicit solution';
    % Create figure
    figure('unit','inches','name',fig_name,'PaperOrientation','portrait','position',[0.3 0.5 6.2 6.5]);
    % Define index
    index = {'(A)','(B)','(C)','(D)','(E)','(F)','(G)','(H)','(I)','(J)','(K)','(L)','(M)'};
    % Plot 3 figures
    zz = 1;
    ax_str = strcat('ax',num2str(zz));
    evalstr = strcat('ax',num2str(zz),{' '},'= axes(''units'',''inches'')'); eval(char(evalstr));
    evalstr = strcat('axpos',num2str(zz),{' '},'= [ 1.35 + (zz-1)*5.9 , 1.0 4.5 4.5 ]'); eval(char(evalstr));
    evalstr = strcat('set(','ax',num2str(zz),',','''position''',',axpos',num2str(zz),')'); eval(char(evalstr));
    plot1 = plot(eval(ax_str),I,r(I,t),'k','linewidth',2);
    ylabel(eval(ax_str),'$R(I,t)\;\; {\rm (cm)}$','interpreter',...
        'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
        [-0.18, 0.5, 0],'fontweight','normal');
    % Now plot the root
    hold on;
    plot2 = line(eval(ax_str),[0 Iroot],[0 0],'color','r','LineStyle','--','linewidth',1);
    plot3 = line(eval(ax_str),[Iroot Iroot],[-8.0 0],'color','r','LineStyle','--','linewidth',1);
    % Now plot root itself
    plot4 = plot(eval(ax_str),Iroot,0,'rx','linewidth',2,'markersize',12);
    axis([0 x_max y_min(zz) y_max(zz) ]);
    % Add (A), (B) and (C)
    set(gca,'fontsize',fontsize_axes);
    % Add labels on x and y-axis
    xlabel(eval(ax_str),'${\rm Cumulative\;infiltration}, I \; {\rm (\,{\rm cm}\,)}$','interpreter','latex',...
        'fontsize',fontsize_labels,'Units','normalized','Position',[0.5, -0.14, 0],'fontweight','normal');
    a = axis;
    set(eval(ax_str),'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 1]); %
    set(eval(ax_str),'TickDir','out');
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
    axis(eval(ax_str),[0 x_max y_min(zz) y_max(zz) ]);
    % set box property to off and remove background color
    set(eval(ax_str),'box','off','color','none')
    % create new, empty axes with box but without ticks
    h2 = axes('Position',get(eval(ax_str),'Position'),'box','on','xtick',[],'ytick',[]);
    % link axes in case of zooming and make rounding box square
    linkaxes([eval(ax_str) h2]);
    axis([0 x_max y_min(zz) y_max(zz) ]);
    a = axis;
    % Add legend
    [lega,legb,legc,legd] = legend(eval(ax_str),[plot1 plot4],{'$\;{\rm Residual\;function}, R(I,t)$','$\;{\rm Root}, I_{\rm r}$'},...
        'interpreter','latex','box','off','location','northeast','fontsize',fontsize_labels);
    %set(legb(1),'Color','k');
    %set(legb(3),'linewidth',2.5); set(legb(2),'Color','r'); set(legb(6),'linewidth',2.5,'markersize',12);
    % set(eval(ax_str),'PaperPositionMode','auto');
    set(gcf, 'Renderer', 'painters');
    set(gcf,'color','w');
    % repeat again
    %    set(gcf,'color','w');
    set(legb(3),'linewidth',2.5); set(legb(2),'Color','r'); set(legb(6),'linewidth',2.5,'markersize',12);
    % addpath to get nicest pdf from printing to pdf within matlab
    %addpath('C:\Jasper Vrugt\Papers\Book\Chapter 1\MATLAB\Figure_resolution');
    % now use export_fig function from Figure_resolution directory
    evalstr = strcat('export_fig figure_residual',{' '},'-','pdf');
    % now save file
    if save_fig == 1, eval(char(evalstr)); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 3a: RESIDUAL FUNCTION OF IMPLICIT SOLUTION (MANY CASES)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if fig_3a == 1
    % Define default parameter values
    N = 20;
    % Create figure name
    fig_name = 'Residual function of implicit solution: many cases';
    % Create figure
    figure('unit','inches','name',fig_name,'PaperOrientation','portrait','position',[0.3 0.5 6.2 6.5]);
    % Define index
    index = {'(A)','(B)','(C)','(D)','(E)','(F)','(G)','(H)','(I)','(J)','(K)','(L)','(M)'};
    zz = 1;
    ax_str = strcat('ax',num2str(zz));
    evalstr = strcat('ax',num2str(zz),{' '},'= axes(''units'',''inches'')'); eval(char(evalstr));
    evalstr = strcat('axpos',num2str(zz),{' '},'= [ 1.35 + (zz-1)*5.9 , 1.0 4.5 4.5 ]'); eval(char(evalstr));
    evalstr = strcat('set(','ax',num2str(zz),',','''position''',',axpos',num2str(zz),')'); eval(char(evalstr));
    % Determine root
    %fun = @(I) r(I,t); Iroot = fzero(fun,2);
    % Now plot
    I = [0.01:0.01:50]';
    y_max = 10;
    y_min = -6;
    x_max = max(I);
    plot2 = line(eval(ax_str),[0 x_max],[0 0],'color',[0.4 0.4 0.4],'LineStyle','-','linewidth',0.5); hold on;
    % Now loop over time
    for zzz = 1:N
        % Draw at random the parameters
        %par(zzz,1:4) = [ 20*rand 50*rand 0 2*rand ];
        %S = par(zzz,1); Ks = par(zzz,2); Ki = par(zzz,3); B = par(zzz,4);
        S = 2; Ks = 1; B = 1.5; Ki = 0;
        dK = Ks - Ki; xi = dK/S^2;
        % Define residual function
        r = @(I,t) I-1/(2*xi)*log(exp(2*xi*B*(I-Ki*t))/B + ...
            (B-1)/B) - (dK*(1-B)+Ki)*t;                     % Residual function
        % %         dr = @(I,t) 1-(B*exp(2*xi*B*(I-Ki*t))) / ...
        % %             (exp(2*xi*B*(I-Ki*t))+B-1);                     % Derivative function
        % Define a time
        t = zzz;
        % Now plot
        colors = (1 - 0.9*(zzz/N)) * ones(1,3);
        plot1 = plot(eval(ax_str),I,r(I,t),'color',colors,'linewidth',1);
        %        if zzz == 1
        %        end
    end
    ylabel(eval(ax_str),'$R(I,t)\;\; {\rm (cm)}$','interpreter',...
        'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
        [-0.18, 0.5, 0],'fontweight','normal');
    % Now plot the root
    hold on;
    %plot2 = line(eval(ax_str),[0 Iroot],[0 0],'color','r','LineStyle','--','linewidth',1);
    %plot3 = line(eval(ax_str),[Iroot Iroot],[-8.0 0],'color','r','LineStyle','--','linewidth',1);
    % Now plot root itself
    %plot4 = plot(eval(ax_str),Iroot,0,'rx','linewidth',2,'markersize',12);
    axis([0 x_max y_min(zz) y_max(zz) ]);
    % Add (A), (B) and (C)
    set(gca,'fontsize',fontsize_axes);
    % Add labels on x and y-axis
    xlabel(eval(ax_str),'${\rm Cumulative\;infiltration}, I \; {\rm (\,{\rm cm}\,)}$','interpreter','latex',...
        'fontsize',fontsize_labels,'Units','normalized','Position',[0.5, -0.14, 0],'fontweight','normal');
    a = axis;
    set(eval(ax_str),'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 1]); %
    set(eval(ax_str),'TickDir','out');
    % Now determine xticks and xticklabels and minorxticks
    xticklabel = get(eval(ax_str),'xticklabel');
    xtick = get(eval(ax_str),'xtick');
    % remove current labels - and replace with own values closer to axis
    set(eval(ax_str),'xticklabel',[]);
    % now plot manually - for all except 12
    dy = a(3)- 0.07*(a(4)-a(3));
    for rr = 1:numel(xtick)
        text(eval(ax_str),xtick(rr),dy,xticklabel(rr),'fontsize',fontsize_axes,'HorizontalAlignment', 'center');
    end
    %        xtickformat('%.0f')
    % %     switch soil_type
    % %         case {4,5,6,8}
    % %             ytickformat('%.1f')
    % %     end
    % Now create minor ticks - y axis
    axis(eval(ax_str),[0 x_max y_min(zz) y_max(zz) ]);
    % set box property to off and remove background color
    set(eval(ax_str),'box','off','color','none')
    % create new, empty axes with box but without ticks
    h2 = axes('Position',get(eval(ax_str),'Position'),'box','on','xtick',[],'ytick',[]);
    % link axes in case of zooming and make rounding box square
    linkaxes([eval(ax_str) h2]);
    axis([0 x_max y_min(zz) y_max(zz) ]);
    a = axis;
    % Add legend
    %    [lega,legb,legc,legd] = legend(eval(ax_str),[plot1 plot4],{'$\;{\rm Residual\;function}, R(I,t)$','$\;{\rm Root}, I_{\rm r}$'},...
    %        'interpreter','latex','box','off','location','northeast','fontsize',fontsize_labels);
    %set(legb(1),'Color','k');
    %set(legb(3),'linewidth',2.5); set(legb(2),'Color','r'); set(legb(6),'linewidth',2.5,'markersize',12);
    % set(eval(ax_str),'PaperPositionMode','auto');
    %    set(gcf, 'Renderer', 'painters');
    set(gcf,'color','w');
    % repeat again
    %    set(gcf,'color','w');
    %    set(legb(3),'linewidth',2.5); set(legb(2),'Color','r'); set(legb(6),'linewidth',2.5,'markersize',12);
    % addpath to get nicest pdf from printing to pdf within matlab
    %addpath('C:\Jasper Vrugt\Papers\Book\Chapter 1\MATLAB\Figure_resolution');
    % now use export_fig function from Figure_resolution directory
    evalstr = strcat('export_fig figure_residual_many',{' '},'-','pdf');
    % now save file
    if save_fig == 1, eval(char(evalstr)); end
end

fontsize_axes = 18; fontsize_labels = 18;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 3b: RESIDUAL FUNCTION OF ALL 12 SOILS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fig_3b == 1
    N = 10;
    CM = cbrewer('seq','Blues',N+1); CM = CM(2:N+1,1:3);
    % Own palette; min_B = [ 0 26 51 ]; max_B = [0 102 204];
    %     N = 20;
    % % #001a33	(0,26,51)
    % % #003366	(0,51,102)
    % % #004080	(0,64,128)
    % % #0059b3	(0,89,179)
    % % #0066cc	(0,102,204)
    %     min_B = [ 0 26 51 ]; max_B = [0 150 250];
    %     for i = 1:N
    %         CM(i,1:3) = min_B + i/N * ( max_B - min_B );
    %     end
    %    CM = CM/255;
    load opt_values_last % Optimum of all twelve soils
    load max_t_HYDR.mat  % Time it takes according to HYDRUS to inf 5 cm
    % Create figure name
    fig_name = 'Plot of residual function of all twelve soils';
    % Create figure
    figure('unit','inches','name',fig_name,'PaperOrientation','portrait','position',[0.1 0.5 19.7 13.9]);
    % Make a plot of cumulative infiltration of each soil
    index = {'(A) Clay','(B) Clay loam','(C) Loam','(D) Loamy sand','(E) Sand','(F) Sandy clay','(G) Sandy clay loam','(H) Sandy loam',...
        '(I) Silt','(J) Silt loam','(K) Silty clay','(L) Silty clay loam'};
    %max_I =   [ 56 75 280 4000 8000 35 350 1200 75 120 8   20 ];
    max_t = max_t_HYDR; %[ 15 9 3.5 0.24 0.10 23.9 3.5 0.75 11 7 140 59 ];
    %max_x = 240;
    % %     max_y = [  3.0  2.0  2.0  0.5  1.0  1.5  5.0  0.2  2.0  2.0  3.0  4.0 ];
    % %     min_y = [ -2.0 -1.0 -1.0 -0.8 -1.0 -1.0 -4.0 -0.2 -1.0 -1.0 -2.0 -3.0 ];
    max_y = [  2.4  1.5  1.8  0.5   1.0  1.3  5.0  0.22  2.0  2.0  2.4  4.0 ];
    min_y = [ -2.0 -1.0 -1.0 -0.65 -1.0 -1.0 -4.0 -0.2 -1.0 -1.0 -2.0 -3.0 ];
    
    %delta_I = [ 14 15  70 1000 2000 7  70  300  15 30  2   5 ];
    %minortick_I = [ 2 3 10 200 400  1  10  50   3  5  0.4  1 ];
    idx = [ 0 1 2 3 0 1 2 3 0 1 2 3 ];
    idy = [ 0 0 0 0 1 1 1 1 2 2 2 2 ];
    for soil_type = 1:12
        % dat = data{soil_type};
        % define max_x
        % max_y = max_I(soil_type);
        % define axis position
        ax_str = strcat('ax',num2str(soil_type));
        evalstr = strcat('ax',num2str(soil_type),{' '},'= axes(''units'',''inches'')'); eval(char(evalstr));
        if soil_type == 1
            evalstr = strcat('axpos',num2str(soil_type),{' '},'= [ 1.25 + idx(soil_type)*4.75 , 7.7-idy(soil_type)*3.4 3.0 2.65 ]'); eval(char(evalstr));
        else
            evalstr = strcat('axpos',num2str(soil_type),{' '},'= [ 1.25 + idx(soil_type)*4.75 , 7.7-idy(soil_type)*3.4 3.8 2.65 ]'); eval(char(evalstr));
        end
        evalstr = strcat('set(','ax',num2str(soil_type),',','''position''',',axpos',num2str(soil_type),')'); eval(char(evalstr));
        % Define range of infiltration values
        I = [0.01:0.01:6]';
        % define maximum of I
        x_max = max(I);
        % Plot y-axis (residual is zero)
        plot2 = line(eval(ax_str),[0 x_max],[0 0],'color',[0.4 0.4 0.4],'LineStyle','-','linewidth',0.5); hold on;
        % Define parameter values
        S = opt(soil_type,1); Ks = opt(soil_type,2); B = opt(soil_type,3); Ki = 0;
        dK = Ks - Ki; xi = dK/S^2;       % Now loop over time
        % Define residual function
        r = @(I,t) I-1/(2*xi)*log(exp(2*xi*B*(I-Ki*t))/B + ...
            (B-1)/B) - (dK*(1-B)+Ki)*t;                     % Residual function
        % Determine at what time we have infiltrated 5 cm of water
        Imax = 5; % hours
        % Determine root
        fun = @(t) r(Imax,t); t5 = fzero(fun,0.1); T(soil_type,1) = t5;
        % y_max and y_min
        y_max = max_y(soil_type); y_min = min_y(soil_type);
        % Loop over N measurement times
        for zzz = 1:N
            % Define a time
            %            max_t = [ 15 9 3.5 0.24 0.10 23.9 3.5 0.75 11 7 140 59 ]
            %t = max_t(soil_type)*(zzz/N);
            t = t5*(zzz/N);
            % Now plot
            % colors = (1 - 0.9*(zzz/N)) * ones(1,3);
            %            plot1 = plot(eval(ax_str),I,r(I,t),'color',colors,'linewidth',1);
            plot1 = plot(eval(ax_str),I,r(I,t),'color',CM(zzz,1:3),'linewidth',1);
            %        plot(eval(ax_str),dat(:,1),dat(:,2),'ro','markersize',5,'linewidth',1.5,'markerfacecolor',[1 1 1]); hold on
        end
        set(gca,'fontsize',fontsize_axes);
        % Add labels on x and y-axis
        if ismember(soil_type,9:12)
            tx = xlabel('${\rm Cumulative\;infiltration}, I \; {\rm (\,cm\,)}$','interpreter','latex',...
                'fontsize',fontsize_labels,'Units','normalized','Position',...
                [0.5, -0.3, 0],'fontweight','normal');
            pos_tx = get(tx,'position'); set(tx,'position',[ pos_tx(1) pos_tx(2)+0.08 pos_tx(3)]);
        end
        if ismember(soil_type,[1 5 9])
            ty = ylabel('$r(I,t,\Delta K,K_{\rm i},\beta,\xi) \; {\rm (cm)}$','interpreter',...
                'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                [-0.2, 0.5, 0],'fontweight','normal');
            pos_ty = get(ty,'position');
            if ismember(soil_type,9)
                set(ty,'position',[ pos_ty(1)+0.04 pos_ty(2:3)]);
            end
            % %             else
            % %                 set(ty,'position',[ pos_ty(1)+0.05 pos_ty(2:3)]);
            % %             end
        end
        text(eval(ax_str),0.02*x_max,0.95*y_max,char(index(soil_type)),'interpreter','latex','fontweight',...
            'normal','fontsize',fontsize_labels); %max_x = 240;
        % Add t_char - 2.5 percentile
        set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 1]); %
        set(gca,'TickDir','out');
        % define xticks
        % %         ytick = 0 : delta_I(soil_type) : max_I(soil_type);
        % %         yminortick = 0:minortick_I(soil_type):max_I(soil_type);
        % %         yminortick = yminortick(~ismember(yminortick,ytick));
        % %         ytickvalue = [];
        % %         for z = 1:numel(ytick)
        % %             if ( rem(delta_I(soil_type),1) == 0 )
        % %                 ytickvalue{z} = num2str(ytick(z),'%4.0f');
        % %             else
        % %                 ytickvalue{z} = num2str(ytick(z),'%4.2f');
        % %             end
        % %         end
        % %         set(eval(ax_str),'ytick',ytick,'yticklabel',ytickvalue);
        % %         % Now create minor ticks - x axis
        % %         eval(strcat(ax_str,'.YAxis.MinorTickValues = yminortick'));
        % define yticks
        xtick = 0 : 1 : x_max;
        xminortick = 0 : 0.2 : x_max;
        xminortick = xminortick(~ismember(xminortick,xtick));
        xtickvalue = [];
        for z = 1:numel(xtick), xtickvalue{z} = num2str(xtick(z),'%4.0f'); end
        set(eval(ax_str),'xtick',xtick,'xticklabel',xtickvalue);
        set(eval(ax_str),'xtick',xtick,'xticklabel',[]);
        axis(eval(ax_str),[0 x_max y_min y_max]);
        a = axis;
        % now plot manually
        dy = a(3)- 0.1*(a(4)-a(3));
        for rr = 1:numel(xtick)
            text(eval(ax_str),xtick(rr),dy,xtickvalue(rr),'fontsize',fontsize_axes,'HorizontalAlignment', 'center');
        end
        %        xtickformat('%.0f')
        switch soil_type
            case {2,4,5,6,8}
                ytickformat('%.1f')
        end
        if soil_type == 1
            cb = colorbar; colormap(CM);
            set(cb,'TickDirection','out','ticklength',0.02,'ticks',[0.15:0.2:1.0],'Ticklabels',[]);
            tickLabelsToUse = {'$0.2T$','0.4','0.6','0.8','1.0'}; %cb.TickLabels;
            %tickLabelsMod = cellfun(@(c) sprintf('%0.1f',str2double(c)) , tickLabelsToUse, 'uni', 0);
            %cb.TickLabels = tickLabelsMod;
            set(cb,'TickDirection','out','ticklength',0.02);
            cb_pos = get(cb,'Position');
            for r = 1:5
                switch r
                    case 1
                        evalstr_tick = strcat('$\frac{2}{10}T$');
                    case 2
                        evalstr_tick = strcat('$\frac{4}{10}T$');
                    case 3
                        evalstr_tick = strcat('$\frac{6}{10}T$');
                    case 4
                        evalstr_tick = strcat('$\frac{8}{10}T$');
                    case 5
                        evalstr_tick = strcat('$\,T$');
                end
                text(eval(ax_str),6.9,-1.37+(r-1)*0.875,evalstr_tick,'interpreter','latex',...
                    'HorizontalAlignment', 'left','fontsize',18);
            end
            set(cb,'Position',[0.2272    0.7339    0.0067    0.2512]); %cb_pos(1:2) 0.01 cb_pos(4)]);
        end
        % Now create minor ticks - y axis
        % %         eval(strcat(ax_str,'.XAxis.MinorTickValues = xminortick'));
        % set axis
        axis(eval(ax_str),[0 x_max y_min y_max]);
        % get handle to current axes (next lines remove xticks on top line of figure
        h1 = gca;
        % set box property to off and remove background color
        set(h1,'box','off','color','none')
        % create new, empty axes with box but without ticks
        h2 = axes('Position',get(h1,'Position'),'box','on','xtick',[],'ytick',[]);
        % set original axes as active
        axes(h1)
        % link axes in case of zooming and make rounding box square
        linkaxes([h1 h2]);
        axis(eval(ax_str),[0 x_max y_min y_max]);
    end
    set(gcf,'color','w');
    evalstr = strcat('export_fig figure_regUSDA',{' '},'-','pdf');
    % now save file
    if save_fig == 1, eval(char(evalstr)); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 4: NEWTON'S METHOD APPLIED TO RESIDUAL FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if fig_4 == 1
    % Define default parameter values
    S = 2; Ks = 1; B = 1.5; Ki = 0; dK = Ks - Ki; xi = dK/S^2;
    % Define residual function
    r = @(I,t) I-1/(2*xi)*log(exp(2*xi*B*(I-Ki*t))/B + ...
        (B-1)/B) - (dK*(1-B)+Ki)*t;                     % Residual function
    dr = @(I,t) 1-(B*exp(2*xi*B*(I-Ki*t))) / ...
        (exp(2*xi*B*(I-Ki*t))+B-1);                     % Derivative function
    % Define a time
    t = 1; % hours
    % Determine root
    fun = @(I) r(I,t); Iroot = fzero(fun,2);
    I = [];
    % Determine root - store results
    tolfun = 1e-12; maxiter = 20; I(1) = 0.45;
    for k = 2:maxiter
        dI = r(I(k-1),t)/dr(I(k-1),t);                  % Increment of I
        I(k) = I(k-1) - dI;                             % Update root guess
        if abs(dI) < tolfun                             % Check if converged
            break
        end
    end
    %fun = @(I) r(I,t); Iroot = fzero(fun,2);
    % Now plot
    I2 = [0.01:0.01:5.5]';
    y_max = 1;
    y_min = -1.5;
    x_max = max(I2);
    % Create figure name
    fig_name = 'Residual function of implicit solution';
    % Create figure
    figure('unit','inches','name',fig_name,'PaperOrientation','portrait','position',[0.3 0.5 7.2 6.5]);
    % Define index
    index = {'(A)','(B)','(C)','(D)','(E)','(F)','(G)','(H)','(I)','(J)','(K)','(L)','(M)'};
    % Plot 3 figures
    zz = 1;
    ax_str = strcat('ax',num2str(zz));
    evalstr = strcat('ax',num2str(zz),{' '},'= axes(''units'',''inches'')'); eval(char(evalstr));
    evalstr = strcat('axpos',num2str(zz),{' '},'= [ 1.35 + (zz-1)*5.9 , 1.0 4.5 4.5 ]'); eval(char(evalstr));
    evalstr = strcat('set(','ax',num2str(zz),',','''position''',',axpos',num2str(zz),')'); eval(char(evalstr));
    plot1 = plot(eval(ax_str),I2,r(I2,t),'k','linewidth',2); hold on
    plot2 = line(eval(ax_str),[0 x_max],[0 0],'color',[0.4 0.4 0.4],'LineStyle','-','linewidth',0.1);
    % Now plot root itself
    plot3 = plot(eval(ax_str),Iroot,0,'rx','linewidth',2,'markersize',12);
    % plot3 = line(eval(ax_str),[Iroot Iroot],[-8.0 0],'color','b','LineStyle','--','linewidth',1);
    ylabel(eval(ax_str),'$R(I,t)\;\; {\rm (cm)}$','interpreter',...
        'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
        [-0.18, 0.5, 0],'fontweight','normal');
    % Now plot the root
    hold on;
    % Now plot Newton progression of steps
    for i = 1:4
        % Now derivative derivative plot
        a = dr(I(i),t);
        % Function goes through (I(i),r(I(i,t))
        y = r(I(i),t);
        % Determine intercept, b
        b = y - a*(I(i) - 0);
        % Determine range of I-values for derivative plot
        Id = [ I(i) I(i+1) ];
        % Plot function, y = a*x + b;
        plot4 = plot(eval(ax_str),Id,a*Id+b,'b--','linewidth',1);
        % Plot small dot at intersection of zero
        plot5 = plot(eval(ax_str),Id(2),a*Id(2)+b,'bo','markerfacecolor','b','markersize',3,'linewidth',0.25);
        % Plot the solution
        plot6 = plot(eval(ax_str),I(i),r(I(i),t),'bo','markerfacecolor','w','linewidth',2,'markersize',5);
        % Add text
        evalstr = strcat('$I_{(',num2str(i-1),')}$');
        switch i
            case 3  % immediately before root solution
                text(eval(ax_str),I(i)-0.2,r(I(i),t) - 0.05*(y_max-y_min),char(evalstr),...
                    'interpreter','latex','fontsize',18,'color','b');
            case {1,2,4}
                text(eval(ax_str),I(i)-0.1,r(I(i),t) + 0.05*(y_max-y_min),char(evalstr),...
                    'interpreter','latex','fontsize',18,'color','b');
        end
    end
    axis([0 x_max y_min(zz) y_max(zz) ]);
    % Add (A), (B) and (C)
    set(gca,'fontsize',fontsize_axes);
    % Add labels on x and y-axis
    xlabel(eval(ax_str),'${\rm Cumulative\;infiltration}, I \; {\rm (\,{\rm cm}\,)}$','interpreter','latex',...
        'fontsize',fontsize_labels,'Units','normalized','Position',[0.5, -0.14, 0],'fontweight','normal');
    a = axis;
    set(eval(ax_str),'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 1]); %
    set(eval(ax_str),'TickDir','out');
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
    axis(eval(ax_str),[0 x_max y_min(zz) y_max(zz) ]);
    % set box property to off and remove background color
    set(eval(ax_str),'box','off','color','none')
    % create new, empty axes with box but without ticks
    h2 = axes('Position',get(eval(ax_str),'Position'),'box','on','xtick',[],'ytick',[]);
    % link axes in case of zooming and make rounding box square
    linkaxes([eval(ax_str) h2]);
    axis([0 x_max y_min(zz) y_max(zz) ]);
    a = axis;
    % Add label (A), (B), or (C)
    % Add legend
    [lega,legb,legc,legd] = legend(eval(ax_str),[plot1 plot3 plot4 plot6],...
        {'$\;{\rm Residual\;function}, R(I,t)$','$\;{\rm Root}, I_{\rm r}$',...
        '$\;{\rm Derivative}, R^{\prime}(I,t)$','$\;{\rm Iterates}, I_{(k)}$'},...
        'interpreter','latex','box','off','location','northeast','fontsize',fontsize_labels);
    %set(legb(1),'Color','k');
    % Change position of the legend
    pos_leg = get(lega,'Position');
    set(lega,'Position',[pos_leg(1)+0.2 pos_leg(2)+0.1 pos_leg(3) pos_leg(4) ]);
    set(legb(2),'Color','r'); set(legb(3:4),'Color','b');
    set(legb(5),'linewidth',2.5);
    set(legb(8),'linewidth',2.5,'markersize',12);
    set(legb(9),'linewidth',2.5);
    set(legb(12),'linewidth',2.5,'markersize',6);
    % now use export_fig function from Figure_resolution directory
    evalstr = strcat('export_fig figure_Newton',{' '},'-','pdf');
    % set(eval(ax_str),'PaperPositionMode','auto');
    set(gcf, 'Renderer', 'painters','color','w');
    %try
    %    set(eval(ax_str),'PaperPositionMode','auto');
    %catch
    % repeat again
    set(legb(2),'Color','r'); set(legb(3:4),'Color','b');
    set(legb(5),'linewidth',2.5);
    set(legb(8),'linewidth',2.5,'markersize',12);
    set(legb(9),'linewidth',2.5);
    set(legb(12),'linewidth',2.5,'markersize',6);
    %end
    % addpath to get nicest pdf from printing to pdf within matlab
    %    addpath('C:\Jasper Vrugt\Papers\Book\Chapter 1\MATLAB\Figure_resolution');
    % now save file
    if save_fig == 1, eval(char(evalstr)); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 5: PLOT MEASURED VERSUS FITTED CUMULATIVE INFILTRATION DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    for soil_type = 1:n_soil
        % Create figure name
        if rem(soil_type,rem_fig(fig_num)) == 0
            fig_name = char(strcat('Plot of observed and fitted cumulative infiltration:',{' '},num2str(fig_num)));
            % Create figure
            figure('unit','inches','name',fig_name,'PaperOrientation','portrait','position',[0.1 0.4 19.7 14]);
            % Set cf to 1;
            cf = 1; fig_num = fig_num + 1; fig_num = min(fig_num,ceil(n_soil/n_graphs));
        end
        dat = data_SWIG{soil_type};
        fx_LM = fx_opt_LM_SWIG{soil_type};
        fx_DREAM = fx_opt_DREAM_SWIG{soil_type};
        % define max_x
        x_max = max(dat(:,1)); y_max = max(dat(:,2));
        % define axis position
        ax_str = strcat('ax',num2str(soil_type));
        evalstr = strcat('ax',num2str(soil_type),{' '},'= axes(''units'',''inches'')'); eval(char(evalstr));
        evalstr = strcat('axpos',num2str(soil_type),{' '},'= [ 1.25 + idx(cf)*4.75 , 7.7-idy(cf)*3.4 3.8 2.65 ]'); eval(char(evalstr));
        evalstr = strcat('set(','ax',num2str(soil_type),',','''position''',',axpos',num2str(soil_type),')'); eval(char(evalstr));
        % plot
        plot(eval(ax_str),dat(:,1),dat(:,2),'ro','markersize',5,'linewidth',1.5,'markerfacecolor',[1 1 1]); hold on
        set(gca,'fontsize',fontsize_axes);
        % Now plot optimum of DREAM and LM
        for u = 1:4
            switch u
                case 1 % I optimized DREAM
                    plot(eval(ax_str),dat(:,1),fx_DREAM(:,1),'b','linewidth',1.5);
                case 2 % t optimized DREAM
                    plot(eval(ax_str),fx_DREAM(:,2),dat(:,2),'k','linewidth',1.5);
                case 3 % I optimized DREAM
                    plot(eval(ax_str),dat(:,1),fx_LM(:,1),'b--','linewidth',1.5);
                case 4 % t optimized DREAM
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
        % Add label with name of soil
        str = data_SWIG{soil_type,2};
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
        ytickformat('%.1f')
        
        %       text(eval(ax_str),0.02*x_max,0.95*y_max,char(index(soil_type)),'interpreter','latex','fontweight',...
        %           'normal','fontsize',fontsize_labels); %max_x = 240;
        % Add t_char - 2.5 percentile
        set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 1]); %
        set(gca,'TickDir','out');
        % define xticks
        %      delta_I = floor(y_max)/4;
        %      ytick = 0 : delta_I : y_max;
        %        ytick = 0 : delta_I(soil_type) : max_I(soil_type);
        %         yminortick = 0 : delta_I/8 : y_max;
        %         yminortick = yminortick(~ismember(yminortick,ytick));
        %         ytickvalue = [];
        %         for z = 1:numel(ytick)
        %             %if rem(delta_I == 0 )
        %                 ytickvalue{z} = num2str(ytick(z),'%4.0f');
        %             %else
        %             %    ytickvalue{z} = num2str(ytick(z),'%4.2f');
        %             %end
        %         end
        %xtickformat('%.0f')
        ytickformat('%.0f')
        %         set(eval(ax_str),'ytick',ytick,'yticklabel',ytickvalue);
        % Now create minor ticks - x axis
        %  eval(strcat(ax_str,'.YAxis.MinorTickValues = yminortick'));
        % define yticks
        % % %         delta_x = ceil(x_max)/4;
        % % %         xtick = 0 : 1 : x_max;
        % % %         xminortick = 0 : delta_x : x_max;
        % % %         xminortick = xminortick(~ismember(xminortick,xtick));
        % % %         xtickvalue = [];
        % % %         for z = 1:numel(xtick), xtickvalue{z} = num2str(xtick(z),'%4.0f'); end
        %        set(eval(ax_str),'xtick',xtick,'xticklabel',xtickvalue);
        %        set(eval(ax_str),'xtick',xtick,'xticklabel',[]);
        % now plot manually
        %         dy = -0.11 * y_max
        %         for rr = 1:numel(xtickvalue)
        %             text(xtick(rr),dy,char(xtickvalue{rr}),'fontsize',fontsize_axes,'HorizontalAlignment', 'center');
        %         end
        % Now create minor ticks - y axis
        %%%% eval(strcat(ax_str,'.XAxis.MinorTickValues = xminortick'));
        % set axis
        %axis([0 x_max 0 y_max]);
        % get handle to current axes (next lines remove xticks on top line of figure
        h1 = gca;
        % set box property to off and remove background color
        set(h1,'box','off','color','none')
        % create new, empty axes with box but without ticks
        h2 = axes('Position',get(h1,'Position'),'box','on','xtick',[],'ytick',[]);
        % set original axes as active
        axes(h1)
        % link axes in case of zooming and make rounding box square
        linkaxes([h1 h2]);
        % Update cf
        cf = cf + 1;
        % check before we move to new figure
        if rem(soil_type+1,rem_fig(fig_num)) == 0
            set(gcf,'color','w');
            % set(legb(1),'Color','r'); set(legb(4),'linewidth',1.5,'markersize',8); set(legb(5),'XData',[0.0603 0.1838],'linewidth',3); set(legb(2),'Color','b');
            evalstr = strcat('export_fig figure_5_',num2str(fig_num-1),{' '},'-','pdf');
            % now save file
            if save_fig == 1, eval(char(evalstr)); end
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 6: PLOT MEASURED VERSUS FITTED CUMULATIVE INFILTRATION DATA (SELECTION)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load SWIG_1D
% Rename content of file
%data_SWIG = SWIG_1D;
% Unpack the data
%n_soil = size(data_SWIG,1);
%load SWIG_646_LM_anal_Jac
% Now select
idx_sel_ORIG = [ 1 6 8 12 27 34 45 48 ]; % FROM ORIGINAL 50 SAMPLES
idx_selSWIG = [ 1 2 3 42 43 45 114 115 122 151 164:172 176 191 197 198 201 224 ...
    283 284 285 300 352 253 259 371 379 380 384 387 388 433 434 437 440 444 ...
    463 476 478 518 519 523 524];
idx_selSWIG = idx_selSWIG(idx_sel_ORIG);
% print table
for soil_type = 1:8
    % get soil_type
    str = data_SWIG{idx_selSWIG(soil_type),3};
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
    str_s = data_SWIG{idx_selSWIG(soil_type),2};
    pr1 = opt_LM_SWIG(idx_selSWIG(soil_type),1:3,1);
    pr2 = opt_LM_SWIG(idx_selSWIG(soil_type),1:3,2);
    for j = 1:3
        if pr1(j) < 0.001
            A = num2str(pr1(j),'%3.2e');
            ii = strfind(A,'e'); 
            R{j} = strcat(A(1:ii-1),'$\cdot 10^{-',A(end),'}$');
        else
            R{j} = num2str(pr1(j),'%3.2f');
        end
        if pr2(j) < 0.001
            A = num2str(pr2(j),'%3.2e');
            ii = strfind(A,'e'); 
            R2{j} = strcat(A(1:ii-1),'$\cdot 10^{-',A(end),'}$');
        else
            R2{j} = num2str(pr2(j),'%3.2f');
        end
    end
    fprintf('%10s & %d & %d & %s & %s & %s & & %s & %s & %s \\\\ \n',char(str_f),str_s,...
        Data(idx_selSWIG(soil_type),2),char(R{1}),char(R{2}),char(R{3}),char(R2{1}),char(R2{2}),char(R2{3}));
end
% 371 on our list;
% idx_selSWIG = idx_selSWIG(idx_sel_ORIG);
% 1) Sand
% 6) Sandy loam
% 8) Sandy clay loam
% 12) Loam
% 27) Silt loam
% 34) Clay loam
% 45) Silty clay loam
% 48) Clay
color_scatter = [0 102 255]/255;
color_LMtime = [ 96 168 48]/255; % [113, 123, 40]/255;
fontsize_axes = 21; fontsize_labels = 21;
if fig_6 == 1
    % Create figure name
    fig_name = 'Plot of observed and fitted data';
    index = {'(A)','(B)','(C)','(D)','(E)','(F)','(G)','(H)','(I)','(J)','(K)','(L)'};
    % Define maximum of
    idx = [ 0 1 2 3 0 1 2 3 ];
    idy = [ 0 0 0 0 1 1 1 1 ];
    n_graphs = 8;
    fig_num = 1; rem_fig = 1:n_graphs:n_soil;
    for soil_type = 1:8
        % Create figure name
        if rem(soil_type,rem_fig(fig_num)) == 0
            fig_name = char(strcat('Plot of observed and fitted cumulative infiltration:',{' '},num2str(fig_num)));
            % Create figure
            figure('unit','inches','name',fig_name,'PaperOrientation','portrait','position',[0.1 0.4 19.7 14]);
            % Set cf to 1;
            cf = 1; fig_num = fig_num + 1; fig_num = min(fig_num,ceil(n_soil/n_graphs));
        end
        dat = data_SWIG{idx_selSWIG(soil_type)};
        % Define plugin structure
        n_data = size(dat,1);
        % Store time and infiltration
        t_meas = dat(1:n_data,1); I_meas = dat(1:n_data,2);
        plugin.t = t_meas;
        % Evaluate Haverkamp
        % Next comment not needed if proper bounds are applied
        eta_I = opt_LM_SWIG(idx_selSWIG(soil_type),1:3,1); eta_I = max(eta_I,1e-4);
        I_pred = Haverkamp_I(eta_I,plugin);
        plugin.I = I_meas;
        % Evaluate Haverkamp
        % Next comment not needed if proper bounds are applied
        eta_t = opt_LM_SWIG(idx_selSWIG(soil_type),1:3,2); eta_t = max(eta_t,1e-4);
        t_pred = Haverkamp_t(eta_t,plugin);
        fx_LM = [ t_pred(:) I_pred(:) ];
        %        fx_LM = fx_opt_LM_SWIG{idx_sel(soil_type)};
        %        fx_DREAM = fx_opt_DREAM_SWIG{idx_sel(soil_type)};
        % define max_x
        x_max = max(dat(:,1)); y_max = max(dat(:,2));
        % define axis position
        ax_str = strcat('ax',num2str(idx_selSWIG(soil_type)));
        evalstr = strcat('ax',num2str(idx_selSWIG(soil_type)),{' '},'= axes(''units'',''inches'')'); eval(char(evalstr));
        evalstr = strcat('axpos',num2str(idx_selSWIG(soil_type)),{' '},'= [ 1.25 + idx(cf)*4.75 , 6-idy(cf)*4.8 3.8 4 ]'); eval(char(evalstr));
        evalstr = strcat('set(','ax',num2str(idx_selSWIG(soil_type)),',','''position''',',axpos',num2str(idx_selSWIG(soil_type)),')'); eval(char(evalstr));
        % plot
        plot1 = plot(eval(ax_str),dat(:,1),dat(:,2),'ro','markersize',5,'linewidth',1.5,'markerfacecolor',[1 1 1]); hold on
        set(gca,'fontsize',fontsize_axes);
        % Now plot optimum of DREAM and LM
        for u = 1:4
            switch u
                case 1 % I optimized LM
                    plot2 = plot(eval(ax_str),dat(:,1),fx_LM(:,2),'color',color_scatter,'linewidth',1.5);
                case 2 % t optimized LM
                    plot3 = plot(eval(ax_str),fx_LM(:,1),dat(:,2),'color',color_LMtime,'linewidth',1.5);
                    % %                 case 3 % I optimized DREAM
                    % %                     plot(eval(ax_str),dat(:,1),fx_DREAM(:,2),'b--','linewidth',1.5);
                    % %                 case 4 % t optimized DREAM
                    % %                     plot(eval(ax_str),fx_DREAM(:,1),dat(:,2),'k--','linewidth',1.5);
            end
        end
        % Add labels on x and y-axis
        if ismember(cf,5:8)
            tx = xlabel('${\rm Time}, t \; {\rm (\,h\,)}$','interpreter','latex',...
                'fontsize',fontsize_labels,'Units','normalized','Position',...
                [0.5, -0.23, 0],'fontweight','normal');
            pos_tx = get(tx,'position'); set(tx,'position',[ pos_tx(1) pos_tx(2)+0.08 pos_tx(3)]);
        end
        if ismember(cf,[1 5])
            ty = ylabel('$I \; {\rm (cm)}$','interpreter',...
                'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                [-0.2, 0.5, 0],'fontweight','normal');
            %             ty = ylabel('$\widetilde{I} \; {\rm (cm)}$','interpreter',...
            %                 'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
            %                 [-0.2, 0.5, 0],'fontweight','normal');
            % %             pos_ty = get(ty,'position');
            % %             if ismember(cf,1)
            % %                 set(ty,'position',[ pos_ty(1)+0.08 pos_ty(2:3)]);
            % %             else
            % %                 set(ty,'position',[ pos_ty(1)+0.005 pos_ty(2:3)]);
            % %             end
        end
        % Add label with name of soil
        str = data_SWIG{idx_selSWIG(soil_type),3};
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
        str_p = char(strcat(index(soil_type),{' '},str_f));
        text(eval(ax_str),a(1)+0.02*(a(2)-a(1)),a(3)+0.96*(a(4)-a(3)),str_p,'interpreter','latex','fontweight',...
            'bold','fontsize',fontsize_labels,'HorizontalAlignment', 'left');
        %        text(eval(ax_str),a(1)+0.02*(a(2)-a(1)),a(3)+1.01*(a(4)-a(3)),strcat(num2str(soil_type),':',num2str(cf),')'),'fontsize',fontsize_labels,'interpreter','latex');
        % Add sample number:
        str_s = data_SWIG{idx_selSWIG(soil_type),2};
        str_sample = char(strcat('${\rm Code}:',{' '},num2str(str_s),'$'));
        text(eval(ax_str),a(1)+0.96*(a(2)-a(1)),a(3)+0.05*(a(4)-a(3)),str_sample,'interpreter','latex','fontweight',...
            'bold','fontsize',fontsize_labels,'HorizontalAlignment', 'right');
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
        
        %       text(eval(ax_str),0.02*x_max,0.95*y_max,char(index(soil_type)),'interpreter','latex','fontweight',...
        %           'normal','fontsize',fontsize_labels); %max_x = 240;
        % Add t_char - 2.5 percentile
        set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 1]); %
        set(gca,'TickDir','out');
        % define xticks
        %      delta_I = floor(y_max)/4;
        %      ytick = 0 : delta_I : y_max;
        %        ytick = 0 : delta_I(soil_type) : max_I(soil_type);
        %         yminortick = 0 : delta_I/8 : y_max;
        %         yminortick = yminortick(~ismember(yminortick,ytick));
        %         ytickvalue = [];
        %         for z = 1:numel(ytick)
        %             %if rem(delta_I == 0 )
        %                 ytickvalue{z} = num2str(ytick(z),'%4.0f');
        %             %else
        %             %    ytickvalue{z} = num2str(ytick(z),'%4.2f');
        %             %end
        %         end
        %xtickformat('%.0f')
        ytickformat('%.0f')
        %axis([0 x_max 0 y_max]);
        % add legend
        if soil_type == 1
            [lega,legb,legc,legd] = legend([plot1 plot2 plot3],{'${\rm \;\;Measured\;data}$',...
                '${\rm \;\;Infiltration\;form}$','${\rm \;\;Time\;form}$'},...
                'interpreter','latex','box','off','location','southeast','fontsize',fontsize_labels);
            set(legb(1),'Color','r'); set(legb(2),'Color',color_scatter);
            set(legb(3),'Color',color_LMtime);
            set(legb(5),'linewidth',1.5,'markersize',6);
            set(legb(6),'linewidth',2.5);
            set(legb(8),'linewidth',2.5);
            leg_pos = get(lega,'Position');
            set(lega,'Position',[leg_pos(1)+0.01 leg_pos(2)+0.03 leg_pos(3:4)]);
            % set(eval(ax_str),'PaperPositionMode','auto');
            set(gcf, 'Renderer', 'painters')
        end
        % get handle to current axes (next lines remove xticks on top line of figure
        h1 = gca;
        % set box property to off and remove background color
        set(h1,'box','off','color','none')
        % create new, empty axes with box but without ticks
        h2 = axes('Position',get(h1,'Position'),'box','on','xtick',[],'ytick',[]);
        % set original axes as active
        axes(h1)
        % link axes in case of zooming and make rounding box square
        linkaxes([h1 h2]);
        % Update cf
        cf = cf + 1;
        % check before we move to new figure
        if rem(soil_type+1,rem_fig(fig_num)) == 0
            set(gcf,'color','w');
            % set(legb(1),'Color','r'); set(legb(4),'linewidth',1.5,'markersize',8); set(legb(5),'XData',[0.0603 0.1838],'linewidth',3); set(legb(2),'Color','b');
            evalstr = strcat('export_fig figure_data_fit',{' '},'-','pdf');
            % now save file
            if save_fig == 1, eval(char(evalstr)); end
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 7: SWIG LM VALUES OF S, Ks AND BETA FROM INFILTRATION AND TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
color_scatter = [0 102 255]/255;
% Selection from data comparison after the fact (LM_SWIG figures)
idx_notsel = [ 5 27 30 48 57 60 61 72 74 76 77 78 80 95 100 109 110 111 113 ...
    116 120 121 126 127 128 131 144 145 147 175 176 177 187 192 193 196  ...
    202 203 204 221 222 225 226 228 275 276 269 302 ...
    303 305 306 371 363 377 385 386 390 393 394 400 410 411 421 424 427 428 435 ...
    445 452 464 465 477 485 513 515 516 520 521 522 540 542 543 544 545 546 547 563 564 566 567 581 592 600 ...
    601:637 639 641 ];
idx_sel = 1:n_soil; idx_sel(idx_notsel) = [];
if fig_7 == 1
    % Load results from LM_SWIG
    %load SWIG_646_LM_anal_Jac.mat
    %load SWIG_646_LM_test_Jac
    eval(char(load_file));
    ii_data = Data(:,3) >= 1; idx_data = Data(ii_data,1);
    ii_not_data = 1:n_soil; ii_not_data(ii_data) = [];
    y_min = [  0  0    0 ];
    y_max = [ 200 250  2 ];
    dy = 0.04*(y_max - y_min);
    y_min = y_min - dy; y_max = y_max + dy;
    % Create figure name
    fig_name = 'Comparison of LM values of infiltration and time form Haverkamp';
    % Create figure
    %    figure('unit','inches','name',fig_name,'PaperOrientation','portrait','position',[0.3 0.5 18.5 6.5]);
    figure('unit','inches','name',fig_name,'PaperOrientation','portrait','position',[0.3 0.5 20 7.3]);
    % Define index
    index = {'(A)','(B)','(C)','(D)','(E)','(F)','(G)','(H)','(I)','(J)','(K)','(L)','(M)'};
    % Plot 3 figures
    for zz = 1:3
        ax_str = strcat('ax',num2str(zz));
        evalstr = strcat('ax',num2str(zz),{' '},'= axes(''units'',''inches'');'); eval(char(evalstr));
        %  evalstr = strcat('axpos',num2str(zz),{' '},'= [ 1.35 + (zz-1)*5.9 , 1.2 4.5 4.5 ]'); eval(char(evalstr));
        evalstr = strcat('axpos',num2str(zz),{' '},'= [ 1.45 + (zz-1)*6.5 , 1.7 5 5 ];'); eval(char(evalstr));
        evalstr = strcat('set(','ax',num2str(zz),',','''position''',',axpos',num2str(zz),');'); eval(char(evalstr));
        % Add one to one plot
        switch zz
            case 1
                %                 plot1 = plot(eval(ax_str),opt_LM_SWIG(ii_not_data,1,1),opt_LM_SWIG(ii_not_data,1,2),'s','color',color_scatter',...
                %                     'linewidth',1.5,'markersize',6,'markerfacecolor','w');
                plot1 = plot(eval(ax_str),opt_LM_SWIG(:,1,1),opt_LM_SWIG(:,1,2),'s','color',color_scatter',...
                    'linewidth',1.5,'markersize',6,'markerfacecolor','w'); hold on
                % %                 ylabel(eval(ax_str),'${\rm Soil\;sorptivity},\;S\;{\rm (cm/h^{1/2})}$','interpreter',...
                % %                     'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                % %                     [-0.16, 0.5, 0],'fontweight','normal');
                % %                 xlabel(eval(ax_str),'${\rm Soil\;sorptivity},\;S\;{\rm (cm/h^{1/2})}$','interpreter',...
                % %                     'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                % %                     [0.5, -0.13, 0],'fontweight','normal');
                ylabel(eval(ax_str),'$S_{\rm T}\;{\rm (cm/h^{1/2})}$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [-0.16, 0.5, 0],'fontweight','normal');
                xlabel(eval(ax_str),'$S_{\rm I}\;{\rm (cm/h^{1/2})}$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [0.5, -0.13, 0],'fontweight','normal');
            case 2
                plot1 = plot(eval(ax_str),opt_LM_SWIG(:,2,1),opt_LM_SWIG(:,2,2),'s','color',color_scatter',...
                    'linewidth',1.5,'markersize',6,'markerfacecolor','w'); hold on
                % %                 ylabel(eval(ax_str),'${\rm Saturated\;hydraulic\;conductivity},\;K_{\rm s}\;{\rm (cm/h)}$','interpreter',...
                % %                     'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                % %                     [-0.16, 0.5, 0],'fontweight','normal');
                % %                 xlabel(eval(ax_str),'${\rm Saturated\;hydraulic\;conductivity},\;K_{\rm s}\;{\rm (cm/h)}$','interpreter',...
                % %                     'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                % %                     [ 0.5, -0.13, 0],'fontweight','normal');
                ylabel(eval(ax_str),'$K_{\rm s,T}\;{\rm (cm/h)}$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [-0.16, 0.5, 0],'fontweight','normal');
                xlabel(eval(ax_str),'$K_{\rm s,I}\;{\rm (cm/h)}$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [ 0.5, -0.13, 0],'fontweight','normal');
            case 3
                plot1 = plot(eval(ax_str),opt_LM_SWIG(:,3,1),opt_LM_SWIG(:,3,2),'s','color',color_scatter',...
                    'linewidth',1.5,'markersize',6,'markerfacecolor','w'); hold on
                % %                 ylabel(eval(ax_str),'${\rm Coefficient},\;\beta\;(-)$','interpreter',...
                % %                     'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                % %                     [-0.14, 0.5, 0],'fontweight','normal');
                % %                 xlabel(eval(ax_str),'${\rm Coefficient},\;\beta\;(-)$','interpreter',...
                % %                     'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                % %                     [ 0.5, -0.13, 0],'fontweight','normal');
                ylabel(eval(ax_str),'$\beta_{\rm T}\;(-)$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [-0.14, 0.5, 0],'fontweight','normal');
                xlabel(eval(ax_str),'$\beta_{\rm I}\;(-)$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [ 0.5, -0.13, 0],'fontweight','normal');
        end
        plot2 = line(eval(ax_str),[0 y_max(zz)],[0 y_max(zz)],'color','k','linewidth',1);
        
        %         plot2 = plot(eval(ax_str),opt_LM_SWIG(idx_data,zz,1),opt_LM_SWIG(idx_data,zz,2),'s','color',color_scatter,...
        %             'linewidth',1.5,'markersize',6,'markerfacecolor',color_scatter);
        % plot uncertainty of coefficient beta
        axis(eval(ax_str),[y_min(zz) y_max(zz) y_min(zz) y_max(zz) ]);
        % Add (A), (B) and (C)
        set(gca,'fontsize',fontsize_axes);
        % Add labels on x and y-axis
        a = axis;
        set(eval(ax_str),'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 1]); %
        set(eval(ax_str),'TickDir','out','clipping','off');
        % Now determine xticks and xticklabels and minorxticks
        switch zz
            case 3
                xtickformat('%.1f');
            case {1,2}
                xtickformat('%.0f');
        end
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
        if zz == 3
            ytickformat('%.1f')
        end
        % Now create minor ticks - y axis
        %        eval(strcat(ax_str,'.XAxis.MinorTickValues = xminortick'));
        axis(eval(ax_str),[y_min(zz) y_max(zz) y_min(zz) y_max(zz) ]);
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
            [lega,legb,legc,legd] = legend([plot1 plot2],{'${\rm \;\;SWIG\;data}$','${\rm \;\;1:1\;Line}$'},...
                'interpreter','latex','box','off','location','southeast','fontsize',fontsize_labels);
            set(legb(1),'Color',color_scatter); set(legb(4),'linewidth',1.5,'markersize',8); set(legb(2),'Color','k');
            set(legb(5),'linewidth',2.5);
            set(gcf, 'Renderer', 'painters');
        end
        axis(eval(ax_str),[y_min(zz) y_max(zz) y_min(zz) y_max(zz) ]);
        a = axis;
        % Add label (A), (B), or (C)
        text(eval(ax_str),a(1)-0.04*(a(2)-a(1)),a(3)+1.05*(a(4)-a(3)),char(index(zz)),...
            'fontsize',fontsize_labels+2,'interpreter','latex');
        set(gcf,'color','w');
        % Infiltration and time form on axis
        if zz == 1
            text(eval(ax_str),20,-58.2,'${\rm Infiltration\;form}$','fontsize',fontsize_axes-3,'interpreter','latex',...
                'HorizontalAlignment', 'right');
            h = text(eval(ax_str),-62,-50,'${\rm Time\;form}$','fontsize',fontsize_axes-3,'interpreter','latex',...
                'HorizontalAlignment', 'left');
            set(h,'rotation',90);
            % One arrow from left to right with text over it
            x = [0.11 0.14];   % adjust length and location of arrow
            y = [0.0735 0.0735];
            annotation('textarrow',x,y,'FontSize',13,'Linewidth',0.5)
            % One arrow from left to right with text over it
            x = [0.011 0.011];   % adjust length and location of arrow
            y = [0.273 0.353];
            annotation('textarrow',x,y,'FontSize',13,'Linewidth',0.5)
            %            annotation('textbox',[.6 .3 .7 .27],'EdgeColor','none','String','Growth','FontSize',13,'Linewidth',2)
        end
        % inset
        if zz < 3
            ii = find(opt_LM_SWIG(:,zz,1) < 10);
            ax_str = strcat('ax',num2str(zz+5));
            evalstr = strcat('ax',num2str(zz+5),{' '},'= axes(''units'',''inches'');'); eval(char(evalstr));
            %  evalstr = strcat('axpos',num2str(zz),{' '},'= [ 1.35 + (zz-1)*5.9 , 1.2 4.5 4.5 ]'); eval(char(evalstr));
            evalstr = strcat('axpos',num2str(zz+5),{' '},'= [ 1.9 + (zz-1)*6.5 , 4.7 2 2 ];'); eval(char(evalstr));
            evalstr = strcat('set(','ax',num2str(zz+5),',','''position''',',axpos',num2str(zz+5),');'); eval(char(evalstr));
            plot6 = plot(eval(ax_str),opt_LM_SWIG(ii,zz,1),opt_LM_SWIG(ii,zz,2),'s','color',color_scatter,...
                'linewidth',0.5,'markersize',4,'markerfacecolor','w'); hold on
            plot2 = line(eval(ax_str),[0 10],[0 10],'color','k','linewidth',0.5);
            axis(eval(ax_str),[0 10 0 10]);
            set(gca,'fontsize',12);
            % Add labels on x and y-axis
            a = axis;
            xticklabel = [0:2:10];
            %            xtick = get(eval(ax_str),'xtick');
            xtick = [ 0 2 4 6 8 10 ];
            set(eval(ax_str),'xtick',xtick,'xticklabel',xticklabel);
            xticklabel = get(eval(ax_str),'xticklabel');
            set(eval(ax_str),'xtick',xtick,'xticklabel',[]);
            % remove current labels - and replace with own values closer to axis
            %            set(eval(ax_str),'xticklabel',[]);
            % now plot manually - for all except 12
            dy = a(3) - 0.1*(a(4)-a(3));
            for rr = 1:numel(xtick)
                text(eval(ax_str),xtick(rr),dy,xticklabel(rr,:),'fontsize',12,'HorizontalAlignment', 'center');
            end
            set(eval(ax_str),'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 1]); %
            set(eval(ax_str),'TickDir','out','clipping','on');
            set(eval(ax_str),'box','off','color','none')
            % create new, empty axes with box but without ticks
            h3 = axes('Position',get(eval(ax_str),'Position'),'box','on','xtick',[],'ytick',[]);
            linkaxes([eval(ax_str) h3]);
        end
        
    end
    set(legb(1),'Color',color_scatter); set(legb(4),'linewidth',1.5,'markersize',8); set(legb(2),'Color','k');
    % addpath to get nicest pdf from printing to pdf within matlab
    addpath('C:\Jasper Vrugt\Papers\Book\Chapter 1\MATLAB\Figure_resolution');
    % now use export_fig function from Figure_resolution directory
    evalstr = strcat('export_fig figure_LM_SWIG',{' '},'-','pdf');
    % now save file
    if save_fig == 1, eval(char(evalstr)); end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 8: SWIG BETA FROM INFILTRATION AND TIME + UNCERTAINTY 95
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
color_scatter = [0 102 255]/255;
% Selection from data comparison after the fact (LM_SWIG figures)
color_conf_line = [0.6 0.6 0.6];
fontsize_axes = 21; fontsize_labels = 21;
if fig_8 == 1
    % Load results from LM_SWIG
%    load SWIG_646_LM_anal_Jac.mat
%    load SWIG_646_LM_test_Jac
    eval(char(load_file));
    idx_data = [1:n_soil];
    y_min = [ 0 ];
    y_max = [ 2 ];
    dy = 0.04*(y_max - y_min);
    y_min = y_min - dy; y_max = y_max + dy;
    % Create figure name
    fig_name = 'Comparison of LM values of infiltration and time form Haverkamp';
    % Create figure
    figure('unit','inches','name',fig_name,'PaperOrientation','portrait','position',[0.3 0.5 12 10]);
    zz = 3
    ax_str = strcat('ax',num2str(zz));
    evalstr = strcat('ax',num2str(zz),{' '},'= axes(''units'',''inches'')'); eval(char(evalstr));
    evalstr = strcat('axpos',num2str(zz),{' '},'= [ 1.15 1.2 8 8 ]'); eval(char(evalstr));
    evalstr = strcat('set(','ax',num2str(zz),',','''position''',',axpos',num2str(zz),')'); eval(char(evalstr));
    % plot uncertainty of coefficient beta
    for i = 1:numel(idx_data)
        u = idx_data(i);
        % horizontal line
        plot3 = line(eval(ax_str),[min95(u,3,1) max95(u,3,1)],[opt_LM_SWIG(u,3,2) opt_LM_SWIG(u,3,2)],...
            'color',color_conf_line,'linewidth',0.1);
        % vertical line
        plot4 = line(eval(ax_str),[opt_LM_SWIG(u,3,1) opt_LM_SWIG(u,3,1)],[min95(u,3,2) max95(u,3,2)],...
            'color',color_conf_line,'linewidth',0.1);
        % plot standard deviation instead: horizontal line
        % Extract std_opt
        d = 3;
        
    end
    % Add one to one plot
    plot2 = line(eval(ax_str),[0 2],[0 2],'color','k','linewidth',1); hold on
    % Plot all data of SWIG
    plot1 = plot(eval(ax_str),opt_LM_SWIG(:,3,1),opt_LM_SWIG(:,3,2),'s','color',color_scatter',...
        'linewidth',1.5,'markersize',6,'markerfacecolor','w');
    % %     ylabel(eval(ax_str),'${\rm Coefficient},\;\beta\;(-)$','interpreter',...
    % %         'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
    % %         [-0.1, 0.5, 0],'fontweight','normal');
    % %     xlabel(eval(ax_str),'${\rm Coefficient},\;\beta\;(-)$','interpreter',...
    % %         'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
    % %         [ 0.5, -0.085, 0],'fontweight','normal');
    ylabel(eval(ax_str),'$\beta_{\rm T}\;(-)$','interpreter',...
        'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
        [-0.1, 0.5, 0],'fontweight','normal');
    xlabel(eval(ax_str),'$\beta_{\rm I}\;(-)$','interpreter',...
        'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
        [ 0.5, -0.085, 0],'fontweight','normal');
    axis(eval(ax_str),[y_min y_max y_min y_max ]);
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
    dy = a(3) - 0.045*(a(4)-a(3));
    for rr = 1:numel(xtick)
        text(eval(ax_str),xtick(rr),dy,xticklabel(rr),'fontsize',fontsize_axes,'HorizontalAlignment', 'center');
    end
    % Now use the same labels on y-axis
    set(eval(ax_str),'ytick',xtick,'yticklabel',xticklabel);
    %    ytickformat('%.1f')
    % Now create minor ticks - y axis
    %        eval(strcat(ax_str,'.XAxis.MinorTickValues = xminortick'));
    % set box property to off and remove background color
    set(eval(ax_str),'box','off','color','none')
    % create new, empty axes with box but without ticks
    h2 = axes('Position',get(eval(ax_str),'Position'),'box','on','xtick',[],'ytick',[]);
    % set original axes as active
    % axes(eval(ax_str))
    % link axes in case of zooming and make rounding box square
    linkaxes([eval(ax_str) h2]);
    % add legend
    [lega,legb,legc,legd] = legend([plot1 plot2 plot4],{'${\rm \;\;SWIG\;data}$',...
        '${\rm \;\;1:1\;Line}$','${\rm \;\;95\%\;intervals}$'},'interpreter','latex','box','off','location',...
        'northeast','fontsize',fontsize_labels);
    set(legb(1),'Color',color_scatter); set(legb(2),'Color','k'); set(legb(3),'Color',color_conf_line);
    set(legb(5),'linewidth',1.5,'markersize',8);
    set(legb(6),'linewidth',2.5);
    set(legb(7),'linewidth',1);
    set(gcf, 'Renderer', 'painters');
    leg_pos = get(lega,'Position');
    set(lega,'Position',[leg_pos(1) + 0.22 leg_pos(2:4)]);
    axis(eval(ax_str),[y_min y_max y_min y_max ]);
    a = axis;
    % Add label (A), (B), or (C)
    % text(eval(ax_str),a(1)+0.02*(a(2)-a(1)),a(3)+0.98*(a(4)-a(3)),char(index(zz)),...
    %     'fontsize',fontsize_labels+2,'interpreter','latex');
    set(gcf,'color','w');
    
    set(legb(1),'Color',color_scatter); set(legb(2),'Color','k'); set(legb(3),'Color',color_conf_line);
    set(legb(5),'linewidth',1.5,'markersize',8);
    set(legb(6),'linewidth',2.5);
    set(legb(7),'linewidth',1);
    
    % addpath to get nicest pdf from printing to pdf within matlab
    addpath('C:\Jasper Vrugt\Papers\Book\Chapter 1\MATLAB\Figure_resolution');
    % now use export_fig function from Figure_resolution directory
    evalstr = strcat('export_fig figure_LM_SWIG_beta_conf',{' '},'-','pdf');
    % now save file
    if save_fig == 1, eval(char(evalstr)); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 9: NUMERICAL ERROR OF IMPLICIT SOLUTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
color_scatter = [0 102 255]/255;
% Selection from data comparison after the fact (LM_SWIG figures)
color_conf_line = [0.6 0.6 0.6];
fontsize_axes = 21; fontsize_labels = 21;
if fig_9 == 1
    % Load results from numerical error analysis
    load numerr_results.mat
    % Create figure name
    fig_name = 'Histogram of numerical error of implicit solution';
    % Create figure
    figure('unit','inches','name',fig_name,'PaperOrientation','portrait','position',[0.3 0.5 12 8]);
    % Define index
    index = {'(A)','(B)','(C)','(D)','(E)','(F)','(G)','(H)','(I)','(J)','(K)','(L)','(M)'};
    % Plot 3 figures
    zz = 1;
    x_min = -12; x_max = -4; y_min = 0;
    ax_str = strcat('ax',num2str(zz));
    evalstr = strcat('ax',num2str(zz),{' '},'= axes(''units'',''inches'')'); eval(char(evalstr));
    evalstr = strcat('axpos',num2str(zz),{' '},'= [ 1.5 1.2 8 5 ]'); eval(char(evalstr));
    evalstr = strcat('set(','ax',num2str(zz),',','''position''',',axpos',num2str(zz),')'); eval(char(evalstr));
    % plot uncertainty of coefficient beta
    [N,x] = hist(log10(err_I(idx_I1)),20);
    y_max = 1.1*max(N)/sum(N);
    plot1 = bar(eval(ax_str),x,N/sum(N),'facecolor',color_conf_line,'edgecolor','k','linewidth',1); hold on
    xlabel(eval(ax_str),'${\rm Numerical\;error},\;\log_{10}({\rm err})\;\;{\rm (cm)}$','interpreter',...
        'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
        [ 0.5, -0.13, 0],'fontweight','normal');
    ylabel(eval(ax_str),'${\rm Normalized\;frequency}$','interpreter',...
        'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
        [ -0.13 , 0.5 , 0],'fontweight','normal');
    % Add text with how many not converged runs
    evalstr = strcat('${\rm Not\;converged}',{' '},'=',{' '},num2str(not_converged_I),'$');
    text(eval(ax_str),x_min+0.02*(x_max-x_min),y_max,char(evalstr),'interpreter','latex','fontsize',fontsize_axes);
    axis(eval(ax_str),[x_min x_max y_min y_max ]);
    % Add (A), (B) and (C)
    set(gca,'fontsize',fontsize_axes);
    % Add labels on x and y-axis
    a = axis;
    set(eval(ax_str),'XMinorTick','on','YMinorTick','on','TickLength',[0.02, 1]); %
    set(eval(ax_str),'TickDir','out','clipping','off');
    % Now determine xticks and xticklabels and minorxticks
    ytickformat('%.2f');
    xticklabel = get(eval(ax_str),'xticklabel');
    xtick = get(eval(ax_str),'xtick');
    % remove current labels - and replace with own values closer to axis
    set(eval(ax_str),'xticklabel',[]);
    % now plot manually - for all except 12
    dy = a(3) - 0.072*(a(4)-a(3));
    for rr = 1:numel(xtick)
        text(eval(ax_str),xtick(rr),dy,xticklabel(rr),'fontsize',fontsize_axes,'HorizontalAlignment', 'center');
    end
    % Now use the same labels on y-axis
    %set(eval(ax_str),'ytick',xtick,'yticklabel',xticklabel);
    %    ytickformat('%.1f')
    % Now create minor ticks - y axis
    %        eval(strcat(ax_str,'.XAxis.MinorTickValues = xminortick'));
    % set box property to off and remove background color
    set(eval(ax_str),'box','off','color','none')
    % create new, empty axes with box but without ticks
    h2 = axes('Position',get(eval(ax_str),'Position'),'box','on','xtick',[],'ytick',[]);
    % set original axes as active
    axes(eval(ax_str))
    % link axes in case of zooming and make rounding box square
    linkaxes([eval(ax_str) h2]);
    % add legend
    %     [lega,legb,legc,legd] = legend([plot1 plot2 plot3 plot4],{'${\rm \;\;SWIG\;data}$','$\;\;t_{n} \geq 1\;{\rm h}$',...
    %         '${\rm \;\;1:1\;Line}$','${\rm \;\;95\%\;intervals}$'},'interpreter','latex','box','off','location',...
    %         'northeast','fontsize',fontsize_labels);
    %     set(legb(1),'Color',color_scatter); set(legb(2),'Color',color_scatter); set(legb(4),'Color',color_conf_line);
    %     set(legb(6),'linewidth',2,'markersize',10);
    %     set(legb(8),'linewidth',2,'markersize',10);
    %     set(legb(9),'linewidth',3.5);
    %     set(legb(11),'linewidth',1);
    %     set(gcf, 'Renderer', 'painters');
    %     leg_pos = get(lega,'Position');
    %     set(lega,'Position',[leg_pos(1) + 0.22 leg_pos(2:4)]);
    axis(eval(ax_str),[x_min x_max y_min y_max ]);
    a = axis;
    % Add label (A), (B), or (C)
    % text(eval(ax_str),a(1)+0.02*(a(2)-a(1)),a(3)+0.98*(a(4)-a(3)),char(index(zz)),...
    %     'fontsize',fontsize_labels+2,'interpreter','latex');
    set(gcf,'color','w');
    
    %     set(legb(1),'Color',color_scatter); set(legb(2),'Color',color_scatter); set(legb(4),'Color',color_conf_line);
    %     set(legb(6),'linewidth',2,'markersize',10);
    %     set(legb(8),'linewidth',2,'markersize',10);
    %     set(legb(9),'linewidth',3.5);
    %     set(legb(11),'linewidth',1);
    % addpath to get nicest pdf from printing to pdf within matlab
    addpath('C:\Jasper Vrugt\Papers\Book\Chapter 1\MATLAB\Figure_resolution');
    % now use export_fig function from Figure_resolution directory
    evalstr = strcat('export_fig figure_hist_num_err',{' '},'-','pdf');
    % now save file
    if save_fig == 1, eval(char(evalstr)); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 9b: NUMERICAL ERROR INVESTIGATED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
red_color = [ 239 194 191 ]/255;
if fig_9b == 1
    % Load numerical error from numerr_results of program numerr.m
    load numerr_results.mat
    x_min2 = [ 0 0 0 0 ];
    x_max2 = [ 50 50 100 1.6 ];
    y_min2 = [ 0    0  0 0 ];
    y_max2 = [ 100  2  2 8.45 ];
    x_max = [ 52 52 104 1.6 ];
    x_min = [ -2 -2  -4 0 ];
    y_max = [ 104   2.1   2.1 8.45 ];
    y_min = [  -4  -0.1  -0.1 0 ];
    % idx and idy
    id_x = [ 0 1  0 1 ];
    id_y = [ 0 0  1 1 ];
    % Create figure name
    fig_name = 'Numerical error - as function of Haverkamp parameter values';
    % Create figure
    figure('unit','inches','name',fig_name,'PaperOrientation','portrait','position',[0.3 0.5 12 14]);
    % Define index
    index = {'(A)','(B)','(C)','(D)','(E)','(F)','(G)','(H)','(I)','(J)','(K)','(L)','(M)'};
    % Plot 3 figures
    for zz = 1:4
        ax_str = strcat('ax',num2str(zz));
        evalstr = strcat('ax',num2str(zz),{' '},'= axes(''units'',''inches'')'); eval(char(evalstr));
        evalstr = strcat('axpos',num2str(zz),{' '},'= [ 1.35 + id_x(zz)*5.3 , 6.4-id_y(zz)*5.3 3.8 3.8 ]'); eval(char(evalstr));
        evalstr = strcat('set(','ax',num2str(zz),',','''position''',',axpos',num2str(zz),')'); eval(char(evalstr));
        if zz < 4
            line([x_min2(zz) x_max2(zz)],[y_min2(zz) y_min2(zz)],'color','k','linewidth',0.25); hold on
            line([x_min2(zz) x_max2(zz)],[y_max2(zz) y_max2(zz)],'color','k','linewidth',0.25); hold on
            line([x_min2(zz) x_min2(zz)],[y_min2(zz) y_max2(zz)],'color','k','linewidth',0.25); hold on
            line([x_max2(zz) x_max2(zz)],[y_min2(zz) y_max2(zz)],'color','k','linewidth',0.25); hold on
        end
        switch zz
            case 1
                plot1 = plot(eval(ax_str),Eta(idx_I1,1),Eta(idx_I1,2),'rs','color',red_color,'markerfacecolor',red_color,...
                    'linewidth',0.5,'markersize',4);
                plot2 = plot(eval(ax_str),Eta(idx_I2,1),Eta(idx_I2,2),'bs','color',color_scatter,'markerfacecolor',color_scatter,...
                    'linewidth',2);
                xlabel(eval(ax_str),'$S\;\; {\rm (cm/h^{1/2})}$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [0.5, -0.14, 0],'fontweight','normal');
                ylabel(eval(ax_str),'$K_{\rm s}\;\; {\rm (cm/h)}$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [-0.18, 0.5, 0],'fontweight','normal');
                set(eval(ax_str),'xticklabel',[]);
                xtick = [0 10 20 30 40 50]; xticklabel = {'0','10','20','30','40','50'};
                minorxtick = [0:2:50]; minorxtick = minorxtick(~ismember(minorxtick,xtick));
            case 2
                plot(eval(ax_str),Eta(idx_I1,1),Eta(idx_I1,3),'rs','color',red_color,'markerfacecolor',red_color,...
                    'linewidth',0.5,'markersize',4);
                plot(eval(ax_str),Eta(idx_I2,1),Eta(idx_I2,3),'bs','color',color_scatter,'markerfacecolor',color_scatter,...
                    'linewidth',2);
                xlabel(eval(ax_str),'$S\;\; {\rm (cm/h^{1/2})}$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [0.5, -0.14, 0],'fontweight','normal');
                ylabel(eval(ax_str),'$\beta\;\; {\rm (-)}$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [-0.20, 0.5, 0],'fontweight','normal');
                set(eval(ax_str),'xticklabel',[]);
                xtick = [0 10 20 30 40 50]; xticklabel = {'0','10','20','30','40','50'};
                minorxtick = [0:2:50]; minorxtick = minorxtick(~ismember(minorxtick,xtick));
            case 3
                plot(eval(ax_str),Eta(idx_I1,2),Eta(idx_I1,3),'rs','color',red_color,'markerfacecolor',red_color,...
                    'linewidth',0.5,'markersize',4);
                plot(eval(ax_str),Eta(idx_I2,2),Eta(idx_I2,3),'bs','color',color_scatter,'markerfacecolor',color_scatter,...
                    'linewidth',2);
                xlabel(eval(ax_str),'$K_{\rm s}\;\; {\rm (cm/h)}$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [ 0.5,-0.14,0],'fontweight','normal');
                ylabel(eval(ax_str),'$\beta\;\; {\rm (-)}$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [-0.20, 0.5, 0],'fontweight','normal');
                set(eval(ax_str),'xticklabel',[]);
                xtick = [0 25 50 75 100]; xticklabel = {'0','25','50','75','100'};
                minorxtick = [0:5:100]; minorxtick = minorxtick(~ismember(minorxtick,xtick));
            case 4
                plugin.I = [0.01:0.01:10]';
                for z = 1:numel(idx_I2)
                    t = Haverkamp_t(Eta(idx_I2(z),1:3),plugin);
                    plot(eval(ax_str),t,plugin.I,'color',color_scatter,'linewidth',2); if z == 1; hold on; end
                end
                xlabel(eval(ax_str),'${\rm Time},\;t\;\; {\rm (h)}$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [ 0.5,-0.14,0],'fontweight','normal');
                ylabel(eval(ax_str),'$I\;(\,{\rm cm}\,)$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [-0.18, 0.5, 0],'fontweight','normal');
                xtick = [0 0.4 0.8 1.2 1.6]; xticklabel = {'0.0','0.4','0.8','1.2','1.6'};
                minorxtick = 0:0.1:1.6; minorxtick = minorxtick(~ismember(minorxtick,xtick));
                set(eval(ax_str),'xticklabel',[]);
        end
        hold on
        axis(eval(ax_str),[x_min(zz) x_max(zz) y_min(zz) y_max(zz) ]);
        % Add (A), (B) and (C)
        set(gca,'fontsize',fontsize_axes);
        % Add labels on x and y-axis
        a = axis;
        set(eval(ax_str),'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 1]); %
        set(eval(ax_str),'TickDir','out','clipping','off');
        set(eval(ax_str),'xtick',xtick,'xticklabel',[]);
        % Minor x and yticks
        eval(strcat(ax_str,'.XAxis.MinorTickValues = minorxtick'));
        % Now determine xticks and xticklabels and minorxticks
        % Set ticks
        %        xticklabel = get(eval(ax_str),'xticklabel');
        %        xtick = get(eval(ax_str),'xtick');
        % remove current labels - and replace with own values closer to axis
        %        set(eval(ax_str),'xticklabel',[]);
        % now plot manually - for all except 12
        dy = a(3) - 0.08*(a(4)-a(3));
        %xticklabel = xtick;
        for rr = 1:numel(xtick)
            text(eval(ax_str),xtick(rr),dy,xticklabel(rr),'fontsize',fontsize_axes,'HorizontalAlignment', 'center');
        end
        %        xtickformat('%.0f')
        switch zz
            case {2,3}
                ytickformat('%.1f');
        end
        % Now create minor ticks - y axis
        %        eval(strcat(ax_str,'.XAxis.MinorTickValues = xminortick'));
        %        axis(eval(ax_str),[x_min(zz) x_max(zz) a(3) y_max(zz) ]);
        % set box property to off and remove background color
        set(eval(ax_str),'box','off','color','none')
        % create new, empty axes with box but without ticks
        h2 = axes('Position',get(eval(ax_str),'Position'),'box','on','xtick',[],'ytick',[]);
        % set original axes as active
        axes(eval(ax_str))
        % link axes in case of zooming and make rounding box square
        linkaxes([eval(ax_str) h2]);
        % add legend
        if zz == 1
            [lega,legb,legc,legd] = legend([plot2],{'${\rm Overflow}$'},'interpreter','latex','box','off','location',...
                'northeast','fontsize',fontsize_labels);
            %  set(legb(1),'Color','r'); set(legb(3),'linewidth',1.5,'markersize',8);
            set(legb(1),'Color',color_scatter);
            set(legb(3),'linewidth',1.5,'markersize',8);
            leg_pos = get(lega,'Position');
            set(lega,'position',[leg_pos(1)-0.02 leg_pos(2)-0.02 leg_pos(3:4)]);
            % set(eval(ax_str),'PaperPositionMode','auto');
            set(gcf, 'Renderer', 'painters')
        end
        axis(eval(ax_str),[x_min(zz) x_max(zz) y_min(zz) y_max(zz) ]);
        a = axis;
        %         plot(ax(2)*[1,1],ax(3:4),'color',[0.1500 0.1500 0.1500],'linewidth',0.25)
        %         plot(ax(1:2),ax(4)*[1,1],'color',[0.1500 0.1500 0.1500],'linewidth',0.25)
        % Add label (A), (B), or (C)
        text(eval(ax_str),a(1)+0.03*(a(2)-a(1)),a(3)+1.02*(a(4)-a(3)),char(index(zz)),...
            'fontsize',fontsize_labels+2,'interpreter','latex');
        set(gcf,'color','w');
        
    end
    %set(legb(1),'Color',red_color); set(legb(3),'linewidth',1.5,'markersize',8);
    set(legb(1),'Color',color_scatter);
    set(legb(3),'linewidth',1.5,'markersize',8);
    %     % addpath to get nicest pdf from printing to pdf within matlab
    addpath('C:\Jasper Vrugt\Papers\Book\Chapter 1\MATLAB\Figure_resolution');
    % now use export_fig function from Figure_resolution directory
    evalstr = strcat('export_fig figure_num_err',{' '},'-','pdf');
    % now save file
    if save_fig == 1, eval(char(evalstr)); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 10: HISTOGRAM OF MEASUREMENT TIMES OF SWIG DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
color_scatter = [0 102 255]/255;
% Selection from data comparison after the fact (LM_SWIG figures)
fontsize_axes = 21; fontsize_labels = 21;
if fig_10 == 1
    
    % Create figure name
    fig_name = 'Histogram of numerical error of implicit solution';
    % Create figure
    figure('unit','inches','name',fig_name,'PaperOrientation','portrait','position',[0.3 0.5 12 8]);
    % Define index
    index = {'(A)','(B)','(C)','(D)','(E)','(F)','(G)','(H)','(I)','(J)','(K)','(L)','(M)'};
    % Plot 3 figures
    zz = 1;
    x_min = 0; x_max = 8; y_min = 0; %y_max = 0.25;
    ax_str = strcat('ax',num2str(zz));
    evalstr = strcat('ax',num2str(zz),{' '},'= axes(''units'',''inches'')'); eval(char(evalstr));
    evalstr = strcat('axpos',num2str(zz),{' '},'= [ 1.5 1.2 8 5 ]'); eval(char(evalstr));
    evalstr = strcat('set(','ax',num2str(zz),',','''position''',',axpos',num2str(zz),')'); eval(char(evalstr));
    % plot uncertainty of coefficient beta
    [N,x] = hist(Data(:,3),15); y_max = 1.15*max(N);
    plot1 = bar(eval(ax_str),x,N,'facecolor',color_conf_line,'edgecolor','k','linewidth',1); hold on
    % ,...,'color',color_conf_line,'linewidth',0.1);
    xlabel(eval(ax_str),'${\rm Duration\;of\;experiment},\;\tilde{t}_{n}\;{\rm (hour)}$','interpreter',...
        'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
        [ 0.5, -0.13, 0],'fontweight','normal');
    ylabel(eval(ax_str),'${\rm Number\;of\;SWIG\;experiments}$','interpreter',...
        'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
        [ -0.1 , 0.5 , 0],'fontweight','normal');
    % Add text with how many not converged runs
    %  evalstr = strcat('${\rm Not\;converged}',{' '},'=',{' '},num2str(not_converged),'$');
    %  text(eval(ax_str),-12.8,0.34,char(evalstr),'interpreter','latex','fontsize',fontsize_axes);
    axis(eval(ax_str),[x_min x_max y_min y_max ]);
    % Add (A), (B) and (C)
    set(gca,'fontsize',fontsize_axes);
    % Add labels on x and y-axis
    a = axis;
    set(eval(ax_str),'XMinorTick','on','YMinorTick','on','TickLength',[0.02, 1]); %
    set(eval(ax_str),'TickDir','out','clipping','off');
    % Now determine xticks and xticklabels and minorxticks
    ytickformat('%.0f');
    xticklabel = get(eval(ax_str),'xticklabel');
    xtick = get(eval(ax_str),'xtick');
    % remove current labels - and replace with own values closer to axis
    set(eval(ax_str),'xticklabel',[]);
    % now plot manually - for all except 12
    dy = a(3) - 0.072*(a(4)-a(3));
    for rr = 1:numel(xtick)
        text(eval(ax_str),xtick(rr),dy,xticklabel(rr),'fontsize',fontsize_axes,'HorizontalAlignment', 'center');
    end
    % Now use the same labels on y-axis
    %set(eval(ax_str),'ytick',xtick,'yticklabel',xticklabel);
    %    ytickformat('%.1f')
    % Now create minor ticks - y axis
    %        eval(strcat(ax_str,'.XAxis.MinorTickValues = xminortick'));
    % set box property to off and remove background color
    set(eval(ax_str),'box','off','color','none')
    % create new, empty axes with box but without ticks
    h2 = axes('Position',get(eval(ax_str),'Position'),'box','on','xtick',[],'ytick',[]);
    % set original axes as active
    axes(eval(ax_str))
    % link axes in case of zooming and make rounding box square
    linkaxes([eval(ax_str) h2]);
    % add legend
    %     [lega,legb,legc,legd] = legend([plot1 plot2 plot3 plot4],{'${\rm \;\;SWIG\;data}$','$\;\;t_{n} \geq 1\;{\rm h}$',...
    %         '${\rm \;\;1:1\;Line}$','${\rm \;\;95\%\;intervals}$'},'interpreter','latex','box','off','location',...
    %         'northeast','fontsize',fontsize_labels);
    %     set(legb(1),'Color',color_scatter); set(legb(2),'Color',color_scatter); set(legb(4),'Color',color_conf_line);
    %     set(legb(6),'linewidth',2,'markersize',10);
    %     set(legb(8),'linewidth',2,'markersize',10);
    %     set(legb(9),'linewidth',3.5);
    %     set(legb(11),'linewidth',1);
    %     set(gcf, 'Renderer', 'painters');
    %     leg_pos = get(lega,'Position');
    %     set(lega,'Position',[leg_pos(1) + 0.22 leg_pos(2:4)]);
    axis(eval(ax_str),[x_min x_max y_min y_max ]);
    a = axis;
    % Add label (A), (B), or (C)
    % text(eval(ax_str),a(1)+0.02*(a(2)-a(1)),a(3)+0.98*(a(4)-a(3)),char(index(zz)),...
    %     'fontsize',fontsize_labels+2,'interpreter','latex');
    set(gcf,'color','w');
    
    %     set(legb(1),'Color',color_scatter); set(legb(2),'Color',color_scatter); set(legb(4),'Color',color_conf_line);
    %     set(legb(6),'linewidth',2,'markersize',10);
    %     set(legb(8),'linewidth',2,'markersize',10);
    %     set(legb(9),'linewidth',3.5);
    %     set(legb(11),'linewidth',1);
    % addpath to get nicest pdf from printing to pdf within matlab
    addpath('C:\Jasper Vrugt\Papers\Book\Chapter 1\MATLAB\Figure_resolution');
    % now use export_fig function from Figure_resolution directory
    evalstr = strcat('export_fig figure_hist_meas_time',{' '},'-','pdf');
    % now save file
    if save_fig == 1, eval(char(evalstr)); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 11: SWIG DREAM VALUES OF S, Ks AND BETA FROM INFILTRATION AND TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fontsize_axes = 18; fontsize_labels = 18;
color_scatter = [0 102 255]/255;
if fig_11 == 1
    % Load results from LM_SWIG
%    load SWIG_646_LM_anal_Jac.mat
%    load SWIG_646_LM_test_Jac
    eval(char(load_file));
   % load SWIG_646_LM_anal_Jac.mat
    load DREAM_SWIG.mat
    %ii_data = Data(:,2)>=10 & Data(:,3)>= 2; idx_data = Data(ii_data,1);
    ii_data = Data(:,3) >= 1; idx_data = Data(ii_data,1);
    ii_not_data = 1:n_soil; ii_not_data(ii_data) = [];
    y_min = [  -1 -1 -0.08 ];
    y_max = [ 51  51  2.08 ];
    % Create figure name
    fig_name = 'Comparison of DREAM values of infiltration and time form Haverkamp';
    % Create figure
    %    figure('unit','inches','name',fig_name,'PaperOrientation','portrait','position',[0.3 0.5 18.5 6.5]);
    figure('unit','inches','name',fig_name,'PaperOrientation','portrait','position',[0.3 0.5 20 6.5]);
    % Define index
    index = {'(A)','(B)','(C)','(D)','(E)','(F)','(G)','(H)','(I)','(J)','(K)','(L)','(M)'};
    % Plot 3 figures
    for zz = 1:3
        ax_str = strcat('ax',num2str(zz));
        evalstr = strcat('ax',num2str(zz),{' '},'= axes(''units'',''inches'')'); eval(char(evalstr));
        %  evalstr = strcat('axpos',num2str(zz),{' '},'= [ 1.35 + (zz-1)*5.9 , 1.2 4.5 4.5 ]'); eval(char(evalstr));
        evalstr = strcat('axpos',num2str(zz),{' '},'= [ 1.15 + (zz-1)*6.4 , 1.2 5 5 ]'); eval(char(evalstr));
        evalstr = strcat('set(','ax',num2str(zz),',','''position''',',axpos',num2str(zz),')'); eval(char(evalstr));
        % Add one to one plot
        switch zz
            case 1
                %                 plot1 = plot(eval(ax_str),opt_LM_SWIG(ii_not_data,1,1),opt_LM_SWIG(ii_not_data,1,2),'s','color',color_scatter',...
                %                     'linewidth',1.5,'markersize',6,'markerfacecolor','w');
                plot1 = plot(eval(ax_str),opt_DREAM_SWIG(:,1,1),opt_DREAM_SWIG(:,1,2),'s','color',color_scatter',...
                    'linewidth',1.5,'markersize',6,'markerfacecolor','w'); hold on
                ylabel(eval(ax_str),'${\rm Soil\;sorptivity},\;S\;{\rm (cm/h^{1/2})}$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [-0.13, 0.5, 0],'fontweight','normal');
                xlabel(eval(ax_str),'${\rm Soil\;sorptivity},\;S\;{\rm (cm/h^{1/2})}$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [0.5, -0.13, 0],'fontweight','normal');
            case 2
                plot1 = plot(eval(ax_str),opt_DREAM_SWIG(:,2,1),opt_DREAM_SWIG(:,2,2),'s','color',color_scatter',...
                    'linewidth',1.5,'markersize',6,'markerfacecolor','w'); hold on
                ylabel(eval(ax_str),'${\rm Saturated\;hydraulic\;conductivity},\;K_{\rm s}\;{\rm (cm/h)}$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [-0.13, 0.5, 0],'fontweight','normal');
                xlabel(eval(ax_str),'${\rm Saturated\;hydraulic\;conductivity},\;K_{\rm s}\;{\rm (cm/h)}$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [ 0.5, -0.13, 0],'fontweight','normal');
            case 3
                plot1 = plot(eval(ax_str),opt_DREAM_SWIG(:,3,1),opt_DREAM_SWIG(:,3,2),'s','color',color_scatter',...
                    'linewidth',1.5,'markersize',6,'markerfacecolor','w'); hold on
                ylabel(eval(ax_str),'${\rm Coefficient},\;\beta\;(-)$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [-0.14, 0.5, 0],'fontweight','normal');
                xlabel(eval(ax_str),'${\rm Coefficient},\;\beta\;(-)$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [ 0.5, -0.13, 0],'fontweight','normal');
        end
        plot2 = line(eval(ax_str),[0 y_max(zz)],[0 y_max(zz)],'color','k','linewidth',1);
        
        %         plot2 = plot(eval(ax_str),opt_LM_SWIG(idx_data,zz,1),opt_LM_SWIG(idx_data,zz,2),'s','color',color_scatter,...
        %             'linewidth',1.5,'markersize',6,'markerfacecolor',color_scatter);
        % plot uncertainty of coefficient beta
        axis(eval(ax_str),[y_min(zz) y_max(zz) y_min(zz) y_max(zz) ]);
        % Add (A), (B) and (C)
        set(gca,'fontsize',fontsize_axes);
        % Add labels on x and y-axis
        a = axis;
        set(eval(ax_str),'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 1]); %
        set(eval(ax_str),'TickDir','out','clipping','off');
        % Now determine xticks and xticklabels and minorxticks
        switch zz
            case 3
                xtickformat('%.1f');
            case {1,2}
                xtickformat('%.0f');
        end
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
        if zz == 3
            ytickformat('%.1f')
        end
        % Now create minor ticks - y axis
        %        eval(strcat(ax_str,'.XAxis.MinorTickValues = xminortick'));
        axis(eval(ax_str),[y_min(zz) y_max(zz) y_min(zz) y_max(zz) ]);
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
            [lega,legb,legc,legd] = legend([plot1 plot2],{'${\rm \;\;SWIG\;data}$','${\rm \;\;1:1\;Line}$'},...
                'interpreter','latex','box','off','location','southeast','fontsize',fontsize_labels);
            set(legb(1),'Color',color_scatter); set(legb(4),'linewidth',1.5,'markersize',8); set(legb(2),'Color','k');
            set(legb(5),'linewidth',2.5);
            set(gcf, 'Renderer', 'painters');
        end
        axis(eval(ax_str),[y_min(zz) y_max(zz) y_min(zz) y_max(zz) ]);
        a = axis;
        % Add label (A), (B), or (C)
        text(eval(ax_str),a(1)+0.02*(a(2)-a(1)),a(3)+0.96*(a(4)-a(3)),char(index(zz)),...
            'fontsize',fontsize_labels+2,'interpreter','latex');
        set(gcf,'color','w');
        
    end
    set(legb(1),'Color',color_scatter); set(legb(3),'linewidth',2.5,'markersize',6); set(legb(2),'Color','k');
    % addpath to get nicest pdf from printing to pdf within matlab
    addpath('C:\Jasper Vrugt\Papers\Book\Chapter 1\MATLAB\Figure_resolution');
    % now use export_fig function from Figure_resolution directory
    evalstr = strcat('export_fig figure_DREAM_SWIG',{' '},'-','pdf');
    % now save file
    if save_fig == 1, eval(char(evalstr)); end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 12: SWIG BETA FROM DREAM: INFILTRATION AND TIME + UNCERTAINTY 95
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
color_scatter = [0 102 255]/255;
% Selection from data comparison after the fact (LM_SWIG figures)
color_conf_line = [0.6 0.6 0.6];
fontsize_axes = 21; fontsize_labels = 21;
if fig_12 == 1
    % Load results from LM_SWIG
    %ii_data = Data(:,2)>=10 & Data(:,3)>= 2; idx_data = Data(ii_data,1);
    %     ii_data = Data(:,3) >= 1; idx_data = Data(ii_data,1);
    %     ii_not_data = 1:n_soil; ii_not_data(ii_data) = [];
    idx_data = [1:n_soil];
    y_min = -0.04;
    y_max =  2.04;
    % Create figure name
    fig_name = 'Comparison of DREAM values of infiltration and time form Haverkamp';
    % Create figure
    figure('unit','inches','name',fig_name,'PaperOrientation','portrait','position',[0.3 0.5 12 10]);
    zz = 3
    ax_str = strcat('ax',num2str(zz));
    evalstr = strcat('ax',num2str(zz),{' '},'= axes(''units'',''inches'')'); eval(char(evalstr));
    evalstr = strcat('axpos',num2str(zz),{' '},'= [ 1.15 1.2 8 8 ]'); eval(char(evalstr));
    evalstr = strcat('set(','ax',num2str(zz),',','''position''',',axpos',num2str(zz),')'); eval(char(evalstr));
    % plot uncertainty of coefficient beta
    for i = 1:numel(idx_data)
        u = idx_data(i);
        plot3 = line(eval(ax_str),[opt_DREAM_SWIG(u,3,1) opt_DREAM_SWIG(u,3,1)],...
            [r_95_DREAM_SWIG(u,1,3,2) r_95_DREAM_SWIG(u,4,3,2)],...
            'color',color_conf_line,'linewidth',0.1);
        plot4 = line(eval(ax_str),[r_95_DREAM_SWIG(u,1,3,1) r_95_DREAM_SWIG(u,4,3,1)],...
            [opt_DREAM_SWIG(u,3,2) opt_DREAM_SWIG(u,3,2)],...
            'color',color_conf_line,'linewidth',0.1);
    end
    % Add one to one plot
    plot2 = line(eval(ax_str),[0 2],[0 2],'color','k','linewidth',1); hold on
    % Plot all data of SWIG
    plot1 = plot(eval(ax_str),opt_DREAM_SWIG(:,3,1),opt_DREAM_SWIG(:,3,2),'s','color',color_scatter',...
        'linewidth',1.5,'markersize',6,'markerfacecolor','w');
    ylabel(eval(ax_str),'${\rm Coefficient},\;\beta\;(-)$','interpreter',...
        'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
        [-0.1, 0.5, 0],'fontweight','normal');
    xlabel(eval(ax_str),'${\rm Coefficient},\;\beta\;(-)$','interpreter',...
        'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
        [ 0.5, -0.085, 0],'fontweight','normal');
    %     plot2 = plot(eval(ax_str),opt_LM_SWIG(idx_data,zz,1),opt_LM_SWIG(idx_data,zz,2),'s','color',color_scatter,...
    %         'linewidth',1.5,'markersize',6,'markerfacecolor',color_scatter);
    axis(eval(ax_str),[y_min y_max y_min y_max ]);
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
    dy = a(3) - 0.045*(a(4)-a(3));
    for rr = 1:numel(xtick)
        text(eval(ax_str),xtick(rr),dy,xticklabel(rr),'fontsize',fontsize_axes,'HorizontalAlignment', 'center');
    end
    % Now use the same labels on y-axis
    set(eval(ax_str),'ytick',xtick,'yticklabel',xticklabel);
    %    ytickformat('%.1f')
    % Now create minor ticks - y axis
    %        eval(strcat(ax_str,'.XAxis.MinorTickValues = xminortick'));
    % set box property to off and remove background color
    set(eval(ax_str),'box','off','color','none')
    % create new, empty axes with box but without ticks
    h2 = axes('Position',get(eval(ax_str),'Position'),'box','on','xtick',[],'ytick',[]);
    % set original axes as active
    % axes(eval(ax_str))
    % link axes in case of zooming and make rounding box square
    linkaxes([eval(ax_str) h2]);
    % add legend
    [lega,legb,legc,legd] = legend([plot1 plot2 plot4],{'${\rm \;\;SWIG\;data}$',...
        '${\rm \;\;1:1\;Line}$','${\rm \;\;95\%\;intervals}$'},'interpreter','latex','box','off','location',...
        'northeast','fontsize',fontsize_labels);
    set(legb(1),'Color',color_scatter); set(legb(2),'Color','k'); set(legb(3),'Color',color_conf_line);
    set(legb(5),'linewidth',1.5,'markersize',8);
    set(legb(6),'linewidth',2.5);
    set(legb(7),'linewidth',1);
    set(gcf, 'Renderer', 'painters');
    leg_pos = get(lega,'Position');
    set(lega,'Position',[leg_pos(1) + 0.22 leg_pos(2:4)]);
    axis(eval(ax_str),[y_min y_max y_min y_max ]);
    a = axis;
    % Add label (A), (B), or (C)
    % text(eval(ax_str),a(1)+0.02*(a(2)-a(1)),a(3)+0.98*(a(4)-a(3)),char(index(zz)),...
    %     'fontsize',fontsize_labels+2,'interpreter','latex');
    set(gcf,'color','w');
    
    set(legb(1),'Color',color_scatter); set(legb(2),'Color','k'); set(legb(3),'Color',color_conf_line);
    set(legb(5),'linewidth',1.5,'markersize',8);
    set(legb(6),'linewidth',2.5);
    set(legb(7),'linewidth',1);
    
    % addpath to get nicest pdf from printing to pdf within matlab
    addpath('C:\Jasper Vrugt\Papers\Book\Chapter 1\MATLAB\Figure_resolution');
    % now use export_fig function from Figure_resolution directory
    evalstr = strcat('export_fig figure_DREAM_SWIG_beta_conf',{' '},'-','pdf');
    % now save file
    if save_fig == 1, eval(char(evalstr)); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 13: CREATE HISTOGRAMS OF S, Ks and beta (SELECTION)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%idx_sel = [ 1 6 8 12 27 34 45 48 ];
% 1) Sand
% 6) Sandy loam
% 8) Sandy clay loam
% 12) Loam
% 27) Silt loam
% 34) Clay loam
% 45) Silty clay loam
% 48) Clay
%idx_sel2 = [ 1 8 27 45 48 ];

% Load results from LM_SWIG
%load SWIG_646_LM2.mat
%load SWIG_646_LM_test_Jac
eval(char(load_file));
load DREAM_SWIG.mat
% Now select
idx_sel_ORIG = [ 1 6 8 12 27 34 45 48 ]; % FROM ORIGINAL 50 SAMPLES
idx_sel_ORIG2 = [ 1 8 27 45 48 ];
idx_selSWIG = [ 1 2 3 42 43 45 114 115 122 151 164:172 176 191 197 198 201 224 ...
    283 284 285 300 352 253 259 371 379 380 384 387 388 433 434 437 440 444 ...
    463 476 478 518 519 523 524];
% 371 on our list;
idx_selSWIG = idx_selSWIG(idx_sel_ORIG2);

index2 = {'(A','(B','(C','(D','(E','(F','(G','(H','(I','(J','(K','(L'};
fontsize_axes = 18; fontsize_labels = 18;
if fig_13 == 1
    % define axis labels for all figures
    %     x_min = [1.03 0.209 1.52 1.478 0.28 1.25 9.24 31 0.6 1.69 0.355    1.51 0.512 0.0565 1.90 ];
    %     x_max = [1.05 0.219 1.76 1.502 0.32 1.65 9.44 32.6 0.8 1.72 0.395 1.87 0.526 0.0595 2.00 ];
    %     x_sd = [ 2 2 2 2 2 2 1 1 1 2 2 2 2 3 2 ];
    %index = {'(A','(B','(C','(D','(E','(F','(G','(H'};
    %     x_delta = (x_max - x_min)/2;
    %     x_edge = x_delta/20;
    ct = 1;
    n_graphs = 5;
    fig_num = 1; rem_fig = 1:n_graphs:n_soil;
    for soil_type = 1:n_graphs
        % Create figure name
        if rem(soil_type,rem_fig(fig_num)) == 0
            fig_name = char(strcat('Plot of posterior haverkamp parameters:',{' '},num2str(fig_num)));
            % Create figure
            figure('unit','inches','name',fig_name,'PaperOrientation','portrait','position',[0.1 0.01 19.7 14.1]);
            % Set cf to 1;
            cf = 1; fig_num = fig_num + 1; fig_num = min(fig_num,ceil(n_soil/n_graphs));
        end
        % Loop over parameters
        for par = 1:3
            % define axis position
            ax_str = strcat('ax',num2str(idx_selSWIG(soil_type)),num2str(par));
            evalstr = strcat('ax',num2str(idx_selSWIG(soil_type)),num2str(par),{' '},'= axes(''units'',''inches'')'); eval(char(evalstr));
            evalstr = strcat('axpos',num2str(idx_selSWIG(soil_type)),num2str(par),{' '},'= [ 1.0 + (cf-1)*3.6, 7.8-(par-1)*3.47 2.6 2.5 ]'); eval(char(evalstr));
            % set axes
            evalstr = strcat('set(','ax',num2str(idx_selSWIG(soil_type)),num2str(par),',','''position''',',axpos',...
                num2str(idx_selSWIG(soil_type)),num2str(par),')'); eval(char(evalstr));
            %                    set(ax1,'position',axpos1);
            % plot marginal distribution
            P = P_DREAM_SWIG{idx_selSWIG(soil_type),1};
            [Np,Xp] = hist(P(:,par),13); Nps = Np./trapz(Xp,Np); % Integrates to one
            Npo = Nps/max(Nps); % Now maximum density is 1
            bar(eval(ax_str),Xp,Npo,'edgecolor',[0.4 0.4 0.4],'facecolor','w','linewidth',1); hold on;
            % Add optimum - DREAM
            plot(eval(ax_str),opt_DREAM_SWIG(idx_selSWIG(soil_type),par,1),0,'rx','linewidth',3,'markersize',15);
            % Add optimum - LM
            green_LM = [34 139 34]/255; %[0.4660 0.6740 0.1880];
            plot(eval(ax_str),opt_LM_SWIG(idx_selSWIG(soil_type),par,1),0,'bo','color',green_LM,'linewidth',3,...
                'markersize',10,'markerfacecolor','w');
            
            % set to old axis
            axis([xlim 0 1.08]); a = axis;
            % Add common y-labels
            set(eval(ax_str),'XMinorTick','on','YMinorTick','on','TickLength',[0.02, 1],'Clipping','on'); %
            set(eval(ax_str),'TickDir','out');
            %if soil_type == 1
            set(eval(ax_str),'ytick',0:0.2:1.0,'yticklabel',{'0','0.2','0.4','0.6','0.8','1.0'});
            % Now determine xticks and xticklabels and minorxticks
            if par == 3
                xtickformat('%0.2f');
            end
            xticklabel = get(eval(ax_str),'xticklabel'); xtick_new = '';
            xtick = get(eval(ax_str),'xtick');
            %set(eval(ax_str),'xticklabel',xtick_new);
            set(eval(ax_str),'xticklabel',[]);
            % %             if par == 1
            % %                 switch soil_type
            % %                     case 4
            % %                         xticklabel = {'4.5','5.0','5.5'};
            % %                     case 5
            % %                         xtick = [9 10 11]; xticklabel = {'9','10','11'};
            % %                 end
            % %             end
            % %             if par == 2
            % %                 switch soil_type
            % %                     case 4
            % %                         xticklabel = {'1.0','1.5','2.0','2.5'};
            % %                 end
            % %             end
            % %             if par == 3
            % %                 switch soil_type
            % %                     case 1
            % %                         xticklabel = {'0.0','0.5','1.0','1.5'};
            % %                     case 2
            % %                         xticklabel = {'0.0','0.5','1.0'};
            % %                     case 3
            % %                         xticklabel = {'0.0','0.5','1.0'}; xtick = [0:0.5:1];
            % %                     case 4
            % %                         xticklabel = {'0.0','0.5','1.0','1.5','2.0'};
            % %                     case 5
            % %                         xticklabel = {'0.0','0.5','1.0','1.5','2.0'};
            % %                 end
            % %             end
            % now plot manually
            dy = a(3) - 0.095*(a(4)-a(3));
            
            for rr = 1:numel(xtick)
                text(eval(ax_str),xtick(rr),dy,xticklabel(rr),'fontsize',fontsize_axes,'HorizontalAlignment', 'center');
            end
            eval(strcat(ax_str,'.YAxis.MinorTickValues = 0:0.1:1'));
            % set fontsize of numbers to 15
            set(eval(ax_str),'fontsize',fontsize_axes);
            % Add labels on x and y-axis
            switch par
                case 1
                    tx1 = xlabel(eval(ax_str),'$S\;({\rm cm/h}^{1/2})$','interpreter','latex','fontsize',fontsize_labels,'Units',...
                        'normalized','Position',[0.5, -0.2, 0],'fontweight','normal');
                case 2
                    tx2 = xlabel(eval(ax_str),'$K_{\rm s}\;({\rm cm/h})$','interpreter','latex','fontsize',fontsize_labels,'Units',...
                        'normalized','Position',[0.5, -0.2, 0],'fontweight','normal');
                case 3
                    tx3 = xlabel(eval(ax_str),'$\beta\;(-)$','interpreter','latex','fontsize',fontsize_labels,'Units',...
                        'normalized','Position',[0.5, -0.2, 0],'fontweight','normal');
            end
            if ismember(soil_type,1)
                ty = ylabel(eval(ax_str),'Empirical density','interpreter','latex','fontsize',fontsize_labels,'Units',...
                    'normalized','Position',[-0.23, 0.5, 0],'fontweight','normal');
            else
                % make y-axis invisible
                % set(eval(ax_str),'YColor','none');
            end
            % Add soil type
            if par == 1
                % Add label with name of soil
                str = data_SWIG{idx_selSWIG(soil_type),3};
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
                %                str_p = char(strcat(index(soil_type),{' '},str_f));
                text(eval(ax_str),a(1)+0.5*(a(2)-a(1)),1.08,char(str_f),'interpreter','latex','fontweight',...
                    'bold','fontsize',fontsize_labels,'HorizontalAlignment', 'center');
                %str_f = char(str_f);
                %text(eval(ax_str),a(1)+0.5*(a(2)-a(1)),1.1,str_f,'interpreter','latex','fontweight',...
                %    'bold','fontsize',fontsize_labels,'HorizontalAlignment', 'center');
            end
            str_p = strcat(index2(soil_type),num2str(par),')');
            text(eval(ax_str),a(1)+0.02*(a(2)-a(1)),0.95,str_p,'fontsize',fontsize_labels,'interpreter','latex');
            %           text(eval(ax_str),a(1)+0.02*(a(2)-a(1)),0.95,strcat(num2str(idx_sel2(soil_type)),':',num2str(par),')'),'fontsize',fontsize_labels,'interpreter','latex');
            
            % set box property to off and remove background color
            set(eval(ax_str),'box','off','color','none');
            % adjust xtick for beta
            % % %             if par == 3
            % % %                 switch soil_type
            % % %                     case 1
            % % %                         axis([-0.05 1 0 1.08]);
            % % %                     case 2
            % % %                         axis([-0.05 1.2 0 1.08]);
            % % %                     case 3
            % % %                         axis([-0.05 1 0 1.08]);
            % % %                     case 4
            % % %                         axis([-0.05 2 0 1.08]);
            % % %                     case 5
            % % %                         axis([-0.05 2 0 1.08]);
            % % %                 end
            % % %             end
            % create new, empty axes with box but without ticks
            ct = ct + 1;
        end
        cf = cf + 1;
        % check before we move to new figure
        if rem(soil_type+1,rem_fig(fig_num)) == 0
            set(gcf,'color','w');
            % set(legb(1),'Color','r'); set(legb(4),'linewidth',1.5,'markersize',8); set(legb(5),'XData',[0.0603 0.1838],'linewidth',3); set(legb(2),'Color','b');
            evalstr = strcat('export_fig figure_13_',num2str(fig_num-1),{' '},'-','pdf');
            % now save file
            if save_fig == 1, eval(char(evalstr)); end
        end
    end
    set(gcf,'color','w');
    addpath('C:\Jasper Vrugt\Papers\Book\Chapter 1\MATLAB\Figure_resolution');
    % now use export_fig function from Figure_resolution directory
    evalstr = strcat('export_fig figure_13',{' '},'-','pdf');
    % now save file
    if save_fig == 1, eval(char(evalstr)); end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 14: BIVARIATE SCATTER PLOTS; (S,Ks); (S,beta); (Ks,beta) OF
% HAVERKAMP POSTERIOR DISTRIBUTION        if zz < 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
color_scatter = [0 102 255]/255;
%color_scatter = [100 161 244]/255;
if fig_14 == 1
    % Create figure name
    fig_name = 'Plot of bivariate posterior samples of parameters Haverkamp infiltration equation';
    % Create figure
    figure('unit','inches','name',fig_name,'PaperOrientation','portrait','position',[0.1 0.5 19.7 14]);
    % Now make plot of posterior distribution
    cf = 1; index = {'(A','(B','(C','(D','(E','(F','(G','(H'};
    % define axis labels for all figures
    % x_min = [1.03 0.209 1.52 1.478 0.28 1.25 9.24 31 0.6 3.89 4.41 1.02 1.69 0.355    1.51 0.512 0.0565 1.90 ];
    % x_max = [1.05 0.219 1.76 1.502 0.32 1.65 9.44 32.6 0.8 3.93 4.49 1.08 1.72 0.395 1.87 0.526 0.0595 2.00 ];
    % x_delta = (x_max - x_min)/2;
    % x_edge = x_delta/20;
    x_sd = [ 2 2 2 2 2 2 1 1 0 2 2 2 2 2 3 ];
    y_sd = [ 3 1 1 2 1 1 1 2 2 2 1 1 3 2 2 ];
    % %     x_sd = [ 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 3 3 3 ];
    % %     y_sd = [ 3 1 1 2 1 1 1 2 1 2 2 2 2 1 1 3 2 2 ];
    ct = 1; fig_num = 1; par_pairs = [ 1 2; 1 3; 2 3 ];
    %index = {'(A)','(B)','(C)','(D)','(E)','(F)','(G)','(H)','(I)','(J)','(K)','(L)','(M)'};
    n_graphs = 5;
    fig_num = 1; rem_fig = 1:n_graphs:n_soil;
    for soil_type = 1:n_soil
        % Create figure name
        if rem(soil_type,rem_fig(fig_num)) == 0
            fig_name = char(strcat('Bivariate scatter plots of posterior samples:',{' '},num2str(fig_num)));
            % Create figure
            figure('unit','inches','name',fig_name,'PaperOrientation','portrait','position',[0.1 0.4 19.7 14]);
            % Set cf to 1;
            cf = 1; fig_num = fig_num + 1; fig_num = min(fig_num,ceil(n_soil/n_graphs));
        end
        % Now create axes - for different pars + plot histogram
        for par = 1:3
            % define axis position
            ax_str = strcat('ax',num2str(soil_type),num2str(par));
            evalstr = strcat('ax',num2str(soil_type),num2str(par),{' '},'= axes(''units'',''inches'')'); eval(char(evalstr));
            evalstr = strcat('axpos',num2str(soil_type),num2str(par),{' '},'= [ 1.4 + (cf-1)*3.9, 7.85-(par-1)*3.5 2.4 2.5 ]'); eval(char(evalstr));
            % set axes
            evalstr = strcat('set(','ax',num2str(soil_type),num2str(par),',','''position''',',axpos',num2str(soil_type),num2str(par),')'); eval(char(evalstr));
            % plot bivariate samples
            plot(eval(ax_str),Ppost_SWIG(:,par_pairs(par,1),soil_type,1),Ppost_SWIG(:,par_pairs(par,2),soil_type,1),'s',...
                'color',color_scatter,'markersize',5,'linewidth',1,'markerfacecolor','w'); hold on
            % set fontsize of numbers to 15
            set(eval(ax_str),'fontsize',fontsize_axes);
            % now store the values of a as sometimes plotting goes wrong for unknown reason (MATLAB bug)
            a = axis;
            if par == 1, a1 = a; end
            if par == 2, a(1:2) = a1(1:2); a2 = a; end
            if par == 3, a(3:4) = a2(3:4); a(1:2) = a1(3:4); end
            % RESET AXIS
            axis(a);
            % add linear regression line
            reg = polyfit(Ppost_SWIG(:,par_pairs(par,1),soil_type,1),Ppost_SWIG(:,par_pairs(par,2),soil_type,1),1);
            % now create (x,y) line
            x = [a(1) a(2)]; y = reg(1) * x + reg(2);
            plot(eval(ax_str),x,y,'k--','linewidth',1.5);
            % plot optimum
            plot(eval(ax_str),opt_SWIG(soil_type,par_pairs(par,1),1),opt_SWIG(soil_type,par_pairs(par,2),1),'rx','markersize',13,'linewidth',3,'markerfacecolor','w');
            % set to old axis
            %                    axis([x_min(ct)-x_edge(ct) x_max(ct)+x_edge(ct) 0 1.08]);
            % Add common y-labels
            set(eval(ax_str),'XMinorTick','on','YMinorTick','on','TickLength',[0.02, 1],'Clipping','on'); %
            set(eval(ax_str),'TickDir','out');
            xticklabel = get(eval(ax_str),'xticklabel');
            xtick = get(eval(ax_str),'xtick');
            set(eval(ax_str),'xticklabel',[]);
            % now plot manually - for all except 12
            dy = a(3) - 0.095*(a(4)-a(3));
            for rr = 1:numel(xtick)
                text(eval(ax_str),xtick(rr),dy,xticklabel(rr),'fontsize',fontsize_axes,'HorizontalAlignment', 'center');
            end
            % %             % Now determine yticks and yticklabels and minoryticks
            % %             ytick = get(eval(ax_str),'yticklabel'); ytick_new = '';
            % %             % determine number of significant digits x-axis
            % %             sd_y = strcat('%4.',num2str(y_sd(ct)),'f');
            % %             % now adjust xticks to match significant digits
            % %             for z = 1:numel(ytick), ytick_new{z} = num2str(str2num(char(ytick(z))),sd_y); end
            % %             % now set the xtick values
            % %             set(eval(ax_str),'yticklabel',ytick_new);
            % We make an adjustment for plot 12 - imperfect for
            % some reason
            % %             if ct == 12
            % %                 xtick = 0.34:0.02:0.40;
            % %                 % must specify between {} to get 2 0.40 instead of 0.4
            % %                 % set(eval(ax_str),'xtick',xtick,'xticklabel',{'0.34','0.36','0.38','0.40'});
            % %                 xtick_new = {'0.34','0.36','0.38','0.40'};
            % %                 set(eval(ax_str),'xtick',xtick,'xticklabel',[]);
            % %                 dy = a(3) - 0.095*(a(4)-a(3));
            % %                 for rr = 1:numel(xtick)
            % %                     text(eval(ax_str),xtick(rr),dy,char(xtick_new{rr}),'fontsize',fontsize_axes,'HorizontalAlignment', 'center');
            % %                 end
            % %                 xminortick = 0.34:0.005:0.40; xminortick = xminortick(~ismember(xminortick,xtick));
            % %                 eval(strcat(ax_str,'.XAxis.MinorTickValues = xminortick'));
            % %             end
            yticklabel = get(eval(ax_str),'yticklabel');
            ytick = get(eval(ax_str),'ytick');
            nr_sd = 1;
            for zz = 1:numel(ytick)
                nr_zz = sigdigits(str2num(char(yticklabel(zz))));
                nr_sd = max(nr_zz,nr_sd);
            end
            if max(ytick) >= 100
                sd_pen = 3;
            elseif max(ytick) >= 10
                sd_pen = 2;
            elseif max(ytick) >= 1
                sd_pen = 1;
            else
                sd_pen = 0;
            end
            nr_sd = max(0,nr_sd-sd_pen);
            ytickformat(char(strcat('%.',num2str(nr_sd),'f')));
            
            % Add labels on x and y-axis
            switch par_pairs(par,1)
                case 1
                    xlabel(eval(ax_str),'$S\;({\rm cm/h}^{1/2})$','interpreter','latex','fontsize',fontsize_labels,'Units',...
                        'normalized','Position',[0.5, -0.2, 0],'fontweight','normal');
                case 2
                    xlabel(eval(ax_str),'$K_{\rm s}\;({\rm cm/h})$','interpreter','latex','fontsize',fontsize_labels,'Units',...
                        'normalized','Position',[0.5, -0.2, 0],'fontweight','normal');
                case 3
                    xlabel(eval(ax_str),'$\beta\;(-)$','interpreter','latex','fontsize',fontsize_labels,'Units',...
                        'normalized','Position',[0.5, -0.2, 0],'fontweight','normal');
            end
            % ylabels - placing depends on number of values of
            % yticks - this partly depends on significant digits
            % %             mult = numel(ytick_new{1});
            mult = 2;
            switch par_pairs(par,2)
                case 1
                    ylabel(eval(ax_str),'$S\;({\rm cm/h}^{1/2})$','interpreter','latex','fontsize',fontsize_labels,'Units',...
                        'normalized','Position',[-0.23 - (mult-2)*0.05, 0.5, 0],'fontweight','normal');
                case 2
                    ylabel(eval(ax_str),'$K_{\rm s}\;({\rm cm/h})$','interpreter','latex','fontsize',fontsize_labels,'Units',...
                        'normalized','Position',[-0.23 - (mult-2)*0.05, 0.5, 0],'fontweight','normal');
                case 3
                    ylabel(eval(ax_str),'$\beta\;(-)$','interpreter','latex','fontsize',fontsize_labels,'Units',...
                        'normalized','Position',[-0.23 - (mult-2)*0.05, 0.5, 0],'fontweight','normal');
            end
            % Add soil type
            if par == 1
                % Add label with name of soil
                str = data_SWIG{soil_type,2};
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
                str_f = char(str_f);
                text(eval(ax_str),a(1)+0.2*(a(2)-a(1)),a(3)+0.95*(a(4)-a(3)),str_f,'interpreter','latex','fontweight',...
                    'bold','fontsize',fontsize_labels,'HorizontalAlignment', 'left');
            end
            axis(a);
            if ct <= 12
                text(eval(ax_str),a(1)+0.02*(a(2)-a(1)),a(3)+0.95*(a(4)-a(3)),strcat(char(index(cf)),num2str(par),')'),...
                    'fontsize',fontsize_labels,'interpreter','latex');
            else
                text(eval(ax_str),a(1)+0.02*(a(2)-a(1)),a(3)+0.08*(a(4)-a(3)),strcat(char(index(cf)),num2str(par),')'),...
                    'fontsize',fontsize_labels,'interpreter','latex');
            end
            % set box property to off and remove background color
            set(eval(ax_str),'box','off','color','none')
            % add legend to first plot
            if ct == 2
                [lega,legb,legc] = legend(eval(ax_str),{'','${\rm Regression\;line}$','${\rm Optimum}$'},'interpreter','latex','box','off','color','w',...
                    'location','best','fontsize',fontsize_labels);
                set(legb(1),'color','w');
                set(legb(3),'color','r');
                set(legb(6),'XData',[0.0603 0.1838],'linewidth',3,'linestyle',':')
                set(legb(2),'linewidth',3);
                %set(lega(1),
                legpos = get(lega,'position');
                set(lega,'position',[legpos(1)+0.03 legpos(2)-0.012 legpos(3:4)]);
            end
            ct = ct + 1;
        end
        cf = cf + 1;
        % check before we move to new figure
        if rem(soil_type+1,rem_fig(fig_num)) == 0
            set(gcf,'color','w');
            %set(legb(1),'Color','r'); set(legb(4),'linewidth',1.5,'markersize',8); set(legb(5),'XData',[0.0603 0.1838],'linewidth',3); set(legb(2),'Color','b');
            evalstr = strcat('export_fig figure_14_',num2str(fig_num-1),{' '},'-','pdf');
            % now save file
            if save_fig == 1, eval(char(evalstr)); end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 15: NUMBER OF ITERATIONS NEWTON'S METHOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
color_1 = [204 0 0]/255;
color_2 = [0   0 205]/255;
color_3 = [34 139 34]/255;
fontsize_axes = 21; fontsize_labels = 21;
if fig_15 == 1
    % Define default parameter values
    S = 2; Ks = 1; B = 1.5; Ki = 1e-1; eta = [S Ks B Ki];
    % Now compute number of iterations
    plugin.t = [0.0:0.01:20]'; iter = []; I = [];
    for zz = 1:3
        [I(:,zz),i,flag,iter(:,zz)] = Haverkamp_IK(eta,plugin,1e-12,20,zz);
    end
    if plugin.t(1) == 0
        iter(1,1:3) = nan;
    end
    I
    I = I(:,1);
    def1 = max(I);
    x_max = 20; %max(I);
    y_min = [ 2.8 ];     y_max = [ 7.2 ];
    % Create figure name
    fig_name = 'Number of iterations Newton''s method: infiltration form';
    % Create figure
    figure('unit','inches','name',fig_name,'PaperOrientation','portrait','position',[0.3 0.5 10 6.5]);
    % Define index
    index = {'(A)','(B)','(C)','(D)','(E)','(F)','(G)','(H)','(I)','(J)','(K)','(L)','(M)'};
    % Plot 3 figures
    zz = 1;
    ax_str = strcat('ax',num2str(zz));
    evalstr = strcat('ax',num2str(zz),{' '},'= axes(''units'',''inches'');'); eval(char(evalstr));
    evalstr = strcat('axpos',num2str(zz),{' '},'= [ 1 1.3 8 4 ];'); eval(char(evalstr));
    evalstr = strcat('set(','ax',num2str(zz),',','''position''',',axpos',num2str(zz),');'); eval(char(evalstr));
    plot1 = plot(eval(ax_str),plugin.t,iter(:,1),'color',color_1,'linewidth',3,'linestyle',':'); hold on
    plot2 = plot(eval(ax_str),plugin.t,iter(:,2),'color',color_2,'linewidth',3,'linestyle',':');
    plot3 = plot(eval(ax_str),plugin.t,iter(:,3),'color',color_3,'linewidth',3,'linestyle',':');
    xlabel(eval(ax_str),'${\rm Time},\;t\;\;({\rm h})$','interpreter','latex','fontsize',fontsize_labels,'Units',...
        'normalized','Position',[0.5, -0.16, 0],'fontweight','normal');
    ylabel(eval(ax_str),'${\rm \#\;of\;Newton\;iterations}$','interpreter',...
        'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
        [-0.08, 0.5, 0],'fontweight','normal');
    set(eval(ax_str),'ytick',[3:1:7],'yticklabel',[3:1:7]);
    % %         end
    axis([0 x_max y_min(zz) y_max(zz) ]);
    % Add (A), (B) and (C)
    set(gca,'fontsize',fontsize_axes);
    % Add labels on x and y-axis
    %         xlabel(eval(ax_str),'${\rm Cumulative\;infiltration}, I \; {\rm (\,{\rm cm}\,)}$','interpreter','latex',...
    %             'fontsize',fontsize_labels,'Units','normalized','Position',[0.5, -0.14, 0],'fontweight','normal');
    a = axis;
    set(eval(ax_str),'XMinorTick','on','YMinorTick','off','TickLength',[0.02, 1]); %
    set(eval(ax_str),'TickDir','out','clipping','off');
    % Now determine xticks and xticklabels and minorxticks
    xticklabel = get(eval(ax_str),'xticklabel');
    xtick = get(eval(ax_str),'xtick');
    % remove current labels - and replace with own values closer to axis
    set(eval(ax_str),'xticklabel',[]);
    % now plot manually - for all except 12
    dy = a(3) - 0.085*(a(4)-a(3));
    for rr = 1:numel(xtick)
        text(eval(ax_str),xtick(rr),dy,xticklabel(rr),'fontsize',fontsize_axes,'HorizontalAlignment', 'center');
    end
    %        xtickformat('%.0f')
    ytickformat('%.0f')
    % Now create minor ticks - y axis
    %        eval(strcat(ax_str,'.XAxis.MinorTickValues = xminortick'));
    axis(eval(ax_str),[0 x_max a(3) y_max ]);
    % set box property to off and remove background color
    set(eval(ax_str),'box','off','color','none')
    % create new, empty axes with box but without ticks
    h2 = axes('Position',get(eval(ax_str),'Position'),'box','on','xtick',[],'ytick',[]);
    % set original axes as active
    % link axes in case of zooming and make rounding box square
    linkaxes([eval(ax_str) h2]);
    % add legend
    if zz == 1
        [h_leg,legb,legc,legd] = legend([plot1 plot2 plot3],{'$\;I_{(0)} = I_{\rm low}$',...
            '$\;I_{(0)} = \frac{1}{2}(I_{\rm low} + I_{\rm up})$',...
            '$\;I_{(0)} = I_{\rm up}$'},'interpreter','latex','box','off','location',...
            'northeast','fontsize',fontsize_labels);
        set(legb(1),'color',color_1);
        set(legb(2),'color',color_3);
        set(legb(3),'color',color_2);
        set(legb(4),'linewidth',4);
        set(legb(6),'linewidth',4);
        set(legb(8),'linewidth',4);
        HeightScaleFactor = 1.5;
        NewHeight = h_leg.Position(4) * HeightScaleFactor;
        h_leg.Position(2) = h_leg.Position(2) - (NewHeight - h_leg.Position(4));
        h_leg.Position(4) = NewHeight;
        % set(eval(ax_str),'PaperPositionMode','auto');
        set(gcf, 'Renderer', 'painters')
    end
    axis(eval(ax_str),[0 x_max y_min y_max ]);
    a = axis;
    % Add label (A), (B), or (C)
    %         text(eval(ax_str),a(1)+0.02*(a(2)-a(1)),a(3)+0.955*(a(4)-a(3)),char(index(zz)),...
    %             'fontsize',fontsize_labels+2,'interpreter','latex');
    set(gcf,'color','w');
    
    if zz == 1
        % add 2nd x-axis with infiltration
        yloc = a(3)+1.02*(a(4)-a(3));
        line(eval(ax_str),[0 x_max],[yloc yloc],'color','k','linewidth',0.5);
        % now set vertical lines at right places
        mult = x_max / max(def1); ttick = [ 0 4 8 12 16 20 ]; tticklabel = {'0','4','8','12','16','20'};
        dy2 = yloc + 0.07*(a(4)-a(3));
        for rr = 1:numel(ttick)
            line(eval(ax_str),[ttick(rr) ttick(rr)]*mult,[yloc yloc+0.03*(a(4)-a(3))],'color','k','linewidth',0.5);
            % add label
            text(eval(ax_str),mult*ttick(rr),dy2,tticklabel(rr),'fontsize',fontsize_axes,'HorizontalAlignment', 'center');
        end
        % minor ticks
        minorttick = [0:21]; minorttick = minorttick(~ismember(minorttick,ttick));
        for rr = 1:numel(minorttick)
            line(eval(ax_str),[minorttick(rr) minorttick(rr)]*mult,[yloc yloc+0.015*(a(4)-a(3))],'color','k','linewidth',0.5);
        end
        % add axes label
        dy3 = yloc + 0.18*(a(4)-a(3));
        text(eval(ax_str),10,dy3,'${\rm Cumulative\;infiltration},\;I\;\;({\rm cm})$','fontsize',fontsize_labels,'HorizontalAlignment', 'center','interpreter','latex');
    end
    set(legb(1),'color',color_1);
    set(legb(2),'color',color_2);
    set(legb(3),'color',color_3);
    set(legb(4),'linewidth',4);
    set(legb(6),'linewidth',4);
    set(legb(8),'linewidth',4);
    % addpath to get nicest pdf from printing to pdf within matlab
    addpath('C:\Jasper Vrugt\Papers\Book\Chapter 1\MATLAB\Figure_resolution');
    % now use export_fig function from Figure_resolution directory
    evalstr = strcat('export_fig figure_newton_iter',{' '},'-','pdf');
    % now save file
    if save_fig == 1, eval(char(evalstr)); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 16: NUMBER OF ITERATIONS NEWTON'S METHOD AS FUNCTION OF DELTA T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
color_1 = [204 0 0]/255;
color_2 = [0   0 205]/255;
color_3 = [34 139 34]/255;
fontsize_axes = 21; fontsize_labels = 21;
if fig_16 == 1
    % Load results from numerr.m file
    load iter_Newton.mat
    % Now compute number of iterations
    x_max = max(dt); %max(I);
    y_min = [ 1.5 ];     y_max = [ 3.5 ];
    % Create figure name
    fig_name = 'Number of iterations Newton''s method as function dt: infiltration form';
    % Create figure
    figure('unit','inches','name',fig_name,'PaperOrientation','portrait','position',[0.3 0.5 10 6.5]);
    % Plot 3 figures
    zz = 1;
    ax_str = strcat('ax',num2str(zz));
    evalstr = strcat('ax',num2str(zz),{' '},'= axes(''units'',''inches'');'); eval(char(evalstr));
    evalstr = strcat('axpos',num2str(zz),{' '},'= [ 1.3 1.3 8 4 ];'); eval(char(evalstr));
    evalstr = strcat('set(','ax',num2str(zz),',','''position''',',axpos',num2str(zz),');'); eval(char(evalstr));
    plot1 = plot(eval(ax_str),dt,sumK(:,1),'color',color_1,'linewidth',3,'linestyle',':'); hold on
    plot2 = plot(eval(ax_str),dt,sumK(:,2),'color',color_2,'linewidth',3,'linestyle',':');
    plot3 = plot(eval(ax_str),dt,sumK(:,3),'color',color_3,'linewidth',3,'linestyle',':');
    xlabel(eval(ax_str),'${\rm Time\;step},\;\Delta t\;\;({\rm h})$','interpreter','latex','fontsize',fontsize_labels,'Units',...
        'normalized','Position',[0.5, -0.16, 0],'fontweight','normal');
    ylabel(eval(ax_str),'${\rm Avg.\;\#\;of\;Newton\;iterations}$','interpreter',...
        'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
        [-0.1, 0.5, 0],'fontweight','normal');
    set(eval(ax_str),'ytick',[1:.5:3.5],'yticklabel',[1:.5:3.5]);
    ytickformat('%.1f');
    % %         end
    axis([0 x_max y_min(zz) y_max(zz) ]);
    % Add (A), (B) and (C)
    set(gca,'fontsize',fontsize_axes);
    % Add labels on x and y-axis
    %         xlabel(eval(ax_str),'${\rm Cumulative\;infiltration}, I \; {\rm (\,{\rm cm}\,)}$','interpreter','latex',...
    %             'fontsize',fontsize_labels,'Units','normalized','Position',[0.5, -0.14, 0],'fontweight','normal');
    a = axis;
    set(eval(ax_str),'XMinorTick','on','YMinorTick','on','TickLength',[0.02, 1]); %
    set(eval(ax_str),'TickDir','out','clipping','off');
    % Now determine xticks and xticklabels and minorxticks
    xticklabel = get(eval(ax_str),'xticklabel');
    xtick = get(eval(ax_str),'xtick');
    % remove current labels - and replace with own values closer to axis
    set(eval(ax_str),'xticklabel',[]);
    % now plot manually - for all except 12
    dy = a(3) - 0.085*(a(4)-a(3));
    for rr = 1:numel(xtick)
        text(eval(ax_str),xtick(rr),dy,xticklabel(rr),'fontsize',fontsize_axes,'HorizontalAlignment', 'center');
    end
    %        xtickformat('%.0f')
    %    ytickformat('%.0f')
    % Now create minor ticks - y axis
    %        eval(strcat(ax_str,'.XAxis.MinorTickValues = xminortick'));
    axis(eval(ax_str),[0 x_max a(3) y_max ]);
    % set box property to off and remove background color
    set(eval(ax_str),'box','off','color','none')
    % create new, empty axes with box but without ticks
    h2 = axes('Position',get(eval(ax_str),'Position'),'box','on','xtick',[],'ytick',[]);
    % set original axes as active
    % link axes in case of zooming and make rounding box square
    linkaxes([eval(ax_str) h2]);
    % add legend
    if zz == 1
        [h_leg,legb,legc,legd] = legend([plot1 plot2 plot3],{'$\;I_{(0)} = I_{\rm low}$',...
            '$\;I_{(0)} = \frac{1}{2}(I_{\rm low} + I_{\rm up}\vphantom{\frac{A^{2}}{B}})$',...
            '$\;I_{(0)} = I_{\rm up}$'},'interpreter','latex','box','off','location',...
            'northeast','fontsize',fontsize_labels);
        set(legb(1),'color',color_1);
        set(legb(2),'color',color_3);
        set(legb(3),'color',color_2);
        set(legb(4),'linewidth',4);
        set(legb(6),'linewidth',4);
        set(legb(8),'linewidth',4);
        HeightScaleFactor = 1.5;
        NewHeight = h_leg.Position(4) * HeightScaleFactor;
        h_leg.Position(2) = h_leg.Position(2) - (NewHeight - h_leg.Position(4));
        h_leg.Position(4) = NewHeight;
        % set(eval(ax_str),'PaperPositionMode','auto');
        set(gcf, 'Renderer', 'painters')
    end
    axis(eval(ax_str),[0 x_max y_min y_max ]);
    a = axis;
    % Add label (A), (B), or (C)
    %         text(eval(ax_str),a(1)+0.02*(a(2)-a(1)),a(3)+0.955*(a(4)-a(3)),char(index(zz)),...
    %             'fontsize',fontsize_labels+2,'interpreter','latex');
    set(gcf,'color','w');
    set(legb(1),'color',color_1);
    set(legb(2),'color',color_2);
    set(legb(3),'color',color_3);
    set(legb(4),'linewidth',4);
    set(legb(6),'linewidth',4);
    set(legb(8),'linewidth',4);
    % addpath to get nicest pdf from printing to pdf within matlab
    addpath('C:\Jasper Vrugt\Papers\Book\Chapter 1\MATLAB\Figure_resolution');
    % now use export_fig function from Figure_resolution directory
    evalstr = strcat('export_fig figure_newton_iter_dt',{' '},'-','pdf');
    % now save file
    if save_fig == 1, eval(char(evalstr)); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 17: DREAM UNCERTAINTY VERSUS LM UNCERTAINTY + OPTIMUM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
color_scatter = [0 102 255]/255;
% Selection from data comparison after the fact (LM_SWIG figures)
idx_notsel = [ 5 27 30 48 57 60 61 72 74 76 77 78 80 95 100 109 110 111 113 ...
    116 120 121 126 127 128 131 144 145 147 175 176 177 187 192 193 196  ...
    202 203 204 221 222 225 226 228 275 276 269 302 ...
    303 305 306 371 363 377 385 386 390 393 394 400 410 411 421 424 427 428 435 ...
    445 452 464 465 477 485 513 515 516 520 521 522 540 542 543 544 545 546 547 563 564 566 567 581 592 600 ...
    601:637 639 641 ];
idx_sel = 1:n_soil; idx_sel(idx_notsel) = [];
%idx_not_imag = find(imag(min95(:,3,1))~=0)
if fig_17 == 1
    % Load results from LM_SWIG
    %load SWIG_646_LM_anal_Jac.mat
%    load SWIG_646_LM_test_Jac
    eval(char(load_file));
    load DREAM_SWIG.mat
    %
    y_min = [  0  0    0 ];
    y_max = [ 200 250  2 ];
    dy = 0.04*(y_max - y_min);
    y_min = y_min - dy; y_max = y_max + dy;
    appr = 1;
    % Select only ones that converged
    idx_data = find(converged_DREAM(:,appr) == 1);
    % Create figure name
    fig_name = 'Comparison of DREAM and LM values of infiltration form Haverkamp';
    % Create figure
    %    figure('unit','inches','name',fig_name,'PaperOrientation','portrait','position',[0.3 0.5 18.5 6.5]);
    figure('unit','inches','name',fig_name,'PaperOrientation','portrait','position',[0.3 0.5 20 7.3]);
    % Define index
    index = {'(A)','(B)','(C)','(D)','(E)','(F)','(G)','(H)','(I)','(J)','(K)','(L)','(M)'};
    % Plot 3 figures
    for zz = 1:3
        ax_str = strcat('ax',num2str(zz));
        evalstr = strcat('ax',num2str(zz),{' '},'= axes(''units'',''inches'');'); eval(char(evalstr));
        %  evalstr = strcat('axpos',num2str(zz),{' '},'= [ 1.35 + (zz-1)*5.9 , 1.2 4.5 4.5 ]'); eval(char(evalstr));
        evalstr = strcat('axpos',num2str(zz),{' '},'= [ 1.45 + (zz-1)*6.5 , 1.7 5 5 ];'); eval(char(evalstr));
        evalstr = strcat('set(','ax',num2str(zz),',','''position''',',axpos',num2str(zz),');'); eval(char(evalstr));
        % Add one to one plot
        for i = 1:numel(idx_data)
            u = idx_data(i);
            % horizontal line
            plot3 = line(eval(ax_str),[r_95_DREAM_SWIG(u,1,zz,appr) r_95_DREAM_SWIG(u,4,zz,appr)],...
                [opt_LM_SWIG(u,zz,appr) opt_LM_SWIG(u,zz,appr)],...
                'color',color_conf_line,'linewidth',0.1); if i == 1; hold on; end
            % vertical line
            plot4 = line(eval(ax_str),[opt_DREAM_SWIG(u,zz,appr) opt_DREAM_SWIG(u,zz,appr)],...
                [min95(u,zz,appr) max95(u,zz,appr)],...
                'color',color_conf_line,'linewidth',0.1);
        end
        switch zz
            case 1
                plot1 = plot(eval(ax_str),opt_DREAM_SWIG(idx_data,1,appr),opt_LM_SWIG(idx_data,1,appr),'s','color',color_scatter',...
                    'linewidth',1.5,'markersize',6,'markerfacecolor','w'); hold on
                ylabel(eval(ax_str),'$S_{\rm LM}\;{\rm (cm/h^{1/2})}$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [-0.18, 0.5, 0],'fontweight','normal');
                xlabel(eval(ax_str),'$S_{\rm DREAM}\;{\rm (cm/h^{1/2})}$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [0.5, -0.13, 0],'fontweight','normal');
            case 2
                plot1 = plot(eval(ax_str),opt_DREAM_SWIG(idx_data,2,appr),opt_LM_SWIG(idx_data,2,appr),'s','color',color_scatter',...
                    'linewidth',1.5,'markersize',6,'markerfacecolor','w'); hold on
                ylabel(eval(ax_str),'$K_{\rm s,LM}\;{\rm (cm/h)}$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [-0.18, 0.5, 0],'fontweight','normal');
                xlabel(eval(ax_str),'$K_{\rm s,DREAM}\;{\rm (cm/h)}$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [ 0.5, -0.13, 0],'fontweight','normal');
            case 3
                plot1 = plot(eval(ax_str),opt_DREAM_SWIG(idx_data,3,appr),opt_LM_SWIG(idx_data,3,appr),'s','color',color_scatter',...
                    'linewidth',1.5,'markersize',6,'markerfacecolor','w'); hold on
                ylabel(eval(ax_str),'$\beta_{\rm LM}\;(-)$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [-0.16, 0.5, 0],'fontweight','normal');
                xlabel(eval(ax_str),'$\beta_{\rm DREAM}\;(-)$','interpreter',...
                    'latex','fontsize',fontsize_labels,'Units','normalized','Position',...
                    [ 0.5, -0.13, 0],'fontweight','normal');
        end
        plot2 = line(eval(ax_str),[0 y_max(zz)],[0 y_max(zz)],'color','k','linewidth',1);
        % plot uncertainty of coefficient beta
        axis(eval(ax_str),[y_min(zz) y_max(zz) y_min(zz) y_max(zz) ]);
        % Add (A), (B) and (C)
        set(gca,'fontsize',fontsize_axes);
        % Add labels on x and y-axis
        a = axis;
        set(eval(ax_str),'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 1]); %
        set(eval(ax_str),'TickDir','out','clipping','on');
        % Now determine xticks and xticklabels and minorxticks
        switch zz
            case 3
                xtickformat('%.1f');
            case {1,2}
                xtickformat('%.0f');
        end
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
        if zz == 3
            ytickformat('%.1f')
        end
        % Now create minor ticks - y axis
        %        eval(strcat(ax_str,'.XAxis.MinorTickValues = xminortick'));
        axis(eval(ax_str),[y_min(zz) y_max(zz) y_min(zz) y_max(zz) ]);
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
            [lega,legb,legc,legd] = legend([plot1 plot2 plot4],{'${\rm \;\;SWIG\;data}$',...
                '${\rm \;\;1:1\;Line}$','${\rm \;\;95\%\;intervals}$'},'interpreter','latex','box','off','location',...
                'southeast','fontsize',fontsize_labels);
            set(legb(1),'Color',color_scatter); set(legb(2),'Color','k'); set(legb(3),'Color',color_conf_line);
            set(legb(5),'linewidth',1.5,'markersize',8);
            set(legb(6),'linewidth',2.5);
            set(legb(7),'linewidth',1);
            set(gcf, 'Renderer', 'painters');
        end
        
        % %     set(legb(1),'Color',color_scatter); set(legb(2),'Color','k'); set(legb(3),'Color',color_conf_line);
        % %     set(legb(5),'linewidth',1.5,'markersize',8);
        % %     set(legb(6),'linewidth',2.5);
        % %     set(legb(7),'linewidth',1);
        
        
        axis(eval(ax_str),[y_min(zz) y_max(zz) y_min(zz) y_max(zz) ]);
        a = axis;
        % Add label (A), (B), or (C)
        text(eval(ax_str),a(1)-0.04*(a(2)-a(1)),a(3)+1.05*(a(4)-a(3)),char(index(zz)),...
            'fontsize',fontsize_labels+2,'interpreter','latex');
        set(gcf,'color','w');
        % Infiltration and time form on axis
        % %         if zz == 1
        % %             text(eval(ax_str),-45,-46,'${\rm Infiltration\;form}$','fontsize',fontsize_axes,'interpreter','latex',...
        % %                 'HorizontalAlignment', 'left');
        % %             h = text(eval(ax_str),-50,-37,'${\rm Time\;form}$','fontsize',fontsize_axes,'interpreter','latex',...
        % %                 'HorizontalAlignment', 'left');
        % %             set(h,'rotation',90);
        % %         end
        % %         % inset
        % %         if zz < 3
        % %             ii = find(opt_LM_SWIG(:,zz,1) < 10);
        % %             ax_str = strcat('ax',num2str(zz+5));
        % %             evalstr = strcat('ax',num2str(zz+5),{' '},'= axes(''units'',''inches'');'); eval(char(evalstr));
        % %             %  evalstr = strcat('axpos',num2str(zz),{' '},'= [ 1.35 + (zz-1)*5.9 , 1.2 4.5 4.5 ]'); eval(char(evalstr));
        % %             evalstr = strcat('axpos',num2str(zz+5),{' '},'= [ 1.9 + (zz-1)*6.5 , 4.7 2 2 ];'); eval(char(evalstr));
        % %             evalstr = strcat('set(','ax',num2str(zz+5),',','''position''',',axpos',num2str(zz+5),');'); eval(char(evalstr));
        % %             plot6 = plot(eval(ax_str),opt_LM_SWIG(ii,zz,1),opt_LM_SWIG(ii,zz,2),'s','color',color_scatter,...
        % %                 'linewidth',0.5,'markersize',4,'markerfacecolor','w'); hold on
        % %             plot2 = line(eval(ax_str),[0 10],[0 10],'color','k','linewidth',0.5);
        % %             axis(eval(ax_str),[0 10 0 10]);
        % %             set(gca,'fontsize',12);
        % %             % Add labels on x and y-axis
        % %             a = axis;
        % %             xticklabel = [0:2:10];
        % %             %            xtick = get(eval(ax_str),'xtick');
        % %             xtick = [ 0 2 4 6 8 10 ];
        % %             set(eval(ax_str),'xtick',xtick,'xticklabel',xticklabel);
        % %             xticklabel = get(eval(ax_str),'xticklabel');
        % %             set(eval(ax_str),'xtick',xtick,'xticklabel',[]);
        % %             % remove current labels - and replace with own values closer to axis
        % %             %            set(eval(ax_str),'xticklabel',[]);
        % %             % now plot manually - for all except 12
        % %             dy = a(3) - 0.1*(a(4)-a(3));
        % %             for rr = 1:numel(xtick)
        % %                 text(eval(ax_str),xtick(rr),dy,xticklabel(rr,:),'fontsize',12,'HorizontalAlignment', 'center');
        % %             end
        % %             xticklabel
        % %             set(eval(ax_str),'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 1]); %
        % %             set(eval(ax_str),'TickDir','out','clipping','on');
        % %             set(eval(ax_str),'box','off','color','none')
        % %             % create new, empty axes with box but without ticks
        % %             h3 = axes('Position',get(eval(ax_str),'Position'),'box','on','xtick',[],'ytick',[]);
        % %             linkaxes([eval(ax_str) h3]);
        % %         end
        
    end
    set(legb(1),'Color',color_scatter); set(legb(2),'Color','k'); set(legb(3),'Color',color_conf_line);
    set(legb(5),'linewidth',1.5,'markersize',8);
    set(legb(6),'linewidth',2.5);
    set(legb(7),'linewidth',1);
    % addpath to get nicest pdf from printing to pdf within matlab
    addpath('C:\Jasper Vrugt\Papers\Book\Chapter 1\MATLAB\Figure_resolution');
    % now use export_fig function from Figure_resolution directory
    evalstr = strcat('export_fig figure_DREAM_LM_SWIG',{' '},'-','pdf');
    % now save file
    if save_fig == 1, eval(char(evalstr)); end
    
end

% % % soil_type = [1:n_soil]';
% % % eta1 = opt_LM_SWIG(soil_type,1:3,1);
% % % eta2 = opt_LM_SWIG(soil_type,1:3,2);
% % % A = eta1 - eta2; ii = find(abs(A(:,3))>1);
% % % for i = 1:numel(ii)
% % %     dat = data_SWIG{soil_type(ii(i))};
% % %     plugin.t = dat(:,1); plugin.I = dat(:,2);
% % %     I_sim = Haverkamp_I(eta1(ii(i),1:3),plugin);
% % %     t_sim = Haverkamp_t(eta2(ii(i),1:3),plugin);
% % %     figure(100),plot(plugin.t,I_sim,'r'); hold on
% % %     figure(100),plot(t_sim,plugin.I,'b');
% % %     % Compute limitB-->0 and limit of B --> 2
% % %     Ki = 0; S = eta1(ii(i),1); Ks = eta1(ii(i),2); B = eta1(ii(i),3); 
% % %     t_sim = plugin.t;
% % %     limB0 = -t_sim*(Ks-Ki) -1/2*t_sim*S^2./(I_sim - Ki*t_sim) + 1/2*(I_sim-Ki*t_sim)
% % %     limB2 = ( -(S^2*(exp(-(4*(Ki - Ks)*(I_sim - Ki*t_sim))/S^2)/2 - 1/2))/2 - ...
% % %           exp(-(4*(Ki - Ks)*(I_sim - Ki*t_sim))/S^2)*(Ki - Ks).*(I_sim - Ki*t_sim) - ...
% % %           t_sim*(Ki - Ks)^2.*(exp(-(4*(Ki - Ks)*(I_sim - Ki*t_sim))/S^2) + 1) ) ./ ...
% % %           ( (Ki - Ks)*(exp(-(4*(Ki - Ks)*(I_sim - Ki*t_sim))/S^2) - 1) )
% % %     [limB0 limB2]
% % %     [J,J2,J3] = jac_Haverkamp_anal(eta1(ii(i),1:3),plugin,1);
% % %     figure(i),subplot(3,1,1),plot(I_sim,J(:,1));
% % %     figure(i),subplot(3,1,2),plot(I_sim,J(:,2));
% % %     figure(i),subplot(3,1,3),plot(I_sim,J(:,3));
% % %     pause;
% % % end
% % index = {'(A) Clay','(B) Clay loam','(C) Loam','(D) Loamy sand','(E) Sand','(F) Sandy clay','(G) Sandy clay loam','(H) Sandy loam',...
% %     '(I) Silt','(J) Silt loam','(K) Silty clay','(L) Silty clay loam'};

