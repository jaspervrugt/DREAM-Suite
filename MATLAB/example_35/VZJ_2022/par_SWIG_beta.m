% Runs all samples from SWIG data base
close all hidden; clc; clear all

save_fig = 0;
fontsize_axes = 18;
fontsize_labels = 18;
fig_5 = 1;

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
% Define bounds on parameters
eta_min = [ 1e-3 1e-3 ]'; eta_max = [ inf inf ]';
% Number of parameters
d = 2;
% How many trials with LM from different starting points?
Ntrials = 10; nu = 2; %1.3; % nu = 10;
% Define beta values to use
Bfix = [0.001 0.1:0.1:0.9 0.95 1.05 1.1:0.1:1.9 1.999]; NB = numel(Bfix);
% Initialize return matrices
[opt_LM_SWIG_beta,min95_beta,max95_beta,min90_beta,max90_beta] = deal(nan(n_soil,d,NB,2)); 
[n_LM_SWIG_beta,SSR_LM_SWIG_beta] = deal(nan(n_soil,NB,2));
LM_res_SWIG_beta = nan(n_soil,4,NB,2); warning off

% Now loop over each soil and call DREAM Package
for approach = 1:2
    approach
    parfor uu = 1:NB
        uu
        [A,B,C,D,E,F,G,H] = par_SWIG(Bfix(uu),approach)
        LM_res_SWIG_beta(:,:,uu,approach) = A;
        opt_LM_SWIG_beta(:,:,uu,approach) = B;
        SSR_LM_SWIG_beta(:,uu,approach) = C;
        n_LM_SWIG_beta(:,uu,approach) = D;
        min90_beta(:,:,uu,approach) = E;
        max90_beta(:,:,uu,approach) = F;
        min95_beta(:,:,uu,approach) = G;
        max95_beta(:,:,uu,approach) = H;
    end
end
% Save optimal values
save SWIG_646_beta opt_LM_SWIG_beta SSR_LM_SWIG_beta data_SWIG n_soil Data min95_beta max95_beta min90_beta max90_beta LM_res_SWIG_beta

% Plot results of LM
ct = 0;
for uu = 1:NB
    for par = 1:2
        ct = ct + 1; figure(ct)
        % Plot all data
        plot(opt_LM_SWIG_beta(1:n_soil,par,uu,1),opt_LM_SWIG_beta(1:n_soil,par,uu,2),'ro','markersize',5,'markerfacecolor','w','linewidth',2);
        hold on;
        % Now only low RMSE
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
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 5: PLOT MEASURED VERSUS FITTED CUMULATIVE INFILTRATION DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load SWIG_646_LM_test4_Jac.mat;
idx_sel = 1:n_soil;
if fig_5 == 1
    % Create figure name
    fig_name = 'Plot of observed and fitted data';
    % Define maximum of
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
        fx_LM_beta = []; plugin.t = dat(:,1);
        % Approach 1 - predict I - using minimum beta value (beta = 0.001 is optimal)
        optI = [ opt_LM_SWIG_beta(idx_sel(soil_type),1:2,1,1) Bfix(1) ];
        fx_LM_beta(:,1) = Haverkamp_I(optI,plugin);
        % Approach 2 - predict time - using minimum beta value (beta = 0.3 is optimal)
        plugin.I = dat(:,2);
        optt = [ opt_LM_SWIG_beta(idx_sel(soil_type),1:2,4,2) Bfix(4) ];
        fx_LM_beta(:,2) = Haverkamp_t(optt,plugin);
        % Now with optimum value of eta with beta estimated as well
        fx_LM = []; 
        % Approach 1 - predict I - using minimum beta value (beta = 0.001 is optimal)
        optI = opt_LM_SWIG(idx_sel(soil_type),1:3,1);
        fx_LM(:,1) = Haverkamp_I(optI,plugin);
        % Approach 2 - predict time - using minimum beta value (beta = 0.3 is optimal)
        plugin.I = dat(:,2);
        optt = opt_LM_SWIG(idx_sel(soil_type),1:3,2);
        fx_LM(:,2) = Haverkamp_t(optt,plugin);
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
                    plot(eval(ax_str),dat(:,1),fx_LM_beta(:,1),'b--','linewidth',1.5);
                    plot(eval(ax_str),dat(:,1),fx_LM(:,1),'b','linewidth',1.5);
                case 2 % t optimized LM
                    plot(eval(ax_str),fx_LM_beta(:,2),dat(:,2),'k--','linewidth',1.5);
                    plot(eval(ax_str),fx_LM(:,2),dat(:,2),'k','linewidth',1.5);
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
        ytickformat('%.1f');
        % Add t_char - 2.5 percentile
        set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 1]); %
        set(gca,'TickDir','out');
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
