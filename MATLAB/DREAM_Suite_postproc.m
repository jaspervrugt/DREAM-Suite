function [chain,output,FX,Z] = DREAM_Suite_postproc(method,...
    DREAMPar,f_handle,Par_info,Meas_info,Lik_info,options,chain,output,...
    it,iloc,MAP_info)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%             DDDDD     RRRRR     EEEEEEE     AAA     MM    MM            %
%             DDDDDD    RRRRRR    EEEEEEE    AAAAA    MM    MM            %
%             DD  DD    RR   RR   EE        AA   AA   MMM  MMM            %
%             DD   DD   RR  RR    EEEE      AA   AA   MMMMMMMM            %
%             DD   DD   RRRRR     EEEE      AAAAAAA   MMM  MMM            %
%             DD  DD    RR RR     EE        AAAAAAA   MM    MM            %
%             DDDDDD    RR  RR    EEEEEEE   AA   AA   MM    MM            %
%             DDDDD     RR   RR   EEEEEEE   AA   AA   MM    MM            %
%                                                                         %
%              SSSSSSSS  UU    UU   II   TTTTTTTTTT   EEEEEEE             %
%              SSSSSSS   UU    UU   II   TTTTTTTTTT   EEEEEEE             %
%              SS        UU    UU   II       TT       EE                  %
%              SSSS      UU    UU   II       TT       EEEE                %
%                 SSSS   UU    UU   II       TT       EEEE                %
%                   SS   UU    UU   II       TT       EE                  %
%               SSSSSS   UUUUUUUU   II       TT       EEEEEEE             %
%              SSSSSSS   UUUUUUUU   II       TT       EEEEEEE             %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% This function generates tables and figures of DREAM-Suite               %
%                                                                         %
% SYNOPSIS: [chain,output,FX,Z] = DREAM_Suite_postproc(method,...         %
%               DREAMPar,f_handle,Par_info,Meas_info,Lik_info,...         %
%               options,chain,output,it,iloc,MAP_info)                    %
%                                                                         %
% © Written by Jasper A. Vrugt, Feb 2007                                  %
% Los Alamos National Laboratory 			        	                  %
% University of California Irvine                                         %
%                                                                         %
% Release Version Jan. 2018     Publication of DREAM manual in EMS        %
% Update          Nov. 2021     Changed output and DREAM_Package_end      %
% Update          June 2023     Included new likelihood functions         %
% Update          July 2023     Scoring rules and other additions         %
% Version 2.0     July 2024     Major changes throughout code             %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

n = it - 1;   % --> Remove unnecessary trailing values in output
output.CR = output.CR(1:n,1:DREAMPar.nCR+1);
output.R_stat = output.R_stat(1:n,1:DREAMPar.d+1);
output.MR_stat = output.MR_stat(1:n,1:2);
output.AR = output.AR(1:n,1:2);
chain = chain(1:iloc,1:DREAMPar.d+2,1:DREAMPar.N);
if isfield(MAP_info,'map')
    sndwch = 1;
else
    sndwch = 0;
end

%% SET LABELS AND COLORS FOR OUTPUT
close all hidden; clc;
set(0,'defaultTextInterpreter','latex');
color_order = [
    0.00  0.00  1.00
    0.00  0.50  0.00
    1.00  0.00  0.00
    0.00  0.75  0.75
    0.75  0.00  0.75
    0.75  0.75  0.00
    0.25  0.25  0.25
    0.75  0.25  0.25
    0.95  0.95  0.00
    0.25  0.25  0.75
    0.75  0.75  0.75
    0.00  1.00  0.00
    0.76  0.57  0.17
    0.54  0.63  0.22
    0.34  0.57  0.92
    1.00  0.10  0.60
    0.88  0.75  0.73
    0.10  0.49  0.47
    0.66  0.34  0.65
    0.99  0.41  0.23 ];

fontsize_labels = 20;
fontsize_title = 20;
fontsize_axis_numbers = 16;
fontsize_A = 18;
fontsize_remark = 15;
fontsize_legend = 16;
thickness_line_in_plot = 2;
linewidth_legend = 3;
linewidth_marker = 3;
dark_gray = [0.25 0.25 0.25];
medium_gray = [0.5 0.5 0.5];
light_gray = [0.75 0.75 0.75];
maxbins = 25;

%% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%%                          DEFAULT VARIABLES                         
%% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

p_alfa = 0.05;          % significance level
alfa = 100*(1-p_alfa);  % confidence interval (%)

%% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%%                          TABULATED OUTPUT                         
%% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

%% 1. OPEN WARNING FILE
fid_w = fopen('warning_file.txt','a+');
% Print whether there is a sandwich correction or not 
switch sndwch 
    case 0
        fprintf(fid_w,'Sandwich correction is INACTIVE\n');
    case 1
        fprintf(fid_w,'Sandwich correction is ACTIVE\n');
end

%% 2A. READ SIMULATED OUTPUT
switch options.modout
    case 'no'
        FX = [];
    case 'yes'
        if ( Meas_info.n > 0 ) || ( Meas_info.n_S > 0 )
            fid_FX = fopen('FX.bin','r','n');
            FX = fread(fid_FX, [ Meas_info.n + Meas_info.n_S , ...
                DREAMPar.N * iloc ],'double')';
            fclose(fid_FX);
        else
            FX = [];
        end
    otherwise
end

%% 2B. READ PAST ARCHIVE OF SAMPLES
switch method
    case {'dream','dream_d'}
        Z = [];
    case {'dream_zs','dream_dzs','dream_kzs','mtdream_zs'}
        fid_Z = fopen('Z.bin','r','n');
        Z = fread(fid_Z, [DREAMPar.d+2,DREAMPar.m],'double')';
        fclose(fid_Z);
end

method = upper(method);
%% 3. TABULATE WITHIN-CHAIN CONVERGENCE DIAGNOSTICS
if strcmp(options.diagnostics,'yes')
    if iloc > 200 % At least 200 samples must be in each chain to run CODA
        tabulate_diagnostics(method,DREAMPar,chain,iloc,sndwch);
    else % Error -- not enough samples in each chain!
        evalstr = char(strcat(['DREAM_Suite WARNING: Cannot ' ...
            'compute coda diagnostics for chains due to inufficient ' ...
            'chain samples -> Use at least DREAMPar.T = '],{' '},...
            num2str((200 * DREAMPar.thinning)),'\n'));
        fprintf(evalstr); fprintf(fid_w,evalstr);
    end
end

%% 4. CHECK WHETHER CHAIN HAS FORMALLY CONVERGED (MULTI-CHAIN STATISTIC)
if ~( sum(output.R_stat(end,2:DREAMPar.d+1) < 1.2 ) == DREAMPar.d )
    evalstr = char(strcat(['DREAM_Suite WARNING: Chains did not ' ...
        'converge according to \\hat{R} scale reduction factor\n']));
    fprintf(evalstr); fprintf(fid_w,evalstr);
end

%% 5. WRITE FINAL LINES OF WARNING_FILE.TXT
fprintf(fid_w,'----------- End of DREAM_Suite warning file ---------\n');
fclose(fid_w);
if ispc || ismac, edit warning_file.txt, end

%% 6. TABULATE POSTERIOR MOMENTS
P = genparset(chain);       % First assemble all chains in one matrix
nP = size(P,1);             % # parameter samples
% Take the last 25% of the posterior samples -- assume that these samples
% are posterior samples (double check that R_stat < 1.2 for all parameters)
P_post = P(floor(3/4*nP):nP,1:DREAMPar.d+2);
if ~ismember(DREAMPar.lik,23) % --> not limits of acceptability
    [~,id] = max(P_post(:,DREAMPar.d+2)); id_ML = id(1);
    ML = P_post(id_ML,1:DREAMPar.d);
    [~,id] = max(sum(P_post(:,DREAMPar.d+1:DREAMPar.d+2),2)); id_MAP = id(1);
    MAP = P_post(id_MAP,1:DREAMPar.d);
    MAP_pr = P_post(id_MAP,DREAMPar.d+1);
    MAP_lik = P_post(id_MAP,DREAMPar.d+2);
else
    ML = nan(1,DREAMPar.d); MAP = nan(1,DREAMPar.d); 
            [MAP_pr,MAP_lik] = deal(nan);
end
MEAN = mean(P_post(:,1:DREAMPar.d));        % Mean posterior parameter values
STD = std(P_post(:,1:DREAMPar.d));          % Posterior parameter standard deviation
CORR = corr(P_post(:,1:DREAMPar.d),...      % Posterior correlation matrix
    'type','Pearson');                      
tabulate_output(method,DREAMPar,ML,...      % Tabulate output
    MAP,MEAN,STD,CORR,sndwch);  

%% <><><><><><><><><><> END OF TABULATED OUTPUT ><><><><><><><><><><><><><>

%% <><><><><><><><><><><><><> FIGURE OUTPUT ><><><><><><><><><><><><><><><>
if strcmp(options.print,'yes') % continue with script below
else % otherwise return back to function header
    return
end

switch sndwch
    case 0
        sndwch_text = ': Sandwich Inactive';
    case 1
        sndwch_text = ': Sandwich Active';
end        

% Check whether model produces simulation or not?
sim_out = []; S_post = []; FX_post = [];
% What are the possibilities:
% model output is likelihood, no Meas_info.Y, or Meas_info.S
% model output is log-likelihood, no Meas_info.Y or Meas_info.S
% model output are summary metrics, Meas_info.S (ABC)
% model output is simulation, Meas_info.Y, and/or Meas_info.S

% First create or unpack the model output simulations of the posterior samples
if ( Meas_info.n > 0 ) || ( Meas_info.n_S > 0 )
    tol = 1e-8;
    switch strcmp(options.modout,'yes')         % Did store output?
        case 1      % If yes, then take from FX
%            sim_out = FX(floor(3/4*size(FX,1)):size(FX,1),1:end);
            sim_out = FX(floor(3/4*nP):nP,1:end);
        otherwise   % Generate posterior simulations
            if DREAMPar.d > 1
                [UP_post,iiUP,idUP] = uniquetol(P_post(:,1:DREAMPar.d),...
                    tol,'byrows',true);
                % Old method: Duplicate can still exist: > 10 places behind
                % comma = writing/storage mistake
                % [UP_post,iiUP,idUP] = unique(P_post(:,1:DREAMPar.d),'rows');               
            else
                [UP_post,iiUP,idUP] = uniquetol(P_post(:,1:DREAMPar.d),...
                    tol);
                % Old method: Duplicate can still exist: > 10 places behind
                % comma = writing/storage mistake
                % [UP_post,iiUP,idUP] = unique(P_post(:,1:DREAMPar.d),...
                %      'stable');
            end
            % Now transform to unnormalized space, if necessary
            UP_postun = X_unnormalize(UP_post,Par_info);
            % Evaluate the model for B_un with verbose = 1
            [sim_out,S_out] = Evaluate_target(UP_postun,DREAMPar, ...
                Meas_info,options,f_handle,1);
            % Duplicate simulations of unique parameter vectors
            sim_out = [ sim_out , S_out ]'; sim_out = sim_out(idUP,:);
            % Clear unused variables
            clear UP_post
    end
end

% Now we need to determine what sim_out actually is
switch strcmp(options.DB,'yes')
    case 1 % sim_out contains simulation and summary metrics
        FX_post = sim_out(:,1:Meas_info.n);
        S_post = sim_out(:,Meas_info.n+1:Meas_info.n+Meas_info.n_S);
    otherwise % Summary metrics or LOA
        if ismember(DREAMPar.lik,21:23)
            S_post = sim_out;
        else
            FX_post = sim_out;
        end
end

if ~isempty(FX_post) % Now get the corresponding model simulation
    FX_MAP = FX_post(id_MAP,1:Meas_info.n)';
    % Compute the RMSE of the maximum aposteriori solution
    RMSE_MAP = sqrt ( sum ( ( FX_MAP - Meas_info.Y(:)).^2) / Meas_info.n );
    % Compute the BIAS of the maximum aposteriori solution
    PBIAS_MAP = 100 * sum ( FX_MAP - Meas_info.Y(:) ) / sum ( Meas_info.Y );
else
    RMSE_MAP = []; PBIAS_MAP = [];
end
% display RMSE and PBIAS
disp(RMSE_MAP), disp(PBIAS_MAP);
% rename method so that subscript is correctly printed
method_old = method;
ii = strfind(method,'_');
if ~isempty(ii)
    method = strcat(method(1:ii-1),'$_{\rm (',method(ii+1:end),')}$');
end
% Create legend/label string for different parameters
str_par = cell(DREAMPar.d,1);
switch isfield(Par_info,'names')
    case 1 % User specified the names of the parameters
        for i = 1:DREAMPar.d
            str_par(i,:) = cellstr(strcat('$\;',Par_info.names(i),'$'));
        end
    otherwise % User did not specify names: use x_1, x_2, etc. instead
        for i = 1:DREAMPar.d
            str_par(i,:) = cellstr(strcat('$\;x_{',num2str(i),'}$'));
        end
end
% How many chains are plotted in relevant figures?
M = min(DREAMPar.N,5);
% Create legend/label string for different chains
str_chain = cell(M,1);
for i = 1:M, str_chain(i,:) = cellstr(strcat('$\;{\rm chain}','\;', ...
        num2str(i),'$')); end
% Create legend/label string for different crossover values
str_cr = cell(DREAMPar.nCR,1);
for i = 1:DREAMPar.nCR, str_cr(i,:) = strcat('$\;n_{\rm CR} = ',{' '},...
        num2str(i/DREAMPar.nCR,'%4.2f'),'$'); end
% Get the screen size
s0 = get(0,'screensize');

% Extract t_start from Lik_info
t_start = Lik_info.t_start;

%% <><><><><><><><><><><><> END OF PRE-PROCESSING ><><><><><><><><><><><><>

%% <><> EVOLUTION OF UNI/MULTIVARIATE R_STATISTIC OF GELMAN AND RUBIN ><><>
figure('units','normalized','position',[0 0 1 1]);
% Define first axis
ax1 = axes('units','normalized'); axpos1 = [ 0.1 0.6 0.8 0.35 ]; set(ax1,...
    'position',axpos1);
% First print the R-statistic of Gelman and Rubin (different colors)
semilogy(ax1,output.R_stat(:,1)/DREAMPar.N,output.R_stat(:,2:DREAMPar.d+1),...
    'linewidth',thickness_line_in_plot); hold on;
ax1.ColorOrder = color_order;
% Now add the theoretical convergence value of 1.2 as horizontal line
semilogy(ax1,[0 DREAMPar.T],[1.2 1.2],'k--','linewidth', ...
    thickness_line_in_plot);
% Define axes
r_max = max(max(output.R_stat(:,2:DREAMPar.d+1))); yl = ylim;
y_min = 0.80; y_max = 10^ceil(log10(min(r_max,yl(2))));
axis(ax1,[ 0 DREAMPar.T y_min y_max]); xylim = [ xlim ylim ];
% Add title
text_title = strcat(method,':',{' '},...
    'Convergence of sampled chains - scale reduction factors',sndwch_text);
text(ax1,(xylim(1)+xylim(2))/2,xylim(3) + 10.^(log10(xylim(3)) + ...
    1.04*(log10(xylim(4))-log10(xylim(3)))),text_title,'fontsize',...
    fontsize_title,'interpreter','latex',...
    'HorizontalAlignment','center');
% Define x-axis ticks and tick direction
set(ax1,'fontsize',fontsize_axis_numbers,'TickDir','out','XMinorTick',...
    'on','YMinorTick','on');
% Set values
xaxis_numbers = 0:floor(DREAMPar.T/5):DREAMPar.T;
set(ax1,'xtick',xaxis_numbers,'xticklabels',[]);
% Add a legend
if DREAMPar.d <= 13
    dx = xylim(2) - xylim(1); %dy = xylim(4) - xylim(3);
    x_loc = xylim(1) + [ 0.94 0.97 ] * dx;
    % Plot legend manually
    for z = 1:DREAMPar.d
        y_loc = 10.^(log10(xylim(3)) + ([0.95 0.95] - ...
            (z-1)/15)*(log10(xylim(4))-log10(xylim(3))));
        line(ax1,x_loc,y_loc,'color',color_order(z,1:3),'linewidth',...
            linewidth_legend);
        text(ax1,x_loc(2)+0.01,y_loc(1),char(str_par(z)),...
            'interpreter','latex','fontsize',fontsize_legend,...
            'color',color_order(z,1:3),'horizontalalignment','left');
    end
end
% Check when convergence has been achieved
check_conv = sum(output.R_stat(:,2:DREAMPar.d+1) > 1.2 , 2);
if any(check_conv > 0) % at least one R_stat has been larger than 1.2
    idx_conv = find(check_conv > 0,1,'last');
    if idx_conv == numel(check_conv)
        idx_conv = [];
    else
        idx_conv = idx_conv + 1;
    end
    % --> idx_conv can be empty --> code has not yet converged
else % R_stat values have always been smaller than 1.2
    idx_conv = 2;
end
switch ~isempty(idx_conv)
    case 1 % did converge on time
        line(ax1,[output.R_stat(idx_conv,1)/DREAMPar.N ...
            output.R_stat(idx_conv,1)/DREAMPar.N],...
            [y_min y_max],'linewidth',thickness_line_in_plot,...
            'color',[0.5 0.5 0.5],'linestyle','--');
        x_loc = output.R_stat(idx_conv,1)/DREAMPar.N + ...
            0.005*output.R_stat(end,1)/DREAMPar.N;
        y_loc = 10.^(log10(xylim(3)) + 0.5*(log10(xylim(4)) - ...
            log10(xylim(3))));
        text(ax1,x_loc,y_loc,'converged','fontsize',fontsize_remark,...
            'rotation',90,'horizontalalignment','center');
    otherwise % did not converge on time --> no plotting
end
ylh = ylabel(ax1,'$\widehat{R}\,-\,{\rm diagnostic}$',...
    'fontsize',fontsize_labels);
ylh.Position(1) = xylim(1) - 0.05 * (xylim(2) - xylim(1));
text(ax1,0.005,0.94,...
    '(A) Univariate $\widehat{R}$\,-\,diagnostic','fontsize',...
    fontsize_A,'units','normalized');
% Plot box around figure;
plot_box(ax1);

% Define first axis
ax2 = axes('units','normalized'); axpos2 = [ 0.1 0.15 0.8 0.35 ]; 
set(ax2,'position',axpos2);
% First print the R-statistic of Gelman and Rubin (different colors)
semilogy(ax2,output.MR_stat(:,1)/DREAMPar.N,output.MR_stat(:,2),...
    'linewidth',thickness_line_in_plot); hold on;
ax2.ColorOrder = color_order;
% Set the the axes
r_max = max(max(output.MR_stat(:,2)));
yl = ylim; %a_m = max(10,1.1*a(4));
y_max = 10^ceil(log10(min(r_max,yl(2))));
axis(ax2,[ 0 DREAMPar.T 0.80 y_max]); xylim = [ xlim ylim ];
% Now add the theoretical convergence value of 1.2 as horizontal line
plot(ax2,[0 DREAMPar.T],[1.2 1.2],'k--','linewidth', ...
    thickness_line_in_plot);
% Set the the axes
set(ax2,'fontsize',fontsize_axis_numbers,'tickdir','out',...
    'XMinorTick','on','YMinorTick','on');
% Add comma's for large numbers
set(ax2,'xtick',xaxis_numbers,'xticklabels',addcomma(xaxis_numbers));
% Check when convergence has been achieved
check_conv = sum(output.MR_stat(:,2) > 1.2 , 2);
if any(check_conv > 0) % at least one R_stat has been larger than 1.2
    idx_conv = find(check_conv > 0,1,'last');
    if idx_conv == numel(check_conv)
        idx_conv = [];
    else
        idx_conv = idx_conv + 1;
    end
    % --> idx_conv can be empty --> code has not yet converged
else % MR_stat values have always been smaller than 1.2
    idx_conv = 2;
end

switch ~isempty(idx_conv)
    case 1 % did converge on time
        line(ax2,[output.MR_stat(idx_conv,1)/DREAMPar.N ...
            output.MR_stat(idx_conv,1)/DREAMPar.N],...
            [y_min y_max],'linewidth',thickness_line_in_plot,'color',...
            [0.5 0.5 0.5],'linestyle','--');
        x_loc = output.MR_stat(idx_conv,1)/DREAMPar.N + ...
            0.005*output.MR_stat(end,1)/DREAMPar.N;
        y_loc = 10.^(log10(xylim(3)) + 0.5*(log10(xylim(4)) - ...
            log10(xylim(3))));
        text(ax2,x_loc,y_loc,'converged','fontsize',fontsize_remark,...
            'rotation',90,'horizontalalignment','center');
    otherwise % did not converge on time --> no plotting
end
% Add labels
xlh = xlabel(ax2,'${\rm Number\;of\;generations}$',...
    'fontsize',fontsize_labels);
xlh.Position(2) = 10.^(log10(xylim(3)) - 0.17*(log10(xylim(4)) - ...
    log10(xylim(3))));
ylh = ylabel(ax2,'$\widehat{R}^{d}\,-\,{\rm diagnostic}$',...
    'fontsize',fontsize_labels);
ylh.Position(1) = xylim(1) - 0.05 * (xylim(2) - xylim(1));
% Assign label to plot
text(ax2,0.005,0.94,...
    '(B) Multivariate $\widehat{R}^{d}$\,-\,diagnostic',...
    'fontsize',fontsize_A,'units','normalized');
% Plot box around figure;
plot_box(ax2); set(gcf,'color','w');

%% <><><><> EVOLUTION OF ACCEPTANCE RATE AND PCR PROBABILITIES <><><><><><>
figure('units','normalized','outerposition',[0 0 1 1]);
% Define first axis
ax3 = axes('units','normalized'); axpos3 = [ 0.1 0.6 0.8 0.35 ]; set(ax3,...
    'position',axpos3);
% Plot the acceptance rate
plot(ax3,output.AR(:,1)/DREAMPar.N,output.AR(:,2),...
    'linewidth',thickness_line_in_plot,'color',[0.5 0.5 0.5]); hold on;
ax3.ColorOrder = color_order;
% Increase fontsize axis numbers
set(ax3,'fontsize',fontsize_axis_numbers,'tickdir','out',...
    'XMinorTick','on','YMinorTick','on');
% Add title
title(ax3,strcat(method,':',{' '},'Evolution of acceptance rate', ...
    sndwch_text),'fontsize',fontsize_title);
% Now add the theoretical convergence value of 1.2 as horizontal line
yl = ylim; axis(ax3,[0 DREAMPar.T 0 yl(2)]);
% define axis values
xaxis_numbers = 0:floor(DREAMPar.T/5):DREAMPar.T;
% Add comma's for large numbers
set(ax3,'xtick',xaxis_numbers,'xticklabels',[]);
% And the significant digits on x-axis numbers
fix_ticklabels(ax3,'y',1); % Add one digit to acceptance rate
% Add labels
ylh = ylabel(ax3,'${\rm Acceptance\;rate}\;(\%)$','fontsize',...
    fontsize_labels,'interpreter','latex');
xylim = [xlim ylim];
ylh.Position(1) = xylim(1) - 0.05 * (xylim(2) - xylim(1));
% Add label(A) on the axis
text(ax3,0.005,0.94,'(A)','fontsize',fontsize_A,'units','normalized');
% Plot box around figure;
plot_box(ax3);

% Define second axis
ax4 = axes('units','normalized'); axpos4 = [ 0.1 0.15 0.8 0.35 ]; 
set(ax4,'position',axpos4);
% First print the pCR values
plot(ax4,output.CR(:,1)/DREAMPar.N,output.CR(:,2:end),...
    'linewidth',thickness_line_in_plot);
set(ax4,'fontsize',fontsize_axis_numbers,'tickdir','out',...
    'XMinorTick','on','YMinorTick','on');
% Add comma's for large numbers
set(ax4,'xtick',xaxis_numbers,'xticklabels',addcomma(xaxis_numbers));
% Adjust yaxis so that we can plot legend
xylim = [xlim ylim];
y_min = xylim(3);
y_max = ceil(10*(xylim(4) + 0.15*(xylim(4) - xylim(3))))/10;
axis(ax4,[0 , DREAMPar.T , y_min y_max]); xylim = [xlim ylim];
% Set colormap
ax4.ColorOrder = color_order;
% Define xlabel
xlh = xlabel(ax4,'${\rm Number\;of\;generations}$','fontsize',...
    fontsize_labels,'interpreter','latex');
xlh.Position(2) = xylim(3) - 0.17*(xylim(4) - xylim(3));
% Make the number of significant digits consistent on y-axis
fix_ticklabels(ax4,'y');
% Add ylabel
ylh = ylabel(ax4,'${\rm Selection\;probability}\;(-)$','fontsize',...
    fontsize_labels,'interpreter','latex');
xylim = [xlim ylim];
ylh.Position(1) = xylim(1) - 0.05 * (xylim(2) - xylim(1));

% Add title
title_plot = strcat(method,':',{' '},...
    'Evolution of selection probabilities crossover values');
if strcmp(DREAMPar.adapt_pCR,'no')
    title_plot = strcat(title_plot,': No adaptation');
end
title(ax4,title_plot,'fontsize',fontsize_title);
% Plot legend manually
dx = xylim(2) - xylim(1); dy = xylim(4) - xylim(3);
x_loc = xylim(1) + [ 0.89 0.91 ] * dx;
y_loc0 = xylim(3) + [ 0.93 0.93 ] * dy;
for i = 1:DREAMPar.nCR
    y_loc = y_loc0 - (i-1)/12 * dy;
    line(ax4,x_loc,y_loc,'color',color_order(i,1:3),'linewidth',...
        linewidth_legend);
    text(ax4,x_loc(2)+0.005,y_loc(1),char(str_cr(i)),...
        'interpreter','latex','fontsize',fontsize_legend,'color',...
        color_order(i,1:3),'horizontalalignment','left');
end
% Add label (B) to plot
text(ax4,0.005,0.94,'(B)','fontsize',fontsize_A,'units','normalized');
% Plot box around figure;
plot_box(ax4); set(gcf,'color','w');

%% <><><>< AUTOCORRELATION PLOTS OF THE SAMPLED PARAMETER SAMPLES ><><><><>
n_row = 4; n_col = 2;
% Now determine row and column number of each parameter
[row_par,col_par,idx_fig] = deal(nan(DREAMPar.d,1));
[row_par(1),col_par(1),idx_fig(1)] = deal(1);
for d = 2:DREAMPar.d
    row_par(d) = row_par(d-1) + 1; col_par(d) = col_par(d-1);
    if row_par(d) == n_row + 1
        row_par(d) = 1; col_par(d) = col_par(d) + 1;
    end
    if col_par(d) == n_col + 1
        col_par(d) = 1; idx_fig(d) = 1;
    else
        idx_fig(d) = 0;
    end
end
T = size(chain,1);                          %
maxlag = min(150,T);                        % Now determine the maximum lag
xaxis_numbers = 0:floor(maxlag/5):maxlag;   % define axis values
% Add a legend
M = min(5,DREAMPar.N);
% Create title
title_str = strcat(char(strcat(method,':',{' '},...
    'Sample autocorrelation function of parameters')),sndwch_text); 
clear ax_acf;
% % % fig_num = 2;
% Must start at figure 3
for par = 1:DREAMPar.d
    if idx_fig(par) == 1            % Open new figure/add title to previous one
        if exist('ax_acf','var')    % Add title to previous figure
            mtit(char(title_str),'fontsize',fontsize_title,...
                'interpreter','latex');
        end
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    % Define first axis
    ax_acf = axes('units','normalized'); axpos_acf = [ 0.075 + ...
        (col_par(par)-1) * 0.48 , ...
        0.78 - (row_par(par)-1) * 0.22 , 0.42 , 0.17 ];
    set(ax_acf,'position',axpos_acf);
    % Plot the ACF of each parameter
    acf_chain = nan(maxlag+1,min(5,DREAMPar.N));
    for z = 1:M
        acf_chain(:,z) = acf(chain(1:T,par,z),maxlag);
        plot(ax_acf,0:maxlag,acf_chain(:,z),'color',color_order(z,1:3),...
            'linewidth',thickness_line_in_plot); if z == 1; hold on; end
    end
    y_min = min(acf_chain(:)); %y_max = max(acf_chain(:));
    % Plot the ACF
    ax_acf.ColorOrder = color_order;
    if y_min < 0
        y_min = 1.1*y_min;
    else
        y_min = 0.9*y_min;
    end
    axis(ax_acf,[0 maxlag y_min 1.02]);
    set(ax_acf,'fontsize',fontsize_axis_numbers,'tickdir','out');
    ytickformat('%2.1f');
    % Add legend first subplot
    switch idx_fig(par)
        case 1
            xylim = [ get(ax_acf,'XLim') get(ax_acf,'YLim') ];
            dx = (xylim(2)-xylim(1)); dy = (xylim(4)-xylim(3));
            x_loc = xylim(1) + [ 0.85 0.89 ] * dx;
            y_loc0 = xylim(3) + [ 0.9 0.9 ] * dy;
            % Plot legend manually
            for z = 1:M
                y_loc = y_loc0 - (z-1)/8 * dy;
                line(ax_acf,x_loc,y_loc,'color',color_order(z,1:3),...
                    'linewidth',linewidth_legend);
                text(ax_acf,x_loc(2)+0.005,y_loc(1),char(str_chain(z)),...
                    'interpreter','latex','fontsize',...
                    fontsize_legend,'color',color_order(z,1:3),...
                    'horizontalalignment','left');
            end
    end
    switch col_par(par)
        case 1  % Add y-label
            ylh = ylabel(ax_acf,'Corr. coeff.','fontsize',fontsize_labels);
            xylim = [xlim ylim];
            ylh.Position(1) = xylim(1) - 0.065 * (xylim(2) - xylim(1));
    end
    if row_par(par) == n_row || (par == DREAMPar.d)
        % Add comma's for large numbers
        set(ax_acf,'xtick',xaxis_numbers,'xticklabels',xaxis_numbers);
        xlabel(ax_acf,'Lag','fontsize',fontsize_labels);
    else
        set(ax_acf,'xtick',xaxis_numbers,'xticklabels',[]);
    end
    fig_code = strcat('(',ranktoletter(par),')');
    text(ax_acf,0.01,0.12,strcat(fig_code,{' '},char(str_par(par))),...
        'fontsize',fontsize_A,'units','normalized','interpreter', ...
        'latex','horizontalalignment','left');
    % Add minor ticks
    ax_acf.XMinorTick = 'on'; ax_acf.YMinorTick = 'on';
    % If DREAMPar.d then add title
    if par == DREAMPar.d
        mtit(char(title_str),'fontsize',fontsize_title,...
            'interpreter','latex');
    end
    % Plot box around figure
    plot_box(ax_acf); set(gcf,'color','w');
end

%% <><><><><><>< CONVERGENCE CHAINS TO TARGET DISTRIBUTION <><><><><><><><>
n_row = 3;
% Now determine row and column number of each parameter
[row_par,idx_fig] = deal(nan(DREAMPar.d,1));
[row_par(1),idx_fig(1)] = deal(1);
for d = 2:DREAMPar.d
    row_par(d) = row_par(d-1) + 1;
    if row_par(d) == n_row + 1
        row_par(d) = 1; idx_fig(d) = 1;
    else
        idx_fig(d) = 0;
    end
end
xaxis_numbers = [1 DREAMPar.T/5 2*DREAMPar.T/5:DREAMPar.T/5:DREAMPar.T];
% Now set x-axis
switch DREAMPar.thinning
    case 1
        xaxis_values = (1:DREAMPar.T)';
    otherwise
        xaxis_values = [1 DREAMPar.thinning:DREAMPar.thinning:DREAMPar.T]';
end
% Add legend entry
M = min(DREAMPar.N,5);
% Now loop over each parameter
for par = 1:DREAMPar.d
    if idx_fig(par) == 1 % Open new figure/add title
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    % Define axis
    ax_conv = axes('units','normalized'); axpos_conv = [ 0.06 , ...
        0.7 - (row_par(par)-1) * 0.28 , 0.9 , 0.25 ]; height_plot = 0.25;
    set(ax_conv,'position',axpos_conv);

    % Now plot a number of chains
    for i = 1:M
        plot(ax_conv,xaxis_values,chain(1:end,par,i),'color',...
            color_order(i,1:3),'linewidth',thickness_line_in_plot);
        if i == 1; hold on; end
    end
    ax_conv.ColorOrder = color_order;
    set(ax_conv,'fontsize',fontsize_axis_numbers,'tickdir','out',...
        'XMinorTick','on','YMinorTick','on');
    if isfield(Par_info,'min') && isfield(Par_info,'boundhandling')
        %if sum(strncmp(Par_info.boundhandling,{'fold','reflect',...
        %          'bound'},inf)) --> gives issue with prior M2016 releases
        if strcmp(Par_info.boundhandling,'fold') || ...
                strcmp(Par_info.boundhandling,'reflect') ...
                || strcmp(Par_info.boundhandling,'bound')
            % Use scaling with prior parameter ranges
            axis(ax_conv,[1 DREAMPar.T Par_info.min(par) ...
                Par_info.max(par)]);
        end
    else
        a = axis; axis(ax_conv,[1 DREAMPar.T min(0.9*a(3),1.1*a(3)) ...
            max(0.9*a(4),1.1*a(4)) ]);
    end
    plot(ax_conv, DREAMPar.T , MAP(par),'kx',... %'color',color_order(end,1:3),...
        'Markersize',fontsize_remark,'linewidth',linewidth_marker);
    ylh = ylabel(ax_conv,str_par(par),'fontsize',fontsize_labels,...
        'interpreter','latex');
    xylim = [xlim ylim];
    ylh.Position(1) = xylim(1) - 0.04 * (xylim(2) - xylim(1));
    if row_par(par) == 1    % Then add title and legend
        switch par
            case DREAMPar.d - 1
                evalstr_t = strcat(method,':',{' '},...
                    'Chain convergence plot of parameters',...
                    str_par(par),',',str_par(par+1));
            case DREAMPar.d
                evalstr_t = strcat(method,':',{' '},...
                    'Chain convergence plot of parameter',...
                    str_par(par));
            otherwise
                evalstr_t = strcat(method,':',{' '},...
                    'Chain convergence plot of parameters',...
                    str_par(par),',',str_par(par+1),',',str_par(par+2));
        end
        % Add title
        title(ax_conv,strcat(evalstr_t,sndwch_text), ...
                'fontsize',fontsize_title,'interpreter','latex');
        % Now add legend - manually to match colors to labels
        xylim = [ get(ax_conv,'XLim') get(ax_conv,'YLim') ];
        dx = (xylim(2)-xylim(1)); dy = (xylim(4)-xylim(3));
        x_loc0 = xylim(1) + [ 0.10 0.12 ] * dx;
        y_loc = xylim(3) + [ 0.9 0.9 ] * dy;
        % Plot legend manually
        for z = 1:M
            x_loc = x_loc0 + (z-1)/12 * dx;
            line(ax_conv,x_loc,y_loc,'color',color_order(z,1:3),'linewidth',...
                linewidth_legend);
            text(ax_conv,x_loc(2)+0.01,y_loc(1),char(str_chain(z)),...
                'interpreter','latex','fontsize',fontsize_legend,...
                'color',color_order(z,1:3),'horizontalalignment','left');
        end
        x_loc = x_loc0 + M/12 * dx;
        % Now add MAP
        hmult = 400*height_plot;
        plot(ax_conv,x_loc(1),y_loc(1),'kx','markersize',10,'linewidth',2);
        text(ax_conv,x_loc(1) + 0.005*dx,y_loc(1)-dy/hmult,'$\;{\rm MAP}$',...
            'interpreter','latex','fontsize',fontsize_legend,'color','k',...
            'horizontalalignment','left','verticalalignment','middle');
    end
    if row_par(par) == 3 || (par == DREAMPar.d)   % Add xlabel and x-values
        set(ax_conv,'xtick',xaxis_numbers,'xticklabels', ...
            addcomma(xaxis_numbers));
        % Add a title
        xlh = xlabel(ax_conv,'${\rm Sample\;number\;of\;chain}$',...
            'fontsize',fontsize_labels,'interpreter','latex');
        xylim = [xlim ylim];
        xlh.Position(2) = xylim(3) - 0.24*(xylim(4) - xylim(3));
    else
        set(ax_conv,'xtick',xaxis_numbers,'xticklabels',[]);
    end
    % Check ticklabels on y-axis
    fix_ticklabels(ax_conv,'y');
    % Add label to each figure
    fig_code = strcat('(',ranktoletter(par),')');
    text(ax_conv,0.003,0.9,fig_code,'units','normalized',...
        'fontsize',fontsize_A,'interpreter','latex',...
        'horizontalalignment','left'); box off;
    % Plot box around figure;
    plot_box(ax_conv); set(gcf,'color','w'); 
end

%% <><><> CONVERGENCE OF LIKELIHOOD AND PRIOR OF TARGET DISTRIBUTION <><><>
figure('units','normalized','outerposition',[0 0 1 1]);
% Now set x-axis
switch DREAMPar.thinning
    case 1
        xaxis_values = (1:DREAMPar.T)';
    otherwise
        xaxis_values = [1 DREAMPar.thinning:DREAMPar.thinning:DREAMPar.T]';
end
xaxis_numbers = [1 DREAMPar.T/5 2*DREAMPar.T/5:DREAMPar.T/5:DREAMPar.T];
M = min(DREAMPar.N,5);
% Now loop
for id_plot = 1:2
    % Define first axis
    switch id_plot
        case 1  % First axes
            ax_gen = axes('units','normalized'); ax_genpos = ...
                [ 0.1 0.6 0.8 0.35 ];
            set(ax_gen,'position',ax_genpos);
        case 2  % Define second axis
            ax_gen = axes('units','normalized'); ax_genpos = ...
                [ 0.1 0.15 0.8 0.35 ];
            set(ax_gen,'position',ax_genpos);
    end
    ax_gen.ColorOrder = color_order; height_plot = 0.35;
    mnmx_chain = nan(M,2);
    % Now plot a number of chains
    for i = 1:M
        switch id_plot
            case 1
                plot_chain = chain(1:end,DREAMPar.d+2,i);
            case 2
                plot_chain = chain(1:end,DREAMPar.d+1,i) + ...
                    chain(1:end,DREAMPar.d+2,i);
        end
        plot(ax_gen,xaxis_values,plot_chain,...
            'color',color_order(i,1:3),'linewidth',thickness_line_in_plot);
        mnmx_chain(i,1:2) = [ min(plot_chain) max(plot_chain) ];
        if i == 1; hold on; end
    end
    % Set font
    set(ax_gen,'fontsize',fontsize_axis_numbers,'tickdir','out',...
        'XMinorTick','on','YMinorTick','on');
    switch id_plot
        case 1
            set(ax_gen,'xtick',xaxis_numbers,'xticklabels',[]);
        case 2
            set(ax_gen,'xtick',xaxis_numbers,'xticklabels',...
                addcomma(xaxis_numbers));
    end
    % Ranges have not been defined -- need to derive them from ParSet
    axis tight; a = axis; da = a(4) - a(3); a_m = a(3) - da/5; a_p = a(4) + da/5;
    axis(ax_gen,[1 DREAMPar.T a_m a_p])
    % Lets add the MAP value
    switch id_plot
        case 1
            plot(ax_gen, DREAMPar.T , MAP_lik,'kx','Markersize',...
                fontsize_remark,'linewidth',3);
            % Then add a y-label
            if ~ismember(DREAMPar.lik,[22 23])
                ylabel(ax_gen,['$\mathcal{L}(\textbf{x}|' ...
                    '\widetilde{\textbf{Y}})$'],...
                    'interpreter','latex','fontsize',fontsize_labels);
                % Then add title
                evalstr_t = strcat(method,':',{' '},...
                    ['Chain convergence plot of log-likelihood, ' ...
                    '$\mathcal{L}(\textbf{x}|\widetilde{\textbf{Y}})$']);
                title(ax_gen,strcat(evalstr_t,sndwch_text),'fontsize', ...
                    fontsize_title,...
                    'interpreter','latex');
            elseif ismember(DREAMPar.lik,22)
                ylabel(ax_gen,'$\hbar(\textbf{x},\mathbf{\varepsilon})$',...
                    'interpreter','latex','fontsize',fontsize_labels);
                % Then add title
                evalstr_t = strcat(method,':',{' '},...
                    ['Chain convergence plot of ABC-fitness function, ' ...
                    '$\hbar(\textbf{x},\mathbf{\varepsilon})$']);
                title(ax_gen,strcat(evalstr_t,sndwch_text),'fontsize', ...
                    fontsize_title,...
                    'interpreter','latex');
            elseif ismember(DREAMPar.lik,23)
                ylabel(ax_gen,'$\hbar(\textbf{x},\mathbf{\Delta})$',...
                    'interpreter','latex','fontsize',fontsize_labels);
                % Then add title
                evalstr_t = strcat(method,':',{' '},...
                    ['Chain convergence plot fitness of limits of ' ...
                    'acceptability, $\hbar(\textbf{x},\mathbf{\Delta})$']);
                title(ax_gen,strcat(evalstr_t,sndwch_text),'fontsize', ...
                    fontsize_title,...
                    'interpreter','latex');
            end
            % Now move ylabel
            xylim = [ get(ax_gen,'XLim') get(ax_gen,'YLim') ];
            % Now add legend - manually to match colors to labels
            dx = (xylim(2)-xylim(1)); dy = (xylim(4)-xylim(3));
            x_loc0 = xylim(1) + [ 0.10 0.12 ] * dx;
            y_loc = xylim(3) + [ 0.94 0.94 ] * dy;
            % Plot legend manually
            for z = 1:M
                x_loc = x_loc0 + (z-1)/12 * dx;
                line(ax_gen,x_loc,y_loc,'color',color_order(z,1:3),...
                    'linewidth',linewidth_legend);
                text(ax_gen,x_loc(2)+0.005*dx,y_loc(1), ...
                    char(str_chain(z)),...
                    'interpreter','latex','fontsize',fontsize_legend,...
                    'color',color_order(z,1:3),...
                    'horizontalalignment','left');
            end
            x_loc = x_loc0 + M/12 * dx;
            % Now add MAP
            hmult = 400*height_plot;
            plot(ax_gen,x_loc(2),y_loc(1),'kx','markersize',10,...
                'linewidth',2);
            text(ax_gen,x_loc(2) + 0.005*dx,y_loc(1)-dy/hmult,...
                '$\;{\rm MAP}$','interpreter','latex','fontsize',...
                fontsize_legend,'color','k','horizontalalignment',...
                'left','verticalalignment','middle');
        case 2
            plot(ax_gen, DREAMPar.T , MAP_lik + MAP_pr,'kx',...
                'Markersize',fontsize_remark,'linewidth',3);
            % Then add a y-label
            if ~ismember(DREAMPar.lik,[22 23])
                ylabel(ax_gen,...
                    ['$\mathcal{P}(\textbf{x}) + ' ...
                    '\mathcal{L}(\textbf{x}|\widetilde{\textbf{Y}})$'],...
                    'interpreter','latex','fontsize',fontsize_labels);
                % Then add title
                evalstr_t = strcat(method,':',{' '},...
                    ['Chain convergence of log-density, ' ...
                    '$\mathcal{P}(\textbf{x}) + ' ...
                    '\mathcal{L}(\textbf{x}|\widetilde{\textbf{Y}})$']);
                title(ax_gen,evalstr_t,'fontsize',fontsize_title, ...
                    'interpreter','latex');
            elseif ismember(DREAMPar.lik,22)
                ylabel(ax_gen,'$\hbar(\textbf{x},\varepsilon)$',...
                    'interpreter','latex','fontsize',fontsize_labels);
                % Then add title
                evalstr_t = strcat(method,':',{' '},...
                    ['Chain convergence plot of ABC-fitness function,' ...
                    ' $\hbar(\textbf{x},\varepsilon)$']);
                title(ax_gen,evalstr_t,'fontsize',fontsize_title,...
                    'interpreter','latex');
            elseif ismember(DREAMPar.lik,23)
                ylabel(ax_gen,'$\hbar(\textbf{x},\mathbf{\Delta})$',...
                    'interpreter','latex','fontsize',fontsize_labels);
                % Then add title
                evalstr_t = strcat(method,':',{' '},...
                    ['Chain convergence plot fitness of limits of ' ...
                    'acceptability, $\hbar(\textbf{x},\mathbf{\Delta})$']);
                title(ax_gen,evalstr_t,'fontsize',fontsize_title,...
                    'interpreter','latex');
            end
            % Add a title
            xlh = xlabel(ax_gen,'${\rm Sample\;number\;of\;chain}$',...
                'fontsize',fontsize_labels,'interpreter','latex');
            xylim = [ get(ax_gen,'XLim') get(ax_gen,'YLim') ];
            xlh.Position(2) = xylim(3) - 0.18 * (xylim(4)-xylim(3));
    end
    % move ylabel
    %ylh.Position(1) = xylim(1) - 0.05 * (xylim(2) - xylim(1));

    % Add label to each figure
    fig_code = strcat('(',ranktoletter(id_plot),')');
    text(ax_gen,0.005,0.94,fig_code,'units','normalized',...
        'fontsize',fontsize_A,'interpreter','latex',...
        'horizontalalignment','left'); box off;
    % Plot box around figure;
    plot_box(ax_gen); set(gcf,'color','w');
end

%% <><><><><>< HISTOGRAMS OF MARGINAL DENSITIES OF PARAMETERS ><><><><><><>
n_row = 2; n_col = 4;
% Now determine row and column number of each parameter
[row_par,col_par,idx_fig] = deal(nan(DREAMPar.d,1));
[row_par(1),col_par(1),idx_fig(1)] = deal(1);
for d = 2:DREAMPar.d
    col_par(d) = col_par(d-1)+1; row_par(d) = row_par(d-1);
    if col_par(d) == n_col + 1
        row_par(d) = row_par(d-1)+1; col_par(d) = 1;
    end
    if row_par(d) == n_row + 1
        row_par(d) = 1; idx_fig(d) = 1;
    else
        idx_fig(d) = 0;
    end
end

% Compute number of bins based on different rules
Nbins = nan(1,DREAMPar.d);
for i = 1:DREAMPar.d
    Nbins(i) = calcnbins(P_post(:,i));
end

% Take the minimum of the number of bins for each parameter
nbins = min(min(Nbins),maxbins);
% Added 2024 as there are too many bins
nbins = max(5,round(nbins/2));
% End added

% Create title
title_str = strcat(char(strcat(method,':',{' '},...
    'Marginal posterior distribution of parameters')),sndwch_text);
% Now loop
for par = 1:DREAMPar.d
    if idx_fig(par) == 1            % Open new figure/add title to previous one
        if exist('ax_hist','var')   % Add title to previous figure
            mtit(char(title_str),'fontsize',fontsize_title,...
                'interpreter','latex');
        end
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    % Define first axis
    ax_hist = axes('units','normalized'); axpos_hist = [ ...
        0.06 + (col_par(par)-1) * 0.24 , ...
        0.56 - (row_par(par)-1) * 0.48 , 0.20 , 0.38 ];
    set(ax_hist,'position',axpos_hist);
    % New implementation using histcounts
    [M,edges] = histcounts(P_post(:,par),nbins,'normalization',...
        'countdensity');
    X = 1/2*(edges(1:nbins)+edges(2:nbins+1)); % midpoint of each bin
    % And plot histogram in red
    h = bar(ax_hist,X,M./max(M),'r'); hold on; % --> can be scaled to 1 if using "trapz(X,N)" instead of "sum(N)"!
    set(h,'Facecolor',[0.5 0.5 0.5],'EdgeColor','w');
    % Now determine the ranges of X
    % %    [xy_values,X_min,X_max] = get_xticks(pars_post(:,par),4)
    % %     set(ax_hist,'fontsize',fontsize_axis_numbers,'tickdir','out',...
    % %         'ticklength',[0.03 0.05],'XTick',xy_values(par,1:nx),...
    % %         'XTickLabel',xy_values(par,1:nx),'XMinorTick','on',...
    % %         'YMinorTick','on','XLim',[min(pars_post(:,par)) max(pars_post(:,par))]);
    % %     set(ax_hist,'fontsize',fontsize_axis_numbers,'tickdir','out',...
    % %         'ticklength',[0.03 0.05],'XTick',xy_values(1,1:nx),...
    % %         'XTickLabel',xy_values(1,1:nx),'XMinorTick','on',...
    % %         'YMinorTick','on','XLim',[X_min X_max]);
    set(ax_hist,'fontsize',fontsize_axis_numbers,'tickdir','out',...
        'ticklength',[0.03 0.05],'XMinorTick','on',...
        'YMinorTick','on'); axis tight;
    % Define xtick values
    xlabel(ax_hist,str_par(par),'fontsize',fontsize_labels,...
        'interpreter','latex');

    % Add legend first subplot
    switch col_par(par)
        case 1  % Add y-label
            set(ax_hist,'ytick',0:0.2:1.0,'yticklabel',...
                {'0.0','0.2','0.4','0.6','0.8','1.0'});
            ylabel(ax_hist,'Empirical density','fontsize',fontsize_labels);
        otherwise
            set(ax_hist,'ytick',0:0.2:1.0,'yticklabel',[]);
    end
    % Now determine the min and max X values of the plot
    %    minX = min(X); maxX = max(X); maxY = 1.02;
    % Now determine appropriate scales
    %    deltaX = 0.1*(maxX - minX);
    % Calculate x_min and x_max
    %   x_min = minX - deltaX; x_max = maxX + deltaX;
    % Lets add the MAP value
    plot(ax_hist,MAP(par),0,'kx','Markersize',15,...
        'linewidth',linewidth_marker);
    fix_ticklabels(ax_hist,'x');
    % Adjust the axis
    %    axis(ax_hist,[x_min x_max 0 maxY]);
    % Add label
    fig_code = strcat('(',ranktoletter(par),')');
    text(ax_hist,0.02,0.94,fig_code,'units','normalized',...
        'fontsize',fontsize_A,'interpreter','latex',...
        'horizontalalignment','left');
    % If DREAMPar.d then add title
    if par == DREAMPar.d
        mtit(char(title_str),'fontsize',fontsize_title,...
            'interpreter','latex');
    end
    % Plot box around figure;
    plot_box(ax_hist); set(gcf,'color','w');
end

%% <><><><><><><><> PLOT MAP WITH CORRELATION ESTIMATES ><><><><><><><><><>
if ismember(DREAMPar.d,2:50)
    % Now plot the R_statistic for each parameter
    figure('units','normalized','outerposition',[0 0 s0(4)/s0(3)+0.05 1])
    % Plot
    ax_corr = axes('units','normalized'); axpos_corr = ...
        [ 0.06 0.06 0.9 0.9 ];
    set(ax_corr,'position',axpos_corr);
    imagesc(ax_corr,CORR); hcol = colorbar; set(hcol,'tickdir','out',...
        'fontsize',fontsize_axis_numbers);
    tickvalues = get(hcol,'Ticks'); tickvalues_fine = cell(1, ...
        numel(tickvalues));
    for i = 1:numel(tickvalues)
        tickvalues_fine{i} = num2str(tickvalues(i),'%3.2f');
    end
    set(hcol,'ticks',tickvalues,'ticklabels',tickvalues_fine,'tickdir', ...
        'out');
    % Now change axis
    set(ax_corr,'fontsize',fontsize_axis_numbers);
    % adjust position of ticks
    set(ax_corr,'XTick',1:DREAMPar.d+1,'xticklabel',[]);
    set(ax_corr,'YTick',1:DREAMPar.d+1,'yticklabel',[]);
    for j = 1 : DREAMPar.d
        % Now add labels as well
        h = text(ax_corr,j,DREAMPar.d + .5 + 0.055*DREAMPar.d,str_par(j),...
            'fontsize',fontsize_labels-2); set(h, 'rotation', 90)
        text(ax_corr,0.5 -0.055*DREAMPar.d,j,str_par(j),'fontsize',...
            fontsize_labels-2);
    end
    % set labels
    set(ax_corr,'xtick', linspace(0.5,DREAMPar.d-0.5,DREAMPar.d), ...
        'ytick',linspace(0.5,DREAMPar.d-0.5,DREAMPar.d));
    set(ax_corr,'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-', ...
        'xcolor', 'k', 'ycolor', 'k');
    % Add title
    title(ax_corr,strcat(method,':',{' '},...
        'Map of correlation coefficients of posterior parameter samples', ...
        sndwch_text),'fontsize',fontsize_title);
    % Set font
    set(ax_corr,'fontsize',fontsize_axis_numbers); set(gcf,'color','w');
    % Plot box around figure;
    %plot_box(ax_corr);
elseif ( DREAMPar.d == 1 )
    fprintf(['DREAM_Suite WARNING: CANNOT PLOT BIVARIATE SCATTER PLOTS' ...
        ' AS DREAMPAR.d = 1\n']);
else
    fprintf(['DREAM_Suite WARNING: CANNOT PLOT MAP WITH CORRELATION ' ...
        'ESTIMATES AS DREAMPAR.d = %1d\n'],DREAMPar.d);
end

%% <><><><>< CORRELATION PLOTS OF THE POSTERIOR PARAMETER SAMPLES ><><><><>
prior_domain = [1 0 0]; % red; used to be blue: [0 0.4470 0.7410];
reg_line = [51 102 255]/255;
% Only plot this matrix if less or equal to 30 parameters
if ismember(DREAMPar.d,2:30)

    % Determine the number of bins for main diagonal (histograms)
    Nbins = nan(1,DREAMPar.d);
    for i = 1:DREAMPar.d
        Nbins(i) = calcnbins(P_post(:,i));
    end
    % Take the minimum of the number of bins for each parameter
    nbins = min(min(Nbins),20);
    % maximum of 20 bins is maximum otherwise may not show in matrixplot
    % Added 2024 as there are too many bins
    nbins = max(5,round(nbins/2));
    % End added

    % How many x-values should be printed on each axis?
    nx = 3;
    [xy_values,x_min,x_max] = det_print_values(P_post(:,1:DREAMPar.d),nx);

    % Create figure
    figure('units','normalized','outerposition',[0 0 1 1]);
    tick_length = [0.02 * DREAMPar.d/5 0.02];
    ax_corr = axes('units','normalized'); axpos_corr = [ 0.05 0.08 0.9 0.85 ];
    set(ax_corr,'position',axpos_corr);
    % Plot a matrix (includes unscaled marginals on main diagonal!
    [H,AX,~,Pl,PAx] = plotmatrix(ax_corr,P_post(:,1:DREAMPar.d),'rs'); 
    hold on;
    % Add title
    title(ax_corr,strcat(method,':',{' '},...
        ['Marginal distribution and bivariate scatter plots of posterior ' ...
        'samples'],sndwch_text),'fontsize',fontsize_title);
    % Now add the prior ranges - if hypercube
    if isfield(Par_info,'boundhandling')
        for i = 1: DREAMPar.d
            for j = 1 : DREAMPar.d
                if i ~= j
                    hold(AX(j,i),'on');
                    % Vertical lines
                    line(AX(j,i),[ Par_info.min(i) Par_info.min(i) ],...
                        [ Par_info.min(j) Par_info.max(j)],...
                        'color',prior_domain);
                    line(AX(j,i),[ Par_info.max(i) Par_info.max(i) ],...
                        [ Par_info.min(j) Par_info.max(j)],...
                        'color',prior_domain);
                    % Horizontal lines
                    line(AX(j,i),[ Par_info.min(i) Par_info.max(i) ],...
                        [ Par_info.min(j) Par_info.min(j)],...
                        'color',prior_domain);
                    line(AX(j,i),[ Par_info.min(i) Par_info.max(i) ],...
                        [ Par_info.max(j) Par_info.max(j)],...
                        'color',prior_domain);
                end
            end
        end
    end
    % Now change all axis to lie between x_min and x_max, and format
    for i = 1 : DREAMPar.d
        for j = 1 : DREAMPar.d
            hold(AX(j,i),'on');
            % Define axis limits
            set(AX(j,i),'XLim',[x_min(i) x_max(i)],'YLim',[x_min(j) ...
                x_max(j)]);
            % Add least squares fit off-diagonal scatter plots
            switch i == j
                case 1      % Histogram: Define limits and bin color/settings
                    hold(PAx(i),'on'); set(PAx(i),'Xlim',[x_min(i) ...
                        x_max(i)]);
                    if verLessThan('matlab','9.1')
                        set(Pl(i),'Facecolor',light_gray,'EdgeColor','w');
                    else
                        set(Pl(i),'NumBins',nbins,'Facecolor',light_gray,...
                            'EdgeColor','w','BinLimits',[x_min(i) x_max(i)]);
                    end
                    plot(PAx(i),MAP(i),0,'kx','linewidth',...
                        thickness_line_in_plot,'markersize',10);
                otherwise   % Scatter plot: Define marker + add least squares line
                    set(H(j,i),'Marker','s','MarkerSize',4,'color', ...
                        light_gray);
                    c = polyfit(P_post(:,i),P_post(:,j),1);
                    x_lim = get(AX(j,i),'XLim');
                    line(AX(j,i),x_lim,c(1)*x_lim + c(2),'color',reg_line,...
                        'linestyle','--','linewidth',1);
                    plot(AX(j,i),MAP(i),MAP(j),'kx','linewidth',...
                        thickness_line_in_plot,'markersize',10);
            end
        end
    end
    % Remove entries of H plot in (1,1) and (DREAMPar.d,DREAMPar.d) as they show underneath
    set(H(1,1),'Marker','none'); set(H(DREAMPar.d,DREAMPar.d), ...
        'Marker','none');
    % Set the y labels and tick direction on all plots left column
    for i = 1:size(AX,1)-1
        set(AX(i,1),'tickdir','out','box','off','xtick',[],...
            'xticklabels',[],'ticklength',tick_length,'ytick',...
            xy_values(i,1:nx),'yticklabels',xy_values(i,1:nx));
        fix_ticklabels(AX(i,1),'y');
        % Now add outer box
        xylim = [ get(AX(i,1),'XLim') get(AX(i,1),'YLim') ];
        line(AX(i,1),[xylim(2) xylim(2)],[xylim(3) xylim(4)],'color','k',...
            'linewidth',0.25);
        line(AX(i,1),[xylim(1) xylim(2)],[xylim(4) xylim(4)],'color','k',...
            'linewidth',0.25);
    end
    % Add ylabels
    ymult = 0.23*DREAMPar.d/10;
    for i = 1:DREAMPar.d
        xylim = [ get(AX(i,1),'XLim') get(AX(i,1),'YLim') ];
        ylh = ylabel(AX(i,1),str_par(i),'fontsize',fontsize_labels);
        ylh.Position(1) = xylim(1) - ymult * (xylim(2)-xylim(1));
    end
    set(AX(DREAMPar.d,1),'tickdir','out','box','off',...
        'xtick',xy_values(1,1:nx),'xticklabels',xy_values(1,1:nx),...
        'ytick',xy_values(DREAMPar.d,1:nx),'yticklabels',...
        xy_values(DREAMPar.d,1:nx),'ticklength',tick_length);
    fix_ticklabels(AX(DREAMPar.d,1),'x'); 
    fix_ticklabels(AX(DREAMPar.d,1),'y');
    % Now add outer box
    xylim = [ get(AX(DREAMPar.d,1),'XLim') get(AX(DREAMPar.d,1),'YLim') ];
    line(AX(DREAMPar.d,1),[xylim(2) xylim(2)],[xylim(3) xylim(4)],...
        'color','k','linewidth',0.25);
    line(AX(DREAMPar.d,1),[xylim(1) xylim(2)],[xylim(4) xylim(4)],...
        'color','k','linewidth',0.25);
    for i = 2:size(AX,2)
        set(AX(DREAMPar.d,i),'tickdir','out','box','off',...
            'xtick',xy_values(i,1:nx),'xticklabels',xy_values(i,1:nx),...
            'ytick',[],'yticklabels',[],'ticklength',tick_length);
        % Fix_ticklabels
        fix_ticklabels(AX(DREAMPar.d,i),'x');
        xylim = [ get(AX(DREAMPar.d,i),'XLim') get(AX(DREAMPar.d,i), ...
            'YLim') ];
        line(AX(DREAMPar.d,i),[xylim(2) xylim(2)],[xylim(3) xylim(4)],...
            'color','k','linewidth',0.25);
        line(AX(DREAMPar.d,i),[xylim(1) xylim(2)],[xylim(4) xylim(4)],...
            'color','k','linewidth',0.25);
    end
    % Add xlabels
    xmult = 0.3*DREAMPar.d/10;
    for i = 1:DREAMPar.d
        xylim = [ get(AX(DREAMPar.d,i),'XLim') get(AX(DREAMPar.d,i), ...
            'YLim') ];
        xlh = xlabel(AX(DREAMPar.d,i),str_par(i),'fontsize',fontsize_labels);
        xlh.Position(2) = xylim(3) - xmult * (xylim(4)-xylim(3));
    end

    % Delete scatter plots on main diagonal (2:DREAMPar.d) and set ticks
    % empty on all other figures
    for i = 1:DREAMPar.d-1
        for j = 2:DREAMPar.d
            if i == j % Delete scatter plot on diagonal
                delete(AX(i,j)); delete(H(i,j));
                set(PAx(i),'xtick',[],'ytick',[]);
            else
                set(AX(i,j),'xtick',[],'ytick',[]);
            end
        end
    end
    % Add legend to plot in top left, i = 1
    i = 1;
    xylim = [ get(PAx(i),'XLim') get(PAx(i),'YLim') ];
    dx = (xylim(2)-xylim(1)); dy = (xylim(4)-xylim(3));
    x_loc = xylim(1) + [ 0.04 0.14 ] * dx;
    y_loc = xylim(3) + [ 0.85 0.85 ] * dy;
    hmult = 400*PAx(1).Position(4);
    % Plot legend manually
    plot(PAx(i),sum(x_loc)/2,y_loc(1),'kx','markersize',...
        max(6,fontsize_remark - DREAMPar.d/2),'linewidth', ...
        thickness_line_in_plot);
    text(PAx(i),x_loc(2) + 0.02*dx,y_loc(1)-dy/hmult,'${\rm MAP}$',...
        'interpreter','latex','fontsize',max(8,fontsize_legend - ...
        DREAMPar.d/2),'color','k',...
        'horizontalalignment','left','verticalalignment','middle');
    set(gcf,'color','w');    
    % Horizontal line (legend for regression line)
    % %     xylim = [ get(AX(1,2),'XLim') get(AX(1,2),'YLim') ];
    % %     dx = (xylim(2)-xylim(1)); dy = (xylim(4)-xylim(3));
    % %     x_loc = xylim(1) + [ 0.04 0.14 ] * dx;
    % %     y_loc = xylim(3) + [ 0.85 0.85 ] * dy;
    % %     line(AX(1,2),x_loc,y_loc,'color',reg_line,'linestyle',':','linewidth',...
    % %             max(1,linewidth_legend - DREAMPar.d/10));
    % %     text(AX(1,2),x_loc(2)+0.02*dx,y_loc(1)-dy/hmult,'$\;{\rm Regression\;line}$',...
    % %             'interpreter','latex','fontsize',max(8,fontsize_legend - DREAMPar.d/3),'color',reg_line,...
    % %             'horizontalalignment','left');
elseif ( DREAMPar.d == 1 )
    fprintf(['DREAM_Suite WARNING: CANNOT PLOT BIVARIATE SCATTER PLOTS ' ...
        'AS DREAMPAR.d = 1\n']);
else
    fprintf(['DREAM_Suite WARNING: CANNOT PLOT BIVARIATE SCATTER PLOTS ' ...
        'AS DREAMPAR.d = %1d (TOO LARGE)\n'],DREAMPar.d);
end

%% <><><><><><>< PLOT HISTOGRAMS OF THE SUMMARY STATISTICS <><><><><><><><>
color_limits = [1 0.5 0.5];
if ~isempty(S_post) && ~ismember(DREAMPar.lik,23) 
    % --> summary statics used but not limits of acceptability
    
    % Plot the sample autocorrelation functions in a layout of nrow,ncol
    n_row = 2; n_col = 4;
    % Now determine row and column number of each parameter
    [row_par,col_par,idx_fig] = deal(nan(DREAMPar.d,1));
    [row_par(1),col_par(1),idx_fig(1)] = deal(1);
    for d = 2:DREAMPar.d
        col_par(d) = col_par(d-1)+1; row_par(d) = row_par(d-1);
        if col_par(d) == n_col + 1
            row_par(d) = row_par(d-1)+1; col_par(d) = 1;
        end
        if row_par(d) == n_row + 1
            row_par(d) = 1; idx_fig(d) = 1;
        else
            idx_fig(d) = 0;
        end
    end
    % How many bins should histogram use?
    Nbins = nan(1,Meas_info.n_S);
    for sm = 1:Meas_info.n_S
        Nbins(sm) = calcnbins(S_post(:,sm));
    end
    % Take the minimum of the number of bins for each parameter
    nbins = min(min(Nbins),maxbins);
    % Added 2024 as there are too many bins
    nbins = max(5,round(nbins/2));
    % End added

    % Maximum value of pdf
    max_pdf = 1.02;
    % Create title
    title_str = strcat(char(strcat(method,':',{' '},...
        'Marginal distribution of summary statistics')),sndwch_text);
    clear ax_Shist
    for sm = 1:Meas_info.n_S
        if idx_fig(sm) == 1            % Open new figure/add title to previous one
            if exist('ax_Shist','var') % Add title to previous figure
                mtit(char(title_str),'fontsize',fontsize_title,...
                    'interpreter','latex');
            end
            figure('units','normalized','outerposition',[0 0 1 1]);
        end
        % Define first axis
        ax_Shist = axes('units','normalized'); axpos_Shist = [ 0.06 ...
            + (col_par(sm)-1) * 0.24 , ...
            0.56 - (row_par(sm)-1) * 0.48 , 0.20 , 0.38 ];
        set(ax_Shist,'position',axpos_Shist);
        % Plot the histogram of each parameter
        [M,edges] = histcounts(S_post(:,sm),nbins,'normalization', ...
            'countdensity');
        Xbin = 1/2*(edges(1:end-1)+edges(2:end)); % midpoint of each bin
        % And plot histogram in red
        h = bar(ax_Shist,Xbin,M./max(M),'r'); hold on; % --> can be
        % scaled to 1 if using "trapz(X,N)" instead of "sum(N)"!
        set(h,'Facecolor',medium_gray,'EdgeColor','w');
        set(ax_Shist,'ytick',0:0.2:1.0,'yticklabel',[]);
        set(ax_Shist,'fontsize',fontsize_axis_numbers,'tickdir','out',...
            'ticklength',[0.03 0.05]);
        % Add x-labels
        s_str = strcat('$S_{',num2str(sm),'}$');
        xlabel(ax_Shist,s_str,'fontsize',fontsize_labels, ...
            'interpreter','latex');
        % Add legend first subplot
        switch col_par(sm)
            case 1  % Add y-label
                ylabel(ax_Shist,'${\rm Empirical\;density}$',...
                    'fontsize',fontsize_labels,'interpreter','latex');
                set(ax_Shist,'ytick',0:0.2:1.0,'yticklabel',...
                    {'0.0','0.2','0.4','0.6','0.8','1.0'});
            otherwise % do nothing
        end
        % Lets add the observed summary statistic
        S_obs = Meas_info.S(sm);
        plot(ax_Shist,S_obs,0,'rx','Markersize',12,'linewidth',2);
        % Now put vertical dotted bars for epsilon
        Sax_left = [S_obs - options.epsilon(sm) S_obs - options.epsilon(sm)];
        % Separate in two lines
        plot(ax_Shist,Sax_left,[0 0.4],'linestyle','--','color',...
            color_limits,'linewidth',thickness_line_in_plot);
        plot(ax_Shist,Sax_left,[0.8 max_pdf],'linestyle','--','color',...
            color_limits,'linewidth',thickness_line_in_plot);
        Sax_right = [S_obs + options.epsilon(sm) S_obs + options.epsilon(sm)];
        plot(ax_Shist,Sax_right,[0 0.4],'linestyle','--','color',...
            color_limits,'linewidth',thickness_line_in_plot);
        plot(ax_Shist,Sax_right,[0.8 max_pdf],'linestyle','--','color',...
            color_limits,'linewidth',thickness_line_in_plot);
        % Adjust the axis
        axis(ax_Shist,[S_obs - 3*options.epsilon(sm) S_obs + ...
            3*options.epsilon(sm) 0 max_pdf ]);
        set(ax_Shist,'XMinorTick','on','YMinorTick','on');
        % Add text on left and right vertical line
        str_left = strcat('$\widetilde{S}_{',num2str(sm),...
            '} - \epsilon_{',num2str(sm),'}$');
        str_right = strcat('$\widetilde{S}_{',num2str(sm),...
            '} + \epsilon_{',num2str(sm),'}$');
        % Get the exact axis to plot text in figure
        dx = diff(xlim); dy = diff(ylim);
        text(ax_Shist,Sax_left(1),0.6,char(str_left),...
            'rotation',90,'interpreter','latex','color',color_limits,...
            'horizontalalignment','center','fontsize',fontsize_legend);
        text(ax_Shist,Sax_right(1),0.6,char(str_right),...
            'rotation',90,'interpreter','latex','color',color_limits,...
            'horizontalalignment','center','fontsize',fontsize_legend);
        % And the significant digits on x-axis numbers
        fix_ticklabels(ax_Shist,'x');
        % Now add legend - manually to match colors to labels
        xylim = [ get(ax_Shist,'XLim') get(ax_Shist,'YLim') ];
        x_loc = xylim(1) + [ 0.82 0.88 ] * dx;
        y_loc = xylim(3) + [ 0.94 0.94 ] * dy;
        % Plot legend manually
        plot(ax_Shist,sum(x_loc)/2,y_loc(1),'rx','markersize',13, ...
            'linewidth',2);
        str_Shist = strcat('$\widetilde{S}_{',num2str(sm),'}$');
        text(ax_Shist,x_loc(2) + 0.02*dx,y_loc(1),char(str_Shist),...
            'interpreter','latex','fontsize',fontsize_legend,'color','r',...
            'horizontalalignment','left','verticalalignment','middle');
        % Add label to each figure
        fig_code = strcat('(',ranktoletter(sm),')');
        text(ax_Shist,0.02,0.94,fig_code,'units','normalized',...
            'fontsize',fontsize_A,'interpreter','latex',...
            'horizontalalignment','left');
        % If sm equal to Meas_info.n_S then add title
        if sm == Meas_info.n_S
            mtit(char(title_str),'fontsize',fontsize_title,...
                'interpreter','latex');
        end
        % Plot box
        plot_box(ax_Shist);
    end
    set(gcf,'color','w');
end

%% <><><><><><>< PLOT TIME SERIES OF THE SUMMARY STATISTICS ><><><><><><><>
if ~isempty(S_post) && ismember(DREAMPar.lik,23) % --> summary statistics:
    % Limits of acceptability
    % Open new figure
    figure('units','normalized','outerposition',[0 0 1 1])
    % Find solutions that satisfy limits of acceptability
    idx_LOA = find ( P_post(:,DREAMPar.d+2) == Meas_info.n_S );
    % Now add the observations
    ax10 = axes('units','normalized'); ax_genpos = [ 0.1 0.65 0.8 0.3 ];
    set(ax10,'position',ax_genpos);
    ax10.ColorOrder = color_order;
    % If isempty(idx) --> no posterior solutions
    if isempty(idx_LOA)
        plot(ax10,1:Meas_info.n_S,Meas_info.S,'r.', ...
            'markersize',fontsize_remark);
        axis tight; a = axis; axis([ a(1:2) min(0.9*a(3),1.1*a(3)) ...
            max(0.9*a(4),1.1*a(4)) ]);
        text(ax10,0.05,0.90,['No solution satisfies limits of ' ...
            'acceptability of all observations'],...
            'fontsize',fontsize_labels,'units','normalized');
    else
        % Derive minimum and maximum simulated values
        tot_S_unc(:,1) = min(S_post(idx_LOA,1:Meas_info.n_S));
        tot_S_unc(:,2) = max(S_post(idx_LOA,1:Meas_info.n_S));
        % Plot total simulation uncertainty
        Fill_Ranges(1:Meas_info.n_S,tot_S_unc(:,1),tot_S_unc(:,2),...
            light_gray,ax10); hold on;
        % Now add the observations
        plot(ax10,1:Meas_info.n_S,Meas_info.S,'r.', ...
            'markersize',fontsize_remark);
        axis tight; a = axis; axis(ax10,[ a(1:2) min(0.9*a(3),1.1*a(3)) ...
            max(0.9*a(4),1.1*a(4)) ]);
        % Add labels
        set(ax10,'fontsize',fontsize_axis_numbers,'tickdir','out');
        xlabel(ax10,'Observation number','fontsize',fontsize_labels);
        ylabel(ax10,'Data/simulated values','fontsize',fontsize_labels);

        a = axis;
        xloc = a(1) + 0.8*(a(2)-a(1)); deltax = 0.02*(a(2)-a(1));
        yloc = a(3) + 0.90*(a(4)-a(3)); deltay = 0.04*(a(4)-a(3));
        patch(ax10,[xloc xloc xloc+deltax xloc+deltax],[yloc yloc+deltay ...
            yloc+deltay yloc],...
            light_gray,'EdgeColor',light_gray);
        text(ax10,xloc+1.5*deltax,yloc + 4/10*deltay, ...
            'Posterior uncertainty',...
            'fontsize',fontsize_legend,'color',light_gray);
        yloc = a(3) + 0.8*(a(4)-a(3)); %deltay = 0.03*(a(4)-a(3));
        plot(ax10,xloc+deltax/2,yloc + deltay/2,'ro','markerfacecolor', ...
            'r','markersize',8);
        text(ax10,xloc+1.5*deltax,yloc + 4/10*deltay,'Measurement data',...
            'fontsize',fontsize_legend,'color','r');

        xt = xticklabels(ax10); xaxis_numbers = nan(1,numel(xt));
        for i = 1:numel(xt)
            xaxis_numbers(i) = str2double(char(xt(i,:)));
        end
        set(ax10,'xtick',xaxis_numbers,'xticklabels', ...
            addcomma(xaxis_numbers));

        % Now calculate percentage inside the PredInt bound defined in Bayes_pdf
        contained_LOA = 100 * (1 - length( find ( Meas_info.S < ...
            tot_S_unc(:,1) | Meas_info.S > tot_S_unc(:,2) ) ) / ...
            Meas_info.n_S );
        % This should be close to the alfa interval that was used
        evalstr = strcat('Limits of acceptability envelope',{' '},...
            num2str(contained_LOA,'%3.2f'),'\%',{' '},'of observations');
        text(ax10,0.05,0.90,evalstr,'interpreter','latex','fontsize',...
            fontsize_labels,'units','normalized','horizontalalignment', ...
            'left');
    end

    title(ax10,strcat(method,':',{' '},...
        ['Limits of acceptability: Ranges of posterior simulations that ' ...
        'satisfy limits of each observation'],sndwch_text),...
        'interpreter','latex','fontsize',fontsize_title);
    % Add label to each figure
    fig_code = strcat('(',ranktoletter(1),')');
    text(ax10,0.01,0.90,fig_code,'units','normalized',...
        'fontsize',fontsize_A,'interpreter','latex',...
        'horizontalalignment','left'); box off;
    % Plot box
    plot_box(ax10);

    % Now second axes
    ax11 = axes('units','normalized'); ax_genpos = [ 0.1 0.15 0.8 0.30 ];
    set(ax11,'position',ax_genpos); height_plot = 0.3;
    ax11.ColorOrder = color_order;
    % Define M
    M = min(DREAMPar.N,5);
    % Define x-axis values and number
    switch DREAMPar.thinning
        case 1
            xaxis_values = 1:DREAMPar.T;
        otherwise
            xaxis_values = [1 DREAMPar.thinning:DREAMPar.thinning:DREAMPar.T];
    end
    xaxis_numbers = [1 DREAMPar.T/5 2*DREAMPar.T/5:DREAMPar.T/5:DREAMPar.T];

    % Now plot a number of chains
    for i = 1:M
        % char(symbol(i)) replaced with color_order(i,1:3);
        plot(ax11,xaxis_values,chain(1:end,DREAMPar.d+2,i),'color', ...
            color_order(i,1:3),...
            'linewidth',thickness_line_in_plot);
        if i == 1; hold on; end
    end
    % Set axis
    set(ax11,'fontsize',fontsize_axis_numbers,'tickdir','out');
    % Ranges have not been defined -- need to derive them from ParSet
    axis tight; a = axis; axis(ax11,[1 DREAMPar.T 0 max(0.9*a(4),1.1*a(4)) ])
    set(ax11,'xtick',xaxis_numbers,'xticklabels',addcomma(xaxis_numbers));
    % Lets add the MAP value
    plot(ax11, DREAMPar.T , max(P_post(:,DREAMPar.d+2)),'kx','Markersize',...
        fontsize_remark,'linewidth',linewidth_marker);
    % Add a title
    xlabel(ax11,'Sample number of chain','fontsize',fontsize_labels);
    % Then add a y-label
    ylabel(ax11,'Number of limits satisfied','fontsize',fontsize_labels);
    % Then add title
    title(ax11,strcat(method,':',{' '},...
        ['Chain convergence plot of no. of limits of acceptability ' ...
        'satisfied'],sndwch_text),...
        'interpreter','latex','fontsize',fontsize_title);

    % Now add legend - manually to match colors to labels
    xylim = [ get(ax11,'XLim') get(ax11,'YLim') ];
    dx = (xylim(2)-xylim(1)); dy = (xylim(4)-xylim(3));
    x_loc0 = xylim(1) + [ 0.10 0.12 ] * dx;
    y_loc = xylim(3) + [ 0.08 0.08 ] * dy;
    % Plot legend manually
    for z = 1:M
        x_loc = x_loc0 + (z-1)/10 * dx;
        line(ax11,x_loc,y_loc,'color',color_order(z,1:3),'linewidth',...
            linewidth_legend);
        text(ax11,x_loc(2)+0.005*dx,y_loc(1),char(str_chain(z)),...
            'interpreter','latex','fontsize',fontsize_legend,'color', ...
            color_order(z,1:3),...
            'horizontalalignment','left');
    end
    x_loc = x_loc0 + M/10 * dx;
    % Now add MAP
    hmult = 400*height_plot;
    plot(ax11,x_loc(2),y_loc(1),'kx','markersize',10,'linewidth',2);
    text(ax11,x_loc(2) + 0.005*dx,y_loc(1)-dy/hmult,'$\;{\rm MAP}$',...
        'interpreter','latex','fontsize',fontsize_legend,'color','k',...
        'horizontalalignment','left','verticalalignment','middle');
    % Add label to each figure
    fig_code = strcat('(',ranktoletter(2),')');
    text(ax11,0.01,0.93,fig_code,'units','normalized',...
        'fontsize',fontsize_A,'interpreter','latex',...
        'horizontalalignment','left'); box off;
    % Plot box
    plot_box(ax11); set(gcf,'color','w');
end

%% <><><><><>< PLOT THE 95% POSTERIOR SIMULATION UNCERTAINTY <><><><><><><>
if ~isempty(FX_post)
    
    % Compute confidence intervals (= due to parameter uncertainty) 
    alfa1 = 100*p_alfa/2;                       % Lower percentile
    alfa2 = 100*(1-p_alfa/2);                   % Upper percentile
    rng(1+round(100*rand),'twister');           % Random seed 
                                                % legacy: randn('state',sum(100*clock));
    par_unc = prctile(FX_post,[alfa1 alfa2])';  % p_alfa/2 & (1-p_alfa/2) confidence limits
    if Meas_info.n == 1, par_unc = par_unc'; end

    %% How to generate total uncertainty?
    % [1] Take MAP parameter values and add structural uncertainty 
    %     --> parameter uncertainty explains part of total uncertainty
    % [2] Sample posterior parameters and add structural uncertainty on top
    %     --> does not guarantee that 95% contains 95%
    unc_method = 1;
    % Method = 2: Does not work yet if modout = 'yes' -> need to add 
    % How: get "unique" samples of P_post
    switch unc_method
        case 1 % Using only MAP value
            % UP_postun = MAP; U_fx = FX_MAP'; Nr = size(FX_post,1); idU = ones(Nr,1); 
            MAP_un = X_unnormalize(MAP,Par_info);   % Transform to unnormalized space
            Nr = size(FX_post,1);                   % # Replicates desired
            [tot_unc,fx_mod] = Bayes_pdf(MAP_un,FX_MAP',ones(Nr,1),Nr,...
                RMSE_MAP,DREAMPar,Meas_info,Lik_info,p_alfa);            
        case 2 % Using posterior parameter samples
            nUP = numel(iiUP);                      % # unique samples   
            Nr = nan(nUP,1);                        % Vector with # replicates 
            for z=1:nUP, Nr(z) = sum(idUP==z); end  % How many replicates each row?
            % NEW: Get parameter and total uncertainty
            [tot_unc,fx_mod] = Bayes_pdf(UP_postun,...
                FX_post(iiUP,1:Meas_info.n),idUP,Nr,...
                RMSE_MAP,DREAMPar,Meas_info,Lik_info,p_alfa);
            % Ufx = FX_post(iiUP,1:Meas_info.n);    % FX of unique posterior samples
                                                    % (= sim_out of UP_postun: CHECKED)
    end
    % Compute spread of 95% confidence intervals
    width_par_unc = mean(par_unc(:,2)-par_unc(:,1));
    % Compute spread of 95% prediction intervals
    width_tot_unc = mean(tot_unc(:,2)-tot_unc(:,1));
    % Open new figure
    figure('units','normalized','outerposition',[0.1 0.1 0.8 0.7])
    % Now add the observations
    ax_unc = axes('units','normalized'); ax_uncpos = [ 0.07 0.13 0.9 0.8 ];
    set(ax_unc,'position',ax_uncpos);
    ax_unc.ColorOrder = color_order;
    % We start with the total uncertainty
    Fill_Ranges(1:size(tot_unc,1),tot_unc(:,1),tot_unc(:,2), ...
        [0.75 0.75 0.75],ax_unc); hold on;
    % And then plot the parameter uncertainty
    Fill_Ranges(1:size(tot_unc,1),par_unc(:,1),par_unc(:,2), ...
        [0.25 0.25 0.25],ax_unc);
    % Now add the observations
    plot(ax_unc,1:size(tot_unc,1),Meas_info.Y,'ro', ...
        'markerfacecolor','r','markersize',3);
    set(ax_unc,'fontsize',fontsize_axis_numbers,'tickdir','out');
    % Fit axes
    axis tight
    % Add labels
    xlabel(ax_unc,'Observation number','fontsize',fontsize_labels);
    ylabel(ax_unc,'Data/simulated values','fontsize',fontsize_labels);
    title(ax_unc,strcat(method,':',{' '},'Predictive uncertainty:',...
        {' '},num2str(alfa),'\% Posterior uncertainty',sndwch_text),...
        'interpreter','latex','fontsize',fontsize_title);
    a = axis;
    xloc = a(1) + 0.8*(a(2)-a(1)); deltax = 0.02*(a(2)-a(1));
    yloc = a(3) + 0.93*(a(4)-a(3)); deltay = 0.03*(a(4)-a(3));
    patch(ax_unc,[xloc xloc xloc+deltax xloc+deltax],[yloc yloc+deltay ...
        yloc+deltay yloc],...
        [0.75 0.75 0.75],'EdgeColor',[0.75 0.75 0.75]);
    text(ax_unc,xloc+1.5*deltax,yloc + 4/10*deltay,'Total uncertainty',...
        'fontsize',fontsize_legend,'color',[0.75 0.75 0.75]);
    yloc = a(3) + 0.87*(a(4)-a(3)); deltay = 0.03*(a(4)-a(3));
    patch(ax_unc,[xloc xloc xloc+deltax xloc+deltax],[yloc yloc+deltay ...
        yloc+deltay yloc],...
        [0.25 0.25 0.25],'EdgeColor',[0.25 0.25 0.25]);
    text(ax_unc,xloc+1.5*deltax,yloc + 4/10*deltay,'Parameter uncertainty',...
        'fontsize',fontsize_legend,'color',[0.25 0.25 0.25]);
    yloc = a(3) + 0.81*(a(4)-a(3));
    plot(ax_unc,xloc+deltax/2,yloc + deltay/2,'ro','markerfacecolor', ...
        'r','markersize',8);
    text(ax_unc,xloc+1.5*deltax,yloc + 4/10*deltay,'Measurement data',...
        'fontsize',fontsize_legend,'color','r');

    % Now calculate percentage inside the PredInt bound defined in Bayes_pdf
%     ctnd_par = 100 * sum( (Meas_info.Y >= par_unc(:,1)) & ...
%         (Meas_info.Y <= par_unc(:,2)) ) / Meas_info.n;
    ctnd_tot = 100 * sum( (Meas_info.Y >= tot_unc(:,1)) & ...
        (Meas_info.Y <= tot_unc(:,2)) ) / Meas_info.n;
    % This should be close to the alfa interval that was used
    evalstr = strcat('(i)',{' '},num2str(alfa),['\% uncertainty interval' ...
        ' envelops'],...
        {' '},num2str(ctnd_tot,'%4.2f'),'\%',{' '},'of observations');
    text(ax_unc,0.02,0.93,evalstr,'fontsize',fontsize_labels,...
        'units','normalized','interpreter','latex');
    % par_unc and tot_unc
    evalstr = strcat('(ii)',{' '},'width of',{' '},num2str(alfa),...
        '\% parameter uncertainty interval:',{' '}, ...
        num2str(width_par_unc,'%4.2f'));
    text(ax_unc,0.02,0.86,evalstr,'fontsize',fontsize_labels, ...
        'units','normalized');
    evalstr = strcat('(iii)',{' '},'width of',{' '},num2str(alfa),...
        '\% total uncertainty interval:',{' '},num2str(width_tot_unc, ...
        '%4.2f'));
    text(ax_unc,0.02,0.79,evalstr,...
        'fontsize',fontsize_labels,'units','normalized');
    ax_unc.XMinorTick = 'on'; ax_unc.YMinorTick = 'on';
    % Plot box
    plot_box(ax_unc); set(gcf,'color','w');

end

%% <><><><><><><><><><><><> NOW DO RESIDUAL ANALYSIS <><><><><><><><><><><>
if exist('FX_MAP','var') && (DREAMPar.lik > 10) && Meas_info.n > 1
    % Unnormalize MAP parameter values 
    MAP_un = X_unnormalize(MAP,Par_info);
    % Compute normalized (partial) residuals and corresponding density
    [~,e,~,X_n,XX_n] = Calc_likelihood(MAP_un,FX_MAP,DREAMPar,Par_info,...
        Meas_info,Lik_info,options,MAP_info);
    switch DREAMPar.lik
        case {13,14,16,17,44,45} % Treat residual correlation (possibly)
            eps_n = X_n; f_eps_n = XX_n;
            par = Lik_info.fpar;                 % all parameters
            par(Lik_info.id_vpar) = MAP_un;      % variable parameters
            nuisvar = par(Lik_info.id_nuisvar);  % isolate nuisance vars
        otherwise
            e_n = X_n; f_e_n = XX_n;    
    end

    % Do we have partial residuals or not - define text/legend on some axes
    switch DREAMPar.lik
        case {13,14,16,17,44,45} % Check if nuisance variables of AR are activated
            switch Lik_info.t_start
                case 1 % No treatment of serial correlation
                    idx_ll = 0;
                otherwise % (2,3,4,5); (2,3) for GL+/SGT/SST and (4,5) for GL
                    idx_ll = 1;
            end
        otherwise % No treatment of serial correlation
            idx_ll = 0;
    end
    switch idx_ll
        case 0
            partial_or_not0 = []; partial_or_not1 = [];
        case 1
            partial_or_not0 = 'part.\;'; partial_or_not1 = 'partial\;';
    end
    text_axis = strcat('${\rm St.\;',partial_or_not0,'residuals}$');
    leg_str2ax = strcat('${\rm Standardized\;',partial_or_not1, ...
        'residuals}$');
    t_end = Meas_info.n;

    %% A: Plot the time series of observed and simulated data
    figure('units','normalized','outerposition',[0 0 1 1]),
    % Define first axis
    ax_res = axes('units','normalized'); ax_respos = [ 0.05 0.72 0.9 0.27 ];
    set(ax_res,'position',ax_respos); height_plot = 0.27;
    plot(ax_res,1:t_end,FX_MAP,'color',color_order(1,1:3),'linewidth',1.5);
    hold on;
    plot(ax_res,1:t_end,Meas_info.Y,'ro','color',color_order(3,1:3),...
        'markerfacecolor',color_order(3,1:3),'markersize',3);
    ax_res.ColorOrder = color_order;
    axis(ax_res,[1 t_end min(min(FX_MAP),min(Meas_info.Y)) ...
        1.05*max(max(FX_MAP),max(Meas_info.Y))]);
    % Now add legend - manually to match colors to labels
    xylim = [ get(ax_res,'XLim') get(ax_res,'YLim') ];
    dx = (xylim(2)-xylim(1)); dy = (xylim(4)-xylim(3));
    x_loc = xylim(1) + [ 0.82 0.84 ] * dx;
    y_loc0 = xylim(3) + [ 0.90 0.90 ] * dy;
    % Plot legend manually
    line(ax_res,x_loc,y_loc0,'color',color_order(1,1:3),'linewidth',...
        linewidth_legend);
    text(ax_res,x_loc(2) + 0.005*dx,y_loc0(1), ...
        strcat('${\rm MAP\;simulation}$',sndwch_text),...
        'interpreter','latex','fontsize',fontsize_legend,'color',...
        color_order(1,1:3),'horizontalalignment','left');
    y_loc = y_loc0 - 1/11 * dy;
    % Now add measurement data
    hmult = 400*height_plot;
    plot(ax_res,1/2*sum(x_loc),y_loc(1),'ro','markersize',8,...
        'linewidth',1,'markerfacecolor','r');
    text(ax_res,x_loc(2) + 0.005*dx,y_loc(1)-dy/hmult,...
        '${\rm Measurement\;data}$',...
        'interpreter','latex','fontsize',fontsize_legend,'color',...
        color_order(3,1:3),'horizontalalignment','left',...
        'verticalalignment','middle');
    set(ax_res,'xtick',[],'xticklabels',[],'tickdir','out',...
        'fontsize',fontsize_axis_numbers,...
        'XMinorTick','on','YMinorTick','on','Ticklength',[0.005 0.005]);
    ylh = ylabel(ax_res,'Model output','fontsize',fontsize_labels);
    xylim = [xlim ylim];
    ylh.Position(1) = xylim(1) - 0.03 * (xylim(2) - xylim(1));
    % Add label
    fig_code = strcat('(',ranktoletter(1),')');
    text(ax_res,0.005,0.92,fig_code,'units','normalized',...
        'fontsize',fontsize_A,'interpreter','latex',...
        'horizontalalignment','left');
    % Plot box
    plot_box(ax_res);

    %% B: Now plot raw and/or studentized partial residuals
    ax_res2 = axes('units','normalized');
    ax_res2pos = [ 0.05 0.56 0.9 0.15 ]; set(ax_res2,...
        'position',ax_res2pos); height_plot = 0.15;
    if ismember(DREAMPar.lik,[13 14 16 17 44 45 52])
        yy2ax = eps_n(t_start:t_end);  % studentized partial residuals
    else
        yy2ax = e_n(t_start:t_end);    % studentized raw residuals 
    end
    yyaxis left  % plot raw residuals
    H1 = plot(ax_res2,1:t_end,e); hold on
    yyaxis right % plot studentized partial or studentized raw residuals
    H2 = plot(ax_res2,t_start:t_end,yy2ax);
    set(ax_res2,'fontsize',fontsize_axis_numbers,'tickdir','out');
    set(H1,'linewidth',thickness_line_in_plot,'color',dark_gray,...
        'linestyle','-','marker','none');
    set(H2,'Marker','s','color',light_gray,'linestyle','none',...
        'linewidth',1,'markersize',5,'markerfacecolor',light_gray);
    yyaxis left
    ylh = ylabel(ax_res2,'Residuals','fontsize',fontsize_labels);
    xylim = [xlim ylim];
    ylh.Position(1) = xylim(1) - 0.03 * (xylim(2) - xylim(1));
    set(ax_res2,'ycolor',dark_gray); %;light_gray});
    yyaxis right
    set(ax_res2,'ycolor',light_gray); %;light_gray});
    ylabel(ax_res2,text_axis,'fontsize',fontsize_labels, ...
        'interpreter','latex');
    set(ax_res2,'XLim',[1 t_end],'XMinorTick','on','YMinorTick','on',...
        'Ticklength',[0.005 0.005]);
    % Now add legend - manually to match colors to labels
    xylim = [ get(ax_res2,'XLim') get(ax_res2,'YLim') ];
    dx = (xylim(2)-xylim(1)); dy = (xylim(4)-xylim(3));
    x_loc0 = xylim(1) + [ 0.05 0.07 ] * dx;
    y_loc = xylim(3) + [ 0.13 0.13 ] * dy;
    % Plot 2nd entry of legend manually
    x_loc = x_loc0 + 1/5 * dx;
    % Now add standardized partial residuals
    hmult = 400*height_plot;
    plot(ax_res2,1/2*sum(x_loc),y_loc(1),'ks','markersize',10,...
        'linewidth',2,'color',light_gray,'markerfacecolor',light_gray);
    text(ax_res2,x_loc(2) + 0.002*dx,y_loc(1)-dy/hmult,leg_str2ax,...
        'interpreter','latex','fontsize',fontsize_legend,'color', ...
        light_gray,...
        'horizontalalignment','left','verticalalignment','middle');
    % Plot legend manually
        xylim = [ get(ax_res2,'XLim') get(ax_res2,'YLim') ];
    dx = (xylim(2)-xylim(1)); dy = (xylim(4)-xylim(3));
    x_loc0 = xylim(1) + [ 0.05 0.07 ] * dx;
    y_loc = xylim(3) + [ 0.13 0.13 ] * dy;
    x_loc = x_loc0;
    % Now add standardized partial residuals
    hmult = 400*height_plot;
    line(ax_res2,x_loc,y_loc,'color',dark_gray,'linewidth',...
        linewidth_legend);
    text(ax_res2,x_loc(2)+0.01*dx,y_loc(1)-dy/hmult, ...
        strcat('${\rm Nonstandardized\;residuals}$',sndwch_text),...
        'interpreter','latex','fontsize',fontsize_legend,'color', ...
        dark_gray,...
        'horizontalalignment','left');
    xlabel(ax_res2,'${\rm Observation\;number}$','fontsize', ...
        fontsize_labels,'interpreter','latex');
    %xlh.Position(2) = xylim(3) - 0.26 * (xylim(4)-xylim(3));
    % Add label
    fig_code = strcat('(',ranktoletter(2),')');
    text(ax_res2,0.005,0.85,fig_code,'units','normalized',...
        'fontsize',fontsize_A,'interpreter','latex', ...
        'horizontalalignment','left');
    % yyaxis left; plot_box(ax_res2); yyaxis right; plot_box(ax_res2);
    plot(ax_res2,xylim(2)*[1,1],xylim(3:4),'color','k','linewidth', ...
        ax_res2.LineWidth);
    plot(ax_res2,xylim(1:2),xylim(4)*[1,1],'color','k','linewidth', ...
        ax_res2.LineWidth);

    %% C: Variance of standardized partial or standardized raw residuals
    fig_size = [ 0.18 0.35 ];
    ax_res3 = axes('units','normalized'); ax_res3pos = [ 0.05 0.08 ...
        fig_size ]; set(ax_res3,...
        'position',ax_res3pos); height_plot = fig_size(2);
    if ismember(DREAMPar.lik,[13 14 16 17 44 45])
        plot(ax_res3,FX_MAP(t_start:t_end),eps_n(t_start:t_end),'rs',...
            'color',light_gray,'markersize',5,'linewidth',1,...
            'markerfacecolor',light_gray); hold on;
        c = polyfit(FX_MAP(t_start:t_end),eps_n(t_start:t_end),1);
    else %if ismember(DREAMPar.lik,[11 12 16])
        plot(ax_res3,FX_MAP(t_start:t_end),e_n(t_start:t_end),'rs',...
            'color',light_gray,'markersize',5,'linewidth',1,...
            'markerfacecolor',light_gray); hold on;
        c = polyfit(FX_MAP(t_start:t_end),e_n(t_start:t_end),1);
    end

    set(ax_res3,'fontsize',fontsize_axis_numbers,'tickdir','out',...
        'XMinorTick','on','YMinorTick','on','Ticklength',[0.02 0.02]);
    title(ax_res3,strcat('Stable variance?',sndwch_text),'fontsize', ...
        fontsize_title);
    xlabel(ax_res3,'Model output','fontsize',fontsize_labels);
    ylabel(ax_res3,text_axis,'fontsize',fontsize_labels, ...
        'interpreter','latex');
    axis tight;
    plot(ax_res3,xlim,c(1)*xlim+c(2),'r','linewidth',...
        thickness_line_in_plot,'linestyle','--');

    xylim = [ get(ax_res3,'XLim') get(ax_res3,'YLim') ];
    dx = (xylim(2)-xylim(1)); dy = (xylim(4)-xylim(3));
    x_loc = xylim(1) + [ 0.13 0.21 ] * dx;
    y_loc0 = xylim(3) + [ 0.94 0.94 ] * dy;
    % Plot legend manually
    plot(ax_res3,1/2*sum(x_loc),y_loc0(1),'rs','color',light_gray,...
        'markersize',8,'linewidth',2,'markerfacecolor',light_gray);
    hmult = 400*height_plot;
    text(ax_res3,x_loc(2) + 0.03*dx,y_loc0(1)-dy/hmult,text_axis,...
        'interpreter','latex','fontsize',fontsize_legend,'color', ...
        light_gray,...
        'horizontalalignment','left','verticalalignment','middle');
    y_loc = y_loc0 - 1/12 * dy;
    % Now add least squares line
    line(ax_res3,x_loc,y_loc,'color','r','linewidth',...
        linewidth_legend,'linestyle',':');
    text(ax_res3,x_loc(2) + 0.03*dx,y_loc(1)-dy/hmult, ...
        '${\rm Least\;squares}$',...
        'interpreter','latex','fontsize',fontsize_legend,'color','r',...
        'horizontalalignment','left');
    % Add label
    fig_code = strcat('(',ranktoletter(3),')');
    text(ax_res3,0.02,0.94,fig_code,'units','normalized',...
        'fontsize',fontsize_A,'interpreter','latex',...
        'horizontalalignment','left');
    % Plot box
    plot_box(ax_res3);

    %% D: Histogram of residuals
    nl = min(25,t_end/20);
    % Define axes
    ax_res4 = axes('units','normalized'); ax_res4pos = [ 0.29 0.08 ...
        fig_size ];
    set(ax_res4,'position',ax_res4pos);
    % July 2024: Determine binning method
    switch DREAMPar.lik
        case {13,14,16,17,44,45}
            % Check whether we have an integer range
            res_plot = eps_n(t_start:t_end);
        otherwise
            res_plot = e_n(t_start:t_end);
    end
    % determine binning method
% %     iq_res = iqr(res_plot);
% %     if iq_res > 1 && iq_res < 10
% %         bin_method = 'integers';
% %     else
% %         bin_method = 'auto';
% %     end
    bin_method = 'auto';
    if nl < 1
        % do nothing
    else
        nbins = min(calcnbins(res_plot),maxbins);
        [n_res,edges] = histcounts(res_plot,nbins, ...
            'normalization','pdf','binmethod',bin_method);
        % Test
        n_res = [0 , n_res , 0]; edges = [-inf, edges,inf];
        % July 11, 2024: determine suitable x-ranges
        r_99 = prctile(res_plot,[0.5 99.5]); 
        % Determine midpoints of bins
        x_res = 1/2*(edges(1:end-1) + edges(2:end));
        stem(ax_res4,x_res,n_res,'filled','ro','color',light_gray,...
            'markersize',5.5,'linewidth',2,'displayname',char(text_axis));
        hold on;
        set(ax_res4,'fontsize',fontsize_axis_numbers,'tickdir','out',...
            'Ticklength',[0.02 0.02],'XMinorTick','on','YMinorTick','on');
        title(ax_res4,strcat('Histograms',sndwch_text),'fontsize', ...
            fontsize_title);
        xlabel(ax_res4,text_axis,'fontsize',fontsize_labels); %axis tight;
        ylabel(ax_res4,'Density','fontsize',fontsize_labels);
    end
    if ismember(DREAMPar.lik,[13 14 16 17 44 45])
        data_res = sortrows([eps_n(t_start:t_end) ...
            f_eps_n(t_start:t_end)],1);
    else %if ismember(DREAMPar.lik,[11 12 15 16])
        data_res = sortrows([e_n(t_start:t_end) f_e_n(t_start:t_end)],1);
    end
    plot(ax_res4,data_res(:,1),data_res(:,2),'r-','linewidth',...
        2,'color','r');
    ytickformat(ax_res4,'%2.1f'); % axis tight

    % Now add legend - manually to match colors to labels
    % July 11, 2024: determine suitable x-ranges
    xylim = nan(1,4); xylim(1:2) = r_99; xylim(3:4) = get(ax_res4,'YLim'); 
    axis(xylim);
    % now get ranges again
%    xylim = [ get(ax_res4,'XLim') get(ax_res4,'YLim') ];
    dx = (xylim(2)-xylim(1)); dy = (xylim(4)-xylim(3));
    x_loc = xylim(1) + [ 0.13 0.21 ] * dx;
    y_loc0 = xylim(3) + [ 0.94 0.94 ] * dy;
    line(ax_res4,x_loc,y_loc0,'color',light_gray,'linewidth',...
        linewidth_legend);
    plot(ax_res4,x_loc(2),y_loc0(1),'ro','color',light_gray,...
        'markersize',8,'markerfacecolor',light_gray);
    text(ax_res4,x_loc(2)+0.03*dx,y_loc0(1)-dy/hmult,text_axis,...
        'interpreter','latex','fontsize',fontsize_legend,'color',...
        light_gray,'horizontalalignment','left');
    % Plot 2nd entry of legend manually
    y_loc = y_loc0 - 1/12 * dy;
    line(ax_res4,x_loc,y_loc,'color','r','linewidth',...
        linewidth_legend);
    text(ax_res4,x_loc(2)+0.03*dx,y_loc(1),'${\rm Likelihood}$',...
        'interpreter','latex','fontsize',fontsize_legend,'color','r',...
        'horizontalalignment','left');
    % Add label
    fig_code = strcat('(',ranktoletter(4),')');
    text(ax_res4,0.02,0.94,fig_code,'units','normalized',...
        'fontsize',fontsize_A,'interpreter','latex',...
        'horizontalalignment','left');
    % Plot box
    plot_box(ax_res4);

    %% E: Autocorrelation function of the residuals
    ax_res5 = axes('units','normalized'); ax_res5pos = [ 0.53 0.08 ...
        fig_size ];
    set(ax_res5,'position',ax_res5pos);
    numLags = min(30,t_end-1);
    if ismember(DREAMPar.lik,[13 14 16 17 44 45])
        [acfunc,lags,bnds] = autocorr(eps_n(t_start:t_end),numLags);
    else %if ismember(DREAMPar.lik,[11 12 16])
        [acfunc,lags,bnds] = autocorr(e_n(t_start:t_end),numLags);
    end
    lineHandles = stem(ax_res5,lags,acfunc,'filled','r-o','color',...
        light_gray,'displayname',char(text_axis)); hold on;
    set(lineHandles(1),'MarkerSize',5.5,'linewidth',2); axis tight;
    set(ax_res5,'fontsize',fontsize_axis_numbers,'tickdir','out');
    xlabel(ax_res5,'Lag','fontsize',fontsize_labels);
    axis tight;
    ylabel(ax_res5,'Autocorrelation','fontsize',fontsize_labels);
    title(ax_res5,strcat('Sample ACF',sndwch_text), ...
        'fontsize',fontsize_title);
    set(ax_res5,'Ticklength',[0.02 0.02],'XMinorTick','on', ...
        'YMinorTick','on');
    numMA = 0;
    plot(ax_res5,[numMA+0.5 numMA+0.5; numLags numLags],...
        [bnds([1 1]) bnds([2 2])],'r--','linewidth', ...
        thickness_line_in_plot);
    a = axis; axis(ax_res5,[0 numLags max(-1,min([-0.22 0.9*a(3) ...
        1.1*a(3)])) min(1,max(1.1*a(4),0.9*a(4)))]);
    set(ax_res5,'ytick',-1.0:0.2:1.0,'yticklabels',{'-1.0','-0.8',...
        '-0.6','-0.4','-0.2','0.0','0.2','0.4','0.6','0.8','1.0'});
    % Now add legend - manually to match colors to labels
    xylim = [ get(ax_res5,'XLim') get(ax_res5,'YLim') ];
    dx = (xylim(2)-xylim(1)); dy = (xylim(4)-xylim(3));
    x_loc = xylim(1) + [ 0.13 0.21 ] * dx;
    y_loc0 = xylim(3) + [ 0.94 0.94 ] * dy;
    line(ax_res5,x_loc,y_loc0,'color',light_gray,'linewidth',...
        linewidth_legend);
    plot(ax_res5,x_loc(2),y_loc0(1),'ro','color',light_gray,...
        'markersize',8,'markerfacecolor',light_gray);
    text(ax_res5,x_loc(2)+0.03*dx,y_loc0(1)-dy/hmult,text_axis,...
        'interpreter','latex','fontsize',fontsize_legend,'color',...
        light_gray,'horizontalalignment','left');
    % Plot 2nd entry of legend manually
    y_loc = y_loc0 - 1/12 * dy;
    % Now add bounds
    line(ax_res5,x_loc,y_loc,'color','r','linewidth',...
        linewidth_legend,'linestyle',':');
    text(ax_res5,x_loc(2) + 0.03*dx,y_loc(1)-dy/hmult, ...
        '${\rm 95\%\;Bounds}$',...
        'interpreter','latex','fontsize',fontsize_legend,'color','r',...
        'horizontalalignment','left');
    % Add label
    fig_code = strcat('(',ranktoletter(5),')');
    text(ax_res5,0.02,0.94,fig_code,'units','normalized',...
        'fontsize',fontsize_A,'interpreter','latex',...
        'horizontalalignment','left');
    % Plot box
    plot_box(ax_res5);

    %% G: Now quantile - quantile plot of residuals
    if ismember(DREAMPar.lik,[13 14 16 17 44 45])
        y = sort(eps_n(t_start:t_end));
    else % if ismember(DREAMPar.lik,[11 12 16])
        y = sort(e_n(t_start:t_end));
    end
    % Compute plotting position of numel(y) numbers
    pp = plotpos(1:numel(y)); 
    switch DREAMPar.lik
        case {14,44} % SEP distribution - GL - partial residuals
            x = SEPinv(pp,nuisvar(3),nuisvar(4)); str_leg = 'SEP';
            x_lab = 'SEP quantiles';
        case {16} % Laplace distribution
            x = LAPinv(pp,0,1); str_leg = 'Laplace';
            x_lab = 'Laplace quantiles';
        case {17} % SST distribution
            x = SSTinv(pp,nuisvar(3),nuisvar(4)); str_leg = 'SST';
            x_lab = 'SST quantiles';
        case {45} % SGT distribution
            x = SGTinv(pp,0,1,nuisvar(3),nuisvar(4),nuisvar(5));
            str_leg = 'SGT'; x_lab = 'SGT quantiles';
        otherwise % 11,12,13 and all others
            x = icdf('normal',pp,0,1); str_leg = 'Normal';
            x_lab = 'Normal quantiles';
    end
    [mx,my] = getmxmy(x,y);
    % Define axes
    ax_res6 = axes('units','normalized');
    ax_res6pos = [ 0.77 0.08 fig_size ];
    set(ax_res6,'position',ax_res6pos);
    plot(ax_res6,x,y,'rs','color',light_gray,'markersize',5,...
        'linewidth',1,'markerfacecolor',light_gray); hold on;
    plot(ax_res6,mx,my,'k--','linewidth',thickness_line_in_plot,...
        'color','r','displayname',str_leg); axis tight;
    set(ax_res6,'fontsize',fontsize_axis_numbers,'tickdir','out');
    set(ax_res6,'Ticklength',[0.02 0.02],'XMinorTick','on',...
        'YMinorTick','on');
    xlabel(ax_res6,x_lab,'fontsize',fontsize_labels);
    ylabel(ax_res6,strcat('Quantiles',{' '},lower(text_axis)),...
        'fontsize',fontsize_labels);
    title(ax_res6,strcat('QQ plot',sndwch_text), ...
        'fontsize',fontsize_title);
    % Add legend manually
    xylim = [ get(ax_res6,'XLim') get(ax_res6,'YLim') ];
    dx = (xylim(2)-xylim(1)); dy = (xylim(4)-xylim(3));
    x_loc = xylim(1) + [ 0.13 0.19 ] * dx;
    y_loc0 = xylim(3) + [ 0.94 0.94 ] * dy;
    % Plot legend manually
    plot(ax_res6,1/2*sum(x_loc),y_loc0(1),'rs','color',light_gray,...
        'markersize',8,'linewidth',2,'markerfacecolor',light_gray);
    hmult = 400*height_plot;
    text(ax_res6,x_loc(2) + 0.03*dx,y_loc0(1)-dy/hmult,text_axis,...
        'interpreter','latex','fontsize',fontsize_legend,'color',...
        light_gray,...
        'horizontalalignment','left','verticalalignment','middle');
    y_loc = y_loc0 - 1/12 * dy;
    % Now add least squares line
    line(ax_res6,x_loc,y_loc,'color','r','linewidth',...
        linewidth_legend,'linestyle',':');
    text(ax_res6,x_loc(2) + 0.03*dx,y_loc(1)-dy/hmult,...
        '${\rm Quantile\;line}$',...
        'interpreter','latex','fontsize',fontsize_legend,'color','r',...
        'horizontalalignment','left');
    % Add label
    fig_code = strcat('(',ranktoletter(6),')');
    text(ax_res6,0.02,0.94,fig_code,'units','normalized',...
        'fontsize',fontsize_A,'interpreter','latex',...
        'horizontalalignment','left');
    % Plot box
    plot_box(ax_res6); set(gcf,'color','w');
end

%% DOUBLE CHECK: STILL REQUIRES MORE CHECKING
if exist('FX_MAP','var') %&& (1<1)
    % Transpose fx_mod and fx_par
    fx_mod = fx_mod'; FX_post = FX_post';

    %% Compute scoring rules of forecast density
    % Logarithmic score [Good, 1952]: maximize
    [LS,~,num_zeroLS] = log_score(fx_mod,Meas_info.Y);
    % Continuous Ranked Probability Score [Matheson and Winkler, 1976]: maximize
    CRPS = CRP_score(fx_mod,Meas_info.Y);
    % Spherical score (zeta = 2): maximize
    [SS,~,num_zero] = spherical_score(fx_mod,Meas_info.Y);
    % Dawid-Sebastiani score: maximize
    [DSS,~,~,m_F,s_F] = dawidsebas_score(fx_mod,Meas_info.Y);
    % Interval score: maximize
    [IS,~] = interval_score(tot_unc,Meas_info.Y,p_alfa);

    %% Compute summary metrics of forecast density
    % Compute p-values (CHECK if < or > than obs)
    p_val = p_values(fx_mod,Meas_info.Y);
    % 1. Compute reliability from p-values [Renard et al., 2011]: maximize
    [RLBL,eCDF,uCDF] = rlbl(p_val);
    % 2. Compute coefficient of variation [Evin et al., 2013]
    CV = 1/Meas_info.n * sum(s_F./m_F);
    % 3a. Mean spread/width of 100alfa confidence interval
    width_par = mean(par_unc(:,2) - par_unc(:,1));
    % 3b. Mean spread/width of 100alfa prediction interval
    width_mod = mean(tot_unc(:,2) - tot_unc(:,1));
    % 4a. Mean precision of parameters [McInerney et al., 2017]
    prec_par = mean(std(FX_post,0,2))./mean(Meas_info.Y);
    % 4b. Mean precision of total [McInerney et al., 2017]
    prec_mod = mean(std(fx_mod,0,2))./mean(Meas_info.Y);
    % 5a. Mean percentage bias [McInerney et al., 2017]
    pbias_par = 100 * sum( mean(FX_post,2) - Meas_info.Y ) / sum(Meas_info.Y);
    %  pbias_par = 100*mean(sum(bsxfun(@minus,fx_post,Meas_info.Y))./sum(Meas_info.Y));
    % 5b. Mean percentage bias [McInerney et al., 2017]
    pbias_mod = 100 * sum( mean(fx_mod,2) - Meas_info.Y ) / sum(Meas_info.Y);
    %  pbias_mod = 100*mean(sum(bsxfun(@minus,fx_mod,Meas_info.Y))./sum(Meas_info.Y));
    % 6a. Mean log-density of posterior realizations
    meanLogP = mean(sum(P_post(:,DREAMPar.d+1:DREAMPar.d+2),2));
    % 6b. Std deviation log-density of posterior realizations
    stdLogP = std(sum(P_post(:,DREAMPar.d+1:DREAMPar.d+2),2));

    %% Compute summary metrics of MAP solution
    % 1a. Max of log-likelihood 
    maxLogL = max(P_post(:,DREAMPar.d+2));
    % 1b. Max of log posterior density
    maxLogP = max(sum(P_post(:,DREAMPar.d+1:DREAMPar.d+2),2));

    %% PLOT SOME OF THE RESULTS
    %% A. Quantile - quantile plot of empirical and theoretical CDF of p-values
    figure('units','normalized','outerposition',[0.05 0.05 s0(4)/s0(3) 1]);
    % Now set the axes
    ax_quant = axes('units','normalized'); ax_quantpos = [ 0.1 0.1 0.8 0.8 ]; 
    set(ax_quant,'position',ax_quantpos);
    ax_quant.ColorOrder = color_order;
    plot(ax_quant, eCDF, uCDF, 'color',light_gray,'displayname',str_leg);
    hold on
    set(ax_quant,'fontsize',fontsize_axis_numbers,'tickdir','out', ...
        'XMinorTick','on','YMinorTick','on');
    line(ax_quant,[0 1],[0 1],'color','r','linestyle','--','linewidth',2, ...
        'displayname','1:1 line'); % 1:1 line
    %    href.Color = [1 0 0];
    ylabel(ax_quant,'Empirical quantiles','interpreter','latex', ...
        'fontsize',fontsize_labels);
    xlabel(ax_quant,'Theoretical quantiles','interpreter','latex', ...
        'fontsize',fontsize_labels);
    title(ax_quant,strcat(['Quantile-quantile plot of empirical and ' ...
        'theoretical CDF of p-values'],sndwch_text),...
        'interpreter','latex','fontsize',fontsize_title);
    [h_leg7,icons7] = legend(ax_quant,'show','location','best', ...
        'box','off','fontsize',fontsize_legend);
    icons7 = findobj(icons7,'Type','line');
    h_leg7.String{1} = ['\color[rgb]{' num2str(light_gray) '} ' h_leg7.String{1}];
    h_leg7.String{2} = ['\color[rgb]{' num2str([1 0 0]) '} ' h_leg7.String{2}];
    set(icons7(1:2),'linewidth',3); set(icons7(3:4),'markersize',10, ...
        'linewidth',2);
    %set(ax_quant,'TickLabelInterpreter','latex','fontsize',fontsize_numbers);
    % Plot box
    plot_box(ax_quant); set(gcf,'color','w');

    %% B. Rank Histogram --> want horizontal line (necessary but not
    % sufficient condition)
    %%[rank_hist,~] = Rank_hist(fx_mod,Meas_info.Y);
    %% Can get rank histogram from p_values: rank_hist =
    % ceil(p_val*size(fx_mod,2)) + 1? (CHECKED: APPROX. TRUE)
    rank_hist = ceil(p_val*size(fx_mod,2)) + 1;
    [nr,xr] = histcounts(rank_hist,'normalization','count');
    % Bin center and normalize with upper bin value to get Xnorm in [0,1]
    Xr = 1/2*(xr(2:end) + xr(1:end-1)); Xr = Xr/max(xr);
    figure('units','normalized','outerposition',[0.05 0.05 s0(4)/s0(3) 1]);
    ax_bins = axes('units','normalized'); ax_binspos = [ 0.1 0.1 0.8 0.8 ];
    set(ax_bins,...
        'position',ax_binspos);
    ax_bins.ColorOrder = color_order;
    bar(ax_bins,Xr,nr/trapz(Xr,nr));
    %bar(ax_bins,Xr,nr);
    xlabel(ax_bins,'Scaled rank of observation','fontsize',fontsize_labels);
    ylabel(ax_bins,'Probability density','fontsize',fontsize_labels);
    set(ax_bins,'fontsize',fontsize_axis_numbers,'tickdir','out', ...
        'XMinorTick','on','YMinorTick','on'); ytickformat(ax_bins,'%3.1f');
    title(ax_bins,strcat(['Rank histogram of the observations given PDF ' ...
        'of total uncertainty'],sndwch_text),'fontsize',fontsize_title);
    % Plot box
    plot_box(ax_bins); set(gcf,'color','w');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    switch sndwch
        case 0
            %% Now print results
            file_metric = strcat(method_old,['_forecast_' ...  % Default name of metrics file
                'metrics.txt']);  
            file_figs = strcat(method_old,'_figures.pdf');     % Default name of figures file
            file_data = strcat(method_old,'_data');            % Default name of data file
        case 1
            %% Now print results
            file_metric = strcat(method_old,['_forecast_' ...  % Default name of metrics file
                'metrics_sandwich.txt']);  
            file_figs = strcat(method_old,['_figures_' ...     % Default name of figures file
                'sandwich.pdf']);             
            file_data = strcat(method_old,'_data_sandwich');   % Default name of data file
    end
    if Lik_info.filename % Filename specified by user in global variable LV
        % Create name of output here
        name_handle = char(f_handle); ii2 = strfind(name_handle,'(x');
        name_mod = name_handle(5:ii2(2)-1);
        str_prt_nuis = [];
        for z = 1:numel(Lik_info.str_nuis)
            str_prt = char(Lik_info.str_nuis(z));
            str_prt_nuis = strcat(str_prt_nuis,'_',str_prt(2:end-1));
        end
        % Remove latest backslash
        if ~isempty(str_prt_nuis)
            str_prt_nuis = erase(str_prt_nuis,'\');
        end
        % Create filename: combine model name, likelihood, nuisance vars
        file_name = strcat(name_mod,'_',Lik_info.name_lik_func, ...
            str_prt_nuis);
        switch sndwch
            case 0        
                % Define name of metrics file
                file_metric = strcat(file_name,['_forecast_metrics_' ...
                    'normalized.txt']);
                % Define name of figures file
                file_figs = strcat(file_name,'_figures_normalized.pdf');
                % Define name of data file
                file_data = strcat(file_name,'_data_normalized.mat');
            case 1
                % Define name of metrics file
                file_metric = strcat(file_name,['_forecast_metrics_' ...
                    'normalized_sandwich.txt']);
                % Define name of figures file
                file_figs = strcat(file_name,['_figures_normalized_' ...
                    'sandwich.pdf']);
                % Define name of data file
                file_data = strcat(file_name,['_data_normalized_' ...
                    'sandwich.mat']);
        end                
    end
    % Print Table to file
    fid_latex = fopen(file_metric,'w');
    fprintf(fid_latex,['=================================================' ...
        '=================================================\n']);
    fprintf(fid_latex,' Likelihood %d \n',DREAMPar.lik);
    switch sndwch
        case 0
            fprintf(fid_latex,' Sandwich correction is INACTIVE \n');
        case 1
            fprintf(fid_latex,' Sandwich correction is ACTIVE \n');
    end            
    fprintf(fid_latex,' %s \n',Lik_info.name_lik);
    if ~isempty(Lik_info.str_nuis)
        str_pr = char(Lik_info.str_nuis(1)); fprintf(fid_latex,[' ' ...
            'Nuisance variables: %s'],str_pr(2:end-1));
        for z = 2:numel(Lik_info.str_nuis)
            str_pr = char(Lik_info.str_nuis(z)); fprintf(fid_latex,', %s ', ...
                str_pr(2:end-1));
        end
        fprintf(fid_latex,'\n');
    end
    fprintf(fid_latex,[' Total number of observations/simulated output : ' ...
        '%d \n'],Meas_info.n);
    fprintf(fid_latex,['__________________________________________________' ...
        '________________________________________________\n']);
    fprintf(fid_latex,[' Score rules                        Abbrev.       ' ...
        '   Value   Unit        Reference                \n']);
    fprintf(fid_latex,['__________________________________________________' ...
        '________________________________________________\n']);
    fprintf(fid_latex,[' Logarithmic Score                  LS          ' ...
        '%10.3f♪  [log(1/y)]  Good 1952                  \n'],LS);
    fprintf(fid_latex,[' Continuous Rnkd Probability Score  CRPS        ' ...
        '%10.3f   [y]         Matheson & Winkler 1972  \n'],CRPS);
    fprintf(fid_latex,[' Spherical Score                    SS          ' ...
        '%10.3f   [-]         Good 1971, Friedman 1983 \n'],SS);
    fprintf(fid_latex,[' Dawid-Sebastiani Score             DSS         ' ...
        '%10.3f   [log(y^2)]  Dawid & Sebastiani 1999  \n'],DSS);
    fprintf(fid_latex,[' Interval Score                     IS          ' ...
        '%10.3f   [y]         Gneiting & Raftery 2007  \n'],IS);
    fprintf(fid_latex,['__________________________________________________' ...
        '________________________________________________\n']);
    fprintf(fid_latex,[' Summary metrics density forecast   Abbrev.       ' ...
        '   Value   Unit        Reference                \n']);
    fprintf(fid_latex,['_________________________________________________' ...
        '_________________________________________________\n']);
    fprintf(fid_latex,[' Reliability                        RLBL        ' ...
        '%10.3f   [-]         Renard et al. 2011         \n'],RLBL);
    fprintf(fid_latex,[' Coefficient of variation           CV          ' ...
        '%10.3f   [-]         Evin et al. 2013           \n'],CV);
    fprintf(fid_latex,[' Containment of 95%% tot interval    ctnd_tot    ' ...
        '%10.3f   [%%]                                     \n'],ctnd_tot);
    fprintf(fid_latex,[' Mean tot. width (spread at 95%%)    wdth_mod    ' ...
        '%10.3f   [y]♫                                    \n'],width_mod);
    fprintf(fid_latex,[' Mean par. width (spread at 95%%)    wdth_par    ' ...
        '%10.3f   [y]♫                                    \n'],width_par);
    fprintf(fid_latex,[' Mean par. precision                prec_par    ' ...
        '%10.3f   [-]         McInerney et al. 2017      \n'],prec_par);
    fprintf(fid_latex,[' Mean tot. precision                prec_mod    ' ...
        '%10.3f   [-]         McInerney et al. 2017      \n'],prec_mod);
    fprintf(fid_latex,[' Mean percentage bias par           pbias_par   ' ...
        '%10.3f   [%%]         McInerney et al. 2017     \n'],pbias_par);
    fprintf(fid_latex,[' Mean percentage bias tot           pbias_mod   ' ...
        '%10.3f   [%%]         McInerney et al. 2017     \n'],pbias_mod);
    fprintf(fid_latex,[' Mean log-posterior density         meanLogP    ' ...
        '%10.3f   [-]                                      \n'],meanLogP);
    fprintf(fid_latex,[' St. dev. of log-posterior density  stdLogP     ' ...
        '%10.3f   [-]                                      \n'],stdLogP);
    fprintf(fid_latex,['_________________________________________________' ...
        '_________________________________________________\n']);
    fprintf(fid_latex,[' Summary metrics of MAP solution    Abbrev.       ' ...
        '   Value   Unit        Reference                \n']);
    fprintf(fid_latex,['_________________________________________________' ...
        '_________________________________________________\n']);
    fprintf(fid_latex,[' Maximum log-posterior density      maxLogP     ' ...
        '%10.3f   [-]                                      \n'],maxLogP);
    fprintf(fid_latex,[' Maximum log-likelihood             maxLogL     ' ...
        '%10.3f   [-]                                      \n'],maxLogL);
    fprintf(fid_latex,[' Root mean square error MAP output  RMSE_MAP    ' ...
        '%10.3f   [y]♫                                     \n'],RMSE_MAP);
    if ~exist('PBIAS_MAP','var')
        PBIAS_MAP = nan;
    end
    fprintf(fid_latex,[' Percentage bias MAP output         PBIAS_MAP   ' ...
        '%10.3f   [%%]                                     \n'],PBIAS_MAP);
    fprintf(fid_latex,['=================================================' ...
        '=================================================\n']);
    fprintf(fid_latex,[' ♪ Of the %d observations %s have density less ' ...
        'than real_min (= 2.22e-308) per forecast pdf \n'],...
        Meas_info.n,num2words(num_zeroLS));
    fprintf(fid_latex,[' ♫ Unit of observations in training data record ' ...
        '\n']);
    fclose(fid_latex);
    % Now open file in text editor of MATLAB
    edit(file_metric)
    % Print Table
    col_val = [LS;CRPS;SS;DSS;IS;RLBL;CV;ctnd_tot;width_mod;width_par; ...
        prec_par;prec_mod;pbias_par;pbias_mod;meanLogP;stdLogP;maxLogP; ...
        maxLogL;RMSE_MAP;PBIAS_MAP];
    col_val_ab = num2str(col_val,'%5.3f');
    tabout = table(char({'Logarithmic Score'; ...
        'Continuous Ranked Probability Score'; ...
        'Spherical Score'; ...
        'Dawid-Sebastiani Score'; ...
        'Interval Score'; ...
        'Reliability'; ...
        'Coefficient of variation'; ...
        'Percentage contained in 95% total uncertainty'; ...
        'Width of 95% total uncertainty'; ...
        'Width of 95% parameter uncertainty'; ...
        'Precision of 95% parameter uncertainty'; ...
        'Precision of 95% total uncertainty'; ...
        'Percentage bias of 95% parameter uncertainty'; ...
        'Percentage bias of 95% total uncertainty'; ...
        'Mean log-posterior density'; ...
        'Std. dev. of log-posterior density'; ...
        'Maximum log-posterior density'; ...
        'Maximum log-likelihood'; ...
        'Root Mean Square Error of MAP solution'; ...
        'Percentage bias of MAP solution'}), ...
        char(col_val_ab),...
        char({'[log(1/y)] ♪';'[y]';'[-]';'[log(y^2)] ♫';'[y] ♫';'[-]'; ...
        '[-]';'[%]';'[y] ♫';'[y] ♫';'[-]';...
        '[-]';'[%]';'[%]';'[-]';'[-]';'[-]';'[-]';'[y] ♫';'[%]'}),...
        char({'Good (1952)';'Matheson and Winkler (1972)';'Good (1971)';...
        'Dawid and Sebastiani (1999)';'Gneiting and Raftery (2007)';...
        'Renard et al. (2010)';'Evin et al. (2013)';'-';...
        '-';'-';'McInerney et al. (2017)';...
        'McInerney et al. (2017)';'McInerney et al. (2017)';...
        'McInerney et al. (2017)';'-';'-';'-';'-';'[y] ♫';'-'}),...
        'VariableNames',{'Metric','Value','Unit','Reference'},...
        'RowNames',{'A';'B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';'M';...
        'N';'O';'P';'Q';'R';'S';'T'});
    % Print table to command window and add footnotes
    disp(tabout)
    fprintf(['    ♪ Of the %d observations %s have a density less than ' ...
        'real_min (= 2.22e-308) according to forecast pdf \n'], ...
        Meas_info.n,num2words(num_zero));
    fprintf('    ♫ Unit of observations in training data record \n');
    % Now print table to figure (cannot use uifigure as does not work
    % with pdf)
    fig = figure('Position',[300 300 920 600]);
    ui_table = uitable(fig,'Data',table2cell(tabout));
    % original column width: {500,150,120,300}
    set(ui_table,'InnerPosition',[0 0 100 100],'OuterPosition', ...
        [0 0 1070 600],'Fontsize',14,'Columnwidth',{400,100,120,300});
    % Removed column with variable name used in code
    ui_table.ColumnName = tabout.Properties.VariableNames;
    %renaming columns to that of table
    ui_table.RowName = [];
    % Now print this to PDF
    print(fig,'-dpdf',file_figs,'-fillpage');

    % Save the memory of MATLAB - only in case LV.filename specified
    if Lik_info.filename
        % MATLAB stop - no file is small?
        eval(char(strcat('save',{' '},file_data,{' '},['-v7.3 chain ' ...
            'output MAP_un FX_MAP DREAMPar Par_info Meas_info ' ...
            'Lik_info options MAP_info par_unc tot_unc'])));
        % Write a Latex table as well with results
        fid_latex = fopen('latex_table.txt','a+','n');
        % Print name of likelihood
        fprintf(fid_latex,' %s & ',Lik_info.name_lik_func);
        % Print nuisance variables - user cannot select any
        if ~isempty(Lik_info.str_nuis)
            str_pr = char(Lik_info.str_nuis(1)); fprintf(fid_latex,'$%s', ...
                str_pr(2:end-1));
            for z = 2:numel(Lik_info.str_nuis) - 1
                str_pr = char(Lik_info.str_nuis(z)); fprintf(fid_latex, ...
                    ', %s ', str_pr(2:end-1));
            end
            str_pr = char(Lik_info.str_nuis(numel(Lik_info.str_nuis)));
            fprintf(fid_latex,', %s$',str_pr(2:end-1));
        else
            fprintf(fid_latex,', %s$',[]);
        end
        % Npw print the various scoring rules/metrics
        fprintf(fid_latex,[' & %10.3f & %10.3f & %10.3f & %10.3f & ' ...
            '%10.3f && %10.3f & %10.3f & %10.3f & %10.3f && %10.3f ' ...
            '& %10.3f & %10.3f \\\\ \n'],...
            LS,CRPS,SS,DSS,IS,RLBL,CV,ctnd_tot/100,width_mod, ...
            maxLogP,RMSE_MAP,PBIAS_MAP);
        % Now close file
        fclose(fid_latex);
    end

    % MATLAB stop
% %     % Get all figures
    figHandles = findall(0,'Type','figure');
% %     % Save files to file name
% %     switch verLessThan('matlab','9.11')
% %         case 1 % Cannot use append option in exportgraphics
% %             % Add path - may change for others
% %             curr_dir = pwd;
% %             idx = strfind(curr_dir,'\');
% %             if isempty(idx) % may not be working in windows
% %                 idx = strfind(curr_dir,'/');
% %             end
% %             new_dir1 = strcat(curr_dir(1:idx(end)),'export');
% %             addpath(new_dir1);
% %             new_dir2 = strcat(curr_dir(1:idx(end)),'num2words');
% %             addpath(new_dir2);
% %             % Save one figure after the next in PDF file
% %             evalstr = strcat('export_fig -pdf -append',{' '},file_figs);
% %             for zz = 2:numel(figHandles)
% %                 figure(figHandles(zz)); eval(char(evalstr));
% %             end
% %         case 0 % Can use append in exportgraphics (MATLAB, 2021+)
            for zz = 2:numel(figHandles)
                figure(figHandles(zz)); set(gcf,'color','w');
                exportgraphics(figHandles(zz), file_figs, 'Append', true);
            end
% %     end

end
%% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

end