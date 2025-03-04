%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                           %
% ---------------------------- Check the following paper ---------------------------------- %
%                                                                                           %
%   Minasny, B., J.A. Vrugt, and A.B. McBratney (2011), Confronting uncertainty in model-   %
%       based geostatistics using Markov chain Monte Carlo simulation, Geoderma, 163,       %
%       150-162, doi:10.1016/j.geoderma.2011.03.011.                                        %
%                                                                                           %
% ----------------------------------------------------------------------------------------- %
%                                                                                           %
% NOTE: Please check output files and figures from built-in postprocessor for results of    %
%       inference                                                                           %
% NOTE: This postprocessor generates an additional plot of the final variogram, its 95%     %
%       uncertainty intervals and the corresponding data                                    %
%                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First assemble all chain in one matrix
ParSet = genparset(chain); DREAMPar.N = size(chain,3);

% Find the maximum aposteriori parameter values (last column of ParSet are log-density values!)
idx = find(ParSet(:,end) == max(ParSet(:,end))); idx = idx(1);

% Take the last 25% of the posterior samples -- assume that these samples
% are posterior samples (double check that R_stat < 1.2 for all parameters)
Pars = ParSet ( floor ( 0.75 * size(ParSet,1) ) : size(ParSet,1), 1 : DREAMPar.d );

% NOTE: Please check results of regular postprocessor
fprintf('POSTPROC_VARIOGRAM: PLEASE RESORT TO TABLES AND FIGURES OF DREAM POSTPROCESSOR FOR RESULTS INFERENCE\n');
fprintf('POSTPROC_VARIOGRAM: THIS POSTPROCESSOR PLOTS THE FINAL VARIOGRAM AND ITS 95%% PREDICTION INTERVALS\n');

% ------------------------------------------------------------------------------------------------------------
% -------------------------------- PLOT THE 95% POSTERIOR SIMULATION UNCERTAINTY -----------------------------
% ------------------------------------------------------------------------------------------------------------

% Load the data
load 'forest.txt';
% Define data to be forest
data = forest;
% Define coordinates
x = data(:,1); y = data(:,2); z = data(:,3);

% Calculate the empirical variogram
[H,G,semv,vcloud] = calc_vario(x,y,z);

% How many posterior samples do we like?
N = size(Pars,1);

% Set the separation distance, hd
hd = (0 : 2 : 800);

% Identify unique posterior solutions
[B,I,J] = unique(Pars(:,1:DREAMPar.d),'rows');
% How many rows does B have?
N_p = size(B,1);
% Loop over each sample
for j = 1 : N_p
    V(j,:) = vmatern(B(j,2),B(j,3),B(j,4),B(j,5),hd);
end
% And create full matrix sim_out
V = V(J,:);

% Calculate the 2.5 percentile posterior variogram
V_low = prctile(V,2.5);

% Calculate the 97.5 percentile posterior variogram
V_high = prctile(V,97.5);

% Now plot the R_statistic for each parameter
figure('units','normalized','outerposition',[0 0 1 1])
% Now plot
plot(hd,mean(V), 'r-', 'color',[0.5 0.5 0.5],'linewidth',3); hold on;
plot(hd,V_low, 'k-.','linewidth',3); plot(hd,V_high, 'k-.','linewidth',3);
plot(semv(:,1),semv(:,2),'r.','markersize',20); set(gca,'fontsize',16);
% Add labels
xlabel('Separation distance (lag) [m]','fontsize',18,'interpreter','latex');
ylabel('Semivariance','fontsize',18,'interpreter','latex');
title('POSTERIOR VARIOGRAM PREDICTION UNCERTAINTY AND EMPIRICAL VARIOGRAM','fontsize',18,'interpreter','latex');
% Add legends
[~,objh] = legend({'mean variogram','2.5\% percentile','97.5\% percentile','empirical variogram'},'interpreter','latex','location','northwest'); legend boxoff; set(objh,'linewidth',3);

% ------------------------------------------------------------------------------------------------------------
% ------------------------------ END PLOT THE 95% POSTERIOR SIMULATION UNCERTAINTY ---------------------------
% ------------------------------------------------------------------------------------------------------------