% Draw from prior ranges
N = 1000;
eta_min = ones(3,1) * 1e-5; eta_max = [ 50 50 2 ]';
% Create initial parameter vectors
Eta = repmat(eta_min,1,N) + lhsdesign(d,N) .* ( repmat(eta_max,1,N) - repmat(eta_min,1,N) ); Eta = Eta';
% Run haverkamp_t for exact solution
plugin.I = [1:1:48]'; err = nan(N,1); flag = nan(N,1);
for i = 1:N
    % Solve Haverkamp explicit equation
    [t,~,flag(i)] = Haverkamp_t(Eta(i,1:3),plugin);
    % Now do haverkamp with newton's method
    %Isec = haverkamp(eta(i,1:3),t);
    plugin.t = t;
    % Solve Haverkamp implicit equation
    I_new = Haverkamp_I(Eta(i,1:3),plugin);
    % Compute error
%    err(i,1) = sum(abs(Isec-plugin.I));
    err(i,1) = sum(abs(I_new - plugin.I));
end
% check which ones are larger than say 0.01
idx = find(err > 0.01);
% How many have flag equals 1
idx_1 = find(flag == 1); 
% Make histogram of those solution
figure(100),hist(log10(err(idx_1)));
% In panel add number of solutions with flag 1
% Of the flag 2 solutions of time, which ones are flag 1 solutions of inf?
idx_2 = find(flag == 2);
F = ones(N,1);
plugin.t = [0.1:0.1:240];
for i = 1:numel(idx_2)
    I_new = Haverkamp_I(Eta(i,1:3),plugin);
    if isnan(I_new)
        F(idx_2(i)) = 2;
    end
end % None

% Next, 

% % % % analyze those curves
% % % for i = 1:numel(idx)
% % %     t = haverkamp_t(eta(idx(i),1:3),plugin);
% % %     % Now do haverkamp
% % %     I = haverkamp(eta(idx(i),1:3),t);
% % %     % Plot
% % %     figure(i),plot(t,plugin.I,'k',t,I,'ro'); 
% % %     % Calculate error
% % % %    err(i,1) = sum(abs(I-plugin.I));
% % % end

% % % % Now check solutions that do not work for t as exp term to infinity
% % % idx = find(isnan(err)); 
% % % % analyze those curves
% % % for i = 1:numel(idx)
% % %     t = haverkamp_t(eta(idx(i),1:3),plugin);
% % %     % Now do haverkamp
% % %     I = haverkamp(eta(idx(i),1:3),t);
% % %     % Plot
% % %     figure(i),plot(t,plugin.I,'k',t,I,'ro'); 
% % %     % Calculate error
% % % %    err(i,1) = sum(abs(I-plugin.I));
% % % end
% % %  
% % % 
    