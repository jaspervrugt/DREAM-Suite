function y = ar2_model(pars,plugin)
% Returns zero value

y = zeros(plugin.N,1);

% % if pars(1) + pars(2) < 1
% %     y(1:2,1) = plugin.IC;
% %     for t = 3:plugin.N
% %         y(t,1) = pars(1) * y(t-1,1) + pars(2) * y(t-2,1);
% %     end
% % else
% %     y = zeros(plugin.N,1);
% % end
