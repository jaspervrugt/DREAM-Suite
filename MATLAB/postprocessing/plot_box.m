function plot_box(ax)
% Plot the box around the figure so that we do not have ticks on right and
% upper side
box off
isholdonque = ishold;
hold on
xylim = [ get(ax,'XLim') get(ax,'YLim') ];
plot(xylim(2)*[1,1],xylim(3:4),'color','k','linewidth',ax.LineWidth);
plot(xylim(1:2),xylim(4)*[1,1],'color','k','linewidth',ax.LineWidth);
if isholdonque == 0
    hold off
end
