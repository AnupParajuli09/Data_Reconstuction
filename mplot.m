function mplot(x1,y1,x2,y2)
x_label = 'Bus number'; % x axis label
y_label = 'Power (kW)'; % y axis label
legend_name = {'Original','Reconstruct'}; % legend names

figure('Renderer', 'painters', 'Position', [10 10 1000 500])
plot(x1,y1,'-ob','LineWidth',1.5)
hold on
plot(x2,y2,'-xr','LineWidth',1.5)
xlabel(x_label,'FontSize',18,'FontName','Times New Roman')
ylabel(y_label,'FontSize',18,'FontName','Times New Roman')
title('Reconstructed power data for dissimilar load profile','FontSize',18,'FontName','Times New Roman')
xlim([1 128])
legend (legend_name,'Location','northeast')
set(gca,'fontsize',16,'Fontname','Calibri','GridAlpha',0.5)
ax = gca

ax.XRuler.Axle.LineWidth = 1.5;
ax.YRuler.Axle.LineWidth = 1.5;
grid
grid minor
% legend (legend_name,'Location','southeast')
saveas(gca,'plot.svg')
end