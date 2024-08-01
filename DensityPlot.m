% Density plot of 2-D arrays

function []=DensityPlot(field,ran,step)
data=evalin('base',field);
fig_name=[field,'_',num2str(step)];
figure(1);
surf(data,data,'EdgeColor','none')
axis equal off; caxis(ran); colorbar; view([0,90]);
set(gcf,'PaperPositionMode','manual',...
    'PaperUnits','points','PaperPosition',...
    [0 0 360 360],'Name',fig_name);
end
