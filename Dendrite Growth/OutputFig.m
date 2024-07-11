% Generating figures of one field variable
% folder1: data folder
% dim: data size
% folder2: figure folder
% field: field names
% cr: color range
% r: numerical range to be processed

function []=OutputFig(folder1,dim,folder2,field1,field2,cr,r,k0)
mainf=cd;

for n=r(1):r(2):r(3)
    %1 Reading data from data_folder
    ss=num2str(n);
    data_name=[field1,'_',ss,'.txt'];
    cd([mainf,'/',folder1]);
    data_id=fopen(data_name,'r');
    data=k0*fscanf(data_id,'%f',dim);
    fclose(data_id);
        
    %2 Outputing figure to figure_folder
    fig_name=[field2,'_',ss];
    surf(data,data,'EdgeColor','none')
    axis equal off; caxis(cr); colorbar; view([0,90]);
    set(gcf,'PaperPositionMode','manual',...
        'PaperUnits','points','PaperPosition',...
        [0 0 360 360],'Name',fig_name);
    cd([mainf,'/',folder2]);
    print(gcf,'-dtiff',fig_name);
    cd(mainf);
end
end
