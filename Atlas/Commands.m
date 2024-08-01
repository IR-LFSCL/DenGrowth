% Commands to be run
time0=clock();

OutputFig('Data',[200,200],'Atlas','phi','phi',[0,1],[0,10000,50000]);
OutputFig('Data',[200,200],'Atlas','cZn','cZn',[0,0.7],[0,10000,50000]);
OutputFig('Data',[200,200],'Atlas','cNa','cNa',[0,0.3],[0,10000,50000]);
OutputFig('Data',[200,200],'Atlas','pot','pot',[-0.35,0],[0,10000,50000]);
% 
% con_time = etime(clock(),time0);
% fprintf('Time Consumed: %10d\n',con_time);

