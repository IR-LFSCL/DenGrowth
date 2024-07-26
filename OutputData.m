% Outputing data at specified time level

% folder: folder name to store data
% field: field names used as file subname 1
% step: time step used as file subname 2

function []=OutputData(folder,field,step)
dim=size(field);

% Setting data folder as current folder
name_m=cd;
cd([name_m,'/',folder]);

% Outputing all fields
for nf=1:dim(2)
    file_name=[field{nf},'_',num2str(step),'.txt'];
    file_id=fopen(file_name,'w');
    data=evalin('base',field{nf});
    fprintf(file_id,'%+.16d\n',data);
    fclose(file_id);
end        
       
% Setting main folder as current folder
cd(name_m);
end
