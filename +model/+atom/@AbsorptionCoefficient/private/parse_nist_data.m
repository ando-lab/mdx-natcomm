% parse material properties

datadir = 'nist_data';
filesToProcess = dir(fullfile(datadir,'*.txt'));

for j=1:numel(filesToProcess)
[~,elname] = fileparts(filesToProcess(j).name);
str = fileread(fullfile(datadir,filesToProcess(j).name));

numstr = '[\d\.E-\+]+';
newlinestr = '(?<=(^|\n))';
searchstr = [newlinestr,'\s*',...
    '(?<edge>[KLM]\d*)?','\s*',...
    '(?<energy>',numstr,')\s+',...
    '(?<mu_abs>',numstr,')\s+',...
    '(?<mu_en>',numstr,')'];

dataTable = struct2table(regexp(str,searchstr,'names'));

% put in matlab table format
fprintf(1,'nist_data_%s = ',elname);
fprintf(1,'table({''%s''},',strjoin(dataTable.edge,''';'''))
fprintf(1,'[%s],',strjoin(dataTable.energy,';'))
fprintf(1,'[%s],',strjoin(dataTable.mu_abs,';'))
fprintf(1,'[%s],',strjoin(dataTable.mu_en,';'))
fprintf(1,'''VariableNames'',{''edge'',''energy'',''mu_atten'',''mu_en''});\n');

end
