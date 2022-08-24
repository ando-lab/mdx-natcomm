function [LogTables,other_text] = read_loggraph(fn)

theFile = readin(fn);
[parse_lines,other_text] = findTables(theFile);
numTables = length(parse_lines);

LogTables = struct('title',{},'type',{},'Graphs',{},'columns',{},...
    'comments',{},'data',{});

for k=1:numTables
    thisTable = parse_lines(k);
    colNames = strsplit(thisTable.column_lines);
    nColumns = length(colNames);
    
    % parse data
    matchstr_nan = '(?<=\s|^)[^\d\s](?=(\s|$))'; % any single non-digit
    data_lines = regexprep(thisTable.data_lines,matchstr_nan,'NaN');
    data = reshape(sscanf(data_lines,'%f'),nColumns,[])';
    
    % parse graphs
    TheseGraphs = parseGraphs(thisTable.graph_lines);
    
    % assign outputs
    LogTables(k).title = thisTable.table_lines;
    LogTables(k).type = lower(thisTable.plottype);
    LogTables(k).Graphs = TheseGraphs;
    LogTables(k).columns = colNames;
    LogTables(k).comments = thisTable.comment_lines;
    LogTables(k).data = data;
    
end
end

function theLines = readin(fn)
fid = fopen(fn,'r');
theLines = textscan(fid,'%s','delimiter','\n');
fclose(fid);
% join the lines again
theLines = strjoin(theLines{1},'\n');
end

function [parse_lines,other_text] = findTables(theFile)
matchstr_table_lines = '\$TABLE\s*:\s*(?<table_lines>[^:]*?)\s*:\s*';
matchstr_graph_lines = '\$(?<plottype>(GRAPHS|SCATTER))\s*:\s*(?<graph_lines>.*?)\s*\$\$';
matchstr_column_lines = '\s*(?<column_lines>.*?)\s*\$\$';
matchstr_comment_lines = '\s*(?<comment_lines>.*?)\s*\$\$';
matchstr_data_lines = '\s*(?<data_lines>.*?)\s*\$\$';

matchstr = [matchstr_table_lines,matchstr_graph_lines,...
    matchstr_column_lines,matchstr_comment_lines,matchstr_data_lines];

[parse_lines,other_text] = regexp(theFile,matchstr,'names','split');
end

function Graphs = parseGraphs(graph_lines)

matchstr_graph = ['\s*(?<graph_name>[^:]*?)',...
    '\s*:\s*',...
    '(?<graph_type>[^:]*?)',...
    '\s*:\s*',...
    '(?<column_list>[\d\,]*?)',...
    '\s*:\s*'];

parse_graphs = regexp(graph_lines,matchstr_graph,'names');

numGraphs = length(parse_graphs);

Graphs = struct('name',{},'cols',{},'axisMode',{},'axisLimits',{});

for j=1:numGraphs
    graphType = parse_graphs(j).graph_type;
    cols = str2num(parse_graphs(j).column_list);
    graphName = parse_graphs(j).graph_name;
    axisLimits = [];
    
    switch lower(graphType)
        case {'a','auto'}
            axisMode = 'auto';
        case {'n','nought'}
            axisMode = 'nought';
        otherwise
            axisMode = 'manual';
            axisLimits = sscanf(graphType,'%f|%fx%f|%f');
    end
    
    Graphs(j) = struct('name',graphName,...
        'cols',cols,'axisMode',axisMode,'axisLimits',axisLimits);
end
end