classdef Log < io.Loggraph
    properties (SetAccess = protected)
        otherText
        LogTables
        summaryText
    end
    methods
        function obj = Log(fileName)
            if nargin==0 || isempty(fileName)
                fileName = 'aimless.log';
            elseif isdir(fileName)
                fileName = fullfile(fileName,'aimless.log');
            end
            % logEntries = read_aimless_log(fileName);
            [obj.LogTables,obj.otherText] = read_loggraph(fileName);
            
            theResult = regexp([obj.otherText{:}],...
                '<!\-\-SUMMARY_BEGIN\-\->(.*)<!\-\-SUMMARY_END\-\->',...
                'tokens');
            obj.summaryText = theResult{1}{1};
            
        end
        
        function LogTable = findTable(obj,titleSearchString)
            % titleSearchString can contain wildcard characters (*)
            isMatching = struct_search(obj.LogTables,'title',titleSearchString);
            ixMatch = find(isMatching,1,'first');
            LogTable= obj.LogTables(ixMatch);
        end
        
        function matchingLines = findInSummary(obj,searchString,numLines)
            % searchstring is used to match the line(s) you are looking
            % for.
            if nargin==2
                numLines = 1;
            end
            
            summaryLines = strsplit(obj.summaryText,'\n');
            regExpStr = regexptranslate('wildcard',searchString);
            theResult = regexp(summaryLines,regExpStr,'start');
            % get rid of lines that do not match
            isMatching = cellfun(@(c) ~isempty(c),theResult); 
            startLineIndex = find(isMatching,1,'first');
            lineIndices =  startLineIndex + (1:numLines) - 1;
            lineIndices = lineIndices(lineIndices<=length(summaryLines));
            
            matchingLines = strjoin(summaryLines(lineIndices),'\n');
        end
    end
    methods(Static)
        
        function plot(LogTable,graphID)
            if isnumeric(graphID)
                graphNumber = graphID;
            else % graphID = graph name search string with wildcards
                isMatching = struct_search(LogTable.Graphs,...
                    'name',graphID);
                graphNumber = find(isMatching,1,'first');
            end
            plot_loggraph(LogTable,graphNumber);
        end
        
        function dataColumns = getData(LogTable,columnNames)
            dataColumns = [];
            for j=1:length(columnNames)
                isMatching = strcmp(LogTable.columns,columnNames{j});
                ixMatch = find(isMatching,1,'first');
                if isempty(ixMatch) % try case insensitive search
                    isMatching = strcmpi(LogTable.columns,columnNames{j});
                    ixMatch = find(isMatching,1,'first');
                    if isempty(ixMatch)
                        error('column %s not found',columnNames{j});
                    end
                end
                dataColumns = [dataColumns, LogTable.data(:,ixMatch)];
            end
        end
    end
end