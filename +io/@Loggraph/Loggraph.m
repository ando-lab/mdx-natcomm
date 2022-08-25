classdef Loggraph
    properties (SetAccess = protected)
        %fileName
        %LogTables
        %otherText
    end
    methods
        function obj = Loggraph()
            %obj.fileName = fileName;
            %obj.LogTables = obj.read(fileName);
        end
%         function plot(LogTables,graphNumber)
%             plot_loggraph(LogTable,graphNumber);
%         end
    end
    methods(Static)
        function [LogTables,otherText] = read(fileName)
            [LogTables,otherText] = read_loggraph(fileName);
        end
        function print(LogTable,graphNumber)
            print_loggraph(LogTable,graphNumber);
        end
        function LogTable = findTable(LogTables,titleSearchString)
            % titleSearchString can contain wildcard characters (*)
            isMatching = struct_search(LogTables,'title',titleSearchString);
            ixMatch = find(isMatching,1,'first');
            LogTable= LogTables(ixMatch);
        end
    end
end