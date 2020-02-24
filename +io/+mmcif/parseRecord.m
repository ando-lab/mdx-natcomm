function recordData = parseRecord(recordStr)

sectionSplit = '((^|\n+)#\s*($|\n+))';
kvLine = '_(?<name>\w+)\.(?<key>\w+)[ ]+(?<value>".+?"|[^ ]+)[ ]*';
loopStart = '^[\n ]*loop_[ ]*\n';
loopLine = '_(?<name>\w+)\.(?<key>\w+)[ ]*\n';
escapeStr = '("(.+?)"|\n;(.+);|(\S+))';

a = regexp(recordStr,sectionSplit,'split');

recordData = [];
remainder = '';

for j=1:length(a)
    if isempty(a{j})
        % do nothing
    elseif strncmp(a{j},'loop_',5) % read as table

        a1 = regexp(a{j},loopStart,'split');
        [a2,a3] = regexp([a1{:}],loopLine,'names','split');
        [a5,a6] = regexp([a3{:}],escapeStr,'match','split');
        cols = reshape(a5,length(a2),[])';
        thisTable = cell2table(cols,'variableNames',{a2.key});
        recordData.(a2(1).name) = thisTable;
        remainder = deblank([a6{:}]);

    elseif strncmp(a{j},'_',1) % read as struct

        [a2,a3] = regexp(a{j},kvLine,'names','split');
        recordData.(a2(1).name) = table2struct(cell2table({a2.value},'variableNames',{a2.key}));

        remainder = deblank([a3{:}]);

    else
        warning('did not recognize record type:');
        disp(a{j});
    end

    if ~isempty(remainder)
        warning('unread text remaining:');
        disp(remainder);
        remainder = '';
    end
end
end
