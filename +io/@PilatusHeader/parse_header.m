function Header = parse_header(obj)

spaced_header_lines = obj.spaced_header_lines();

non_optional_keywords = fieldnames(obj.NON_OPTIONAL_KEYWORDS);
optional_keywords = fieldnames(obj.OPTIONAL_KEYWORDS);
all_keywords = [non_optional_keywords; optional_keywords];

% assign all non-optional keywords to defaults
for iKey = 1:length(non_optional_keywords)
    thisKeyword = non_optional_keywords{iKey};
    thisKey = obj.NON_OPTIONAL_KEYWORDS.(thisKeyword);
    Header.(thisKeyword) = obj.datatype_handling({''},...
                    thisKeyword, thisKey.type);
end


nKey = length(all_keywords);
nLines = length(spaced_header_lines);

for iKey = 1:nKey
    thisKeyword = all_keywords{iKey};
    if isfield(obj.OPTIONAL_KEYWORDS,thisKeyword)
        thisKey = obj.OPTIONAL_KEYWORDS.(thisKeyword);
    else
        thisKey = obj.NON_OPTIONAL_KEYWORDS.(thisKeyword);
    end
    for iLine = 1:nLines
        thisLine = spaced_header_lines{iLine};
        
        if regexp(thisLine,thisKey.pattern)
            try
                % add 1 to value_indices, because python
                % assumes indices start at zero
                if ischar(thisKey.value_indices) && ...
                        strcmp(thisKey.value_indices,'tokens')
                    parsedLine = regexp(thisLine,...
                        thisKey.pattern,'tokens');
                    values = strsplit(strtrim(parsedLine{1}{1}));
                else
                    splitLine = strsplit(strtrim(thisLine));
                    values = splitLine(thisKey.value_indices + 1);
                end
            catch
                error('Failed to parse header line:\n%s',...
                    thisLine);
            end
            if ~isempty(values)
                value = obj.datatype_handling(values,...
                    thisKeyword, thisKey.type);
                if ~isempty(value)
                    Header.(thisKeyword) = value;
                end
            end
        end
    end
end
end
