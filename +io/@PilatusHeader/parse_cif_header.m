function Header = parse_cif_header(obj)

keywords = fieldnames(obj.CIF_BINARY_KEYWORDS);

nKey = length(keywords);
nLines = length(obj.cif_lines);

for iKey = 1:nKey
    thisKeyword = keywords{iKey};
    thisKey = obj.CIF_BINARY_KEYWORDS.(thisKeyword);
    for iLine = 1:nLines
        thisLine = obj.cif_lines{iLine};
        thisParse = regexp(thisLine,thisKey.pattern,'tokens');
        if ~isempty(thisParse)
            switch thisKey.type
                case 'str'
                    value = thisParse{1}{1};
                case 'int'
                    value = str2double(thisParse{1}{1});
            end
            if ~isempty(value)
                Header.(thisKeyword) = value;
            end
        end
    end
end
end