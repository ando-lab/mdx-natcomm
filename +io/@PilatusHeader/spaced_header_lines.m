function value = spaced_header_lines(obj)
%spaced_header_lines Return header_lines with all space equivalent
%charecters converted to space.
value = cell(size(obj.header_lines));
nLines = length(obj.header_lines);

space_equivalent_expression = ['[',regexptranslate('escape',...
    obj.SPACE_EQUIVALENT_CHARACTERS),']'];

for iLine = 1:nLines
    thisLine = obj.header_lines{iLine};
    value{iLine} = regexprep(thisLine, ...
        space_equivalent_expression,' ');
end
end