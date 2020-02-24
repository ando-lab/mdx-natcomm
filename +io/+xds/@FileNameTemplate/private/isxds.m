function tf = isxds(fmtstr)

matchStr = '(?<questionMarks>\?+)';

[fileTemplate,matchedTokens] = regexp(fmtstr,matchStr,...
    'split','names');

nMatches = length(matchedTokens);

tf = nMatches>0;

end