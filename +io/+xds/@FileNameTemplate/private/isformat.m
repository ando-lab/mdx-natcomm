function tf = isformat(fmtstr)

matchStr = '\%(?<zero>[0]*)(?<nDigits>\d+)d';

[fileTemplate,matchedTokens] = regexp(fmtstr,matchStr,...
    'split','names');

nMatches = length(matchedTokens);

tf = nMatches>0;

end