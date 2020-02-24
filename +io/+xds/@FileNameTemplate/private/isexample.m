function tf = isexample(fmtstr)

matchStr = '(?<finalDigits>\d+)(?=\.)';

[fileTemplate,matchedTokens] = regexp(fmtstr,matchStr,...
    'split','names');

nMatches = length(matchedTokens);

tf = nMatches>0;

end