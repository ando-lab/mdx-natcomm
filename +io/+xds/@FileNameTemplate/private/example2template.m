function [fileTemplate,nDigits,isZeroPadded] = example2template(fmtstr)

matchStr = '(?<finalDigits>\d+)(?=\.)';

[fileTemplate,matchedTokens] = regexp(fmtstr,matchStr,...
    'split','names');

nDigits = length(matchedTokens.finalDigits);
isZeroPadded = str2double(matchedTokens.finalDigits(1))==0;

end