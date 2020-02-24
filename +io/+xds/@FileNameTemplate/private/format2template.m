function [fileTemplate,nDigits,isZeroPadded] = format2template(fmtstr)

matchStr = '\%(?<zero>[0]*)(?<nDigits>\d+)d';

[fileTemplate,matchedTokens] = regexp(fmtstr,matchStr,...
    'split','names');

nDigits = str2double({matchedTokens.nDigits});
isZeroPadded = str2double({matchedTokens.zero})==0;

%xdsDigitStr = repmat('?',1,nDigits);
%xdsTemplate = [splitTemplate{1},xdsDigitStr,splitTemplate{2}];
end