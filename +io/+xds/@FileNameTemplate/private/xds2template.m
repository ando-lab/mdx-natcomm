function [fileTemplate,nDigits,isZeroPadded] = xds2template(fmtstr)

matchStr = '(?<questionMarks>\?+)';

[fileTemplate,matchedTokens] = regexp(fmtstr,matchStr,...
    'split','names');

nMatches = length(matchedTokens);

nDigits = zeros([1,nMatches]);
for j = 1:nMatches
    nDigits(j) = length(matchedTokens(j).questionMarks);
end

% XDS file templates always assume zero-padded digits
isZeroPadded = true([1,nMatches]);

end