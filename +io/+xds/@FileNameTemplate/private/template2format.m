function fmtstr = template2format(fileTemplate,nDigits,isZeroPadded)
nDelims = length(nDigits);
delims = cell([1,nDelims]);

for j=1:nDelims
    % %04d
    if isZeroPadded(j)
        delims{j} = ['%0',sprintf('%d',nDigits(j)),'d'];
    else % not isZeroPadded(j)
        delims{j} = ['%',sprintf('%d',nDigits(j)),'d'];
    end
end

fmtstr = strjoin(fileTemplate,delims);

end