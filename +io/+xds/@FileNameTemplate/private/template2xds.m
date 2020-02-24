function xdsstr = template2xds(fileTemplate,nDigits)
nDelims = length(nDigits);
delims = cell([1,nDelims]);

for j=1:nDelims
    delims{j} = repmat('?',1,nDigits(j));
end

xdsstr = strjoin(fileTemplate,delims);

end