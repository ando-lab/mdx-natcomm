function value = datatype_handling(obj,values,keyword,datatype)

% handle rare case of value "not set"
if strcmp(datatype,'float') && strcmp(values{1},'not')
    % NON_OPTIONAL_KEYWORDS should have value, at least NaN
    if isfield(obj.NON_OPTIONAL_KEYWORDS,keyword)
        value = NaN;
    else
        value = [];
    end
end

% do the conversion for standard cases
switch datatype
    case {'float','int'}
        value = str2double(values);
    case 'str'
        value = strjoin(values,' ');
    otherwise
        error('did not recognize datatype %s',datatype);
end
end