function parm = templated_read(fid,template)

% initialize parm to have same fields as template
parm = template;
fn = fieldnames(template);
for j=1:length(fn)
    parm.(fn{j}) = [];
end

% read _all_ fields of template in order
for j=1:length(fn)
    thisFn = fn{j};
    value = [];
    while ~feof(fid) && isempty(value)
        value = sscanf(fgetl(fid),template.(thisFn));
    end
    if feof(fid)
        warning('end of file reached before finishing templated read')
    end
    if isnumeric(value)
        value = value';
    end
    parm.(thisFn) = value;
end
end