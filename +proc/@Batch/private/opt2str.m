function str = opt2str(opt)
str = '';
if isempty(opt)
    return;
end
fn = fieldnames(opt);
for j=1:length(fn)
    val = opt.(fn{j});
    if ischar(val)
        % do nothing
    elseif isnumeric(val) || isa(val,'logical')
        if numel(val)==1
            val = num2str(val);
        else
            val = mat2str(val);
        end
    elseif iscellstr(val) % assume cell string
        val = [strjoin(val,', ')];
    elseif iscell(val)
        for k=1:length(val)
            if isnumeric(val{k}) || isa(val{k},'logical')
                if numel(val{k})==1
                    val{k} = num2str(val{k});
                else
                    val{k} = mat2str(val{k});
                end
            end
        end
        val = [strjoin(val,', ')];
    end
    str = [str, sprintf('%s: %s\n',fn{j},val)];
end
end