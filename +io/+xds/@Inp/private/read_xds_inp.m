function output_struct = read_xds_inp(fn)
    output_struct = initialize_struct();
    
    if nargin==0
        %output_struct = output_struct(1);
        return;
    end
    
    output_cell = read_xds_inp_cell(fn);
    
    for j=1:size(output_cell,1)
        field_name = xds_inp_translate(output_cell{j,1});
        if isempty(output_struct) || isempty(output_struct.(field_name))
            output_struct(1).(field_name) = output_cell{j,2};
        else
            output_struct(1).(field_name) = [output_struct.(field_name);...
                output_cell{j,2}];
        end
    end
end

function struct_init = initialize_struct()
    xds_commands = xds_inp_translate();
    fn = fieldnames(xds_commands);
    c = cell(2,length(fn));
    for j=1:length(fn)
        c{1,j} = fn{j};
        c{2,j} = {};
    end
    struct_init = struct(c{:});
end

function output_cell = read_xds_inp_cell(fn)

fid = fopen(fn,'r');

nocomment_expr = '^\s*(?<expr>[^\!]*)'; % strip comments
varnames_expr = '([A-Za-z_\-\(\)\/\.\'']+)(?=\=)'; % pull command strings
numeric_arg = '(?<=\=\s*)([0-9\.\-]+[0-9\.\s\-]*)$'; % pull numeric arguments
split_arg = '(?<=((^\=)|(\s+)))([a-zA-Z0-9_\-\/\.\?\~]+)'; % pull text arguments

output_cell = cell(0,2);

k = 0;

while ~feof(fid)
    l = fgetl(fid);
    l_nocomment = regexp(l,nocomment_expr,'names');
    if ~isempty(l_nocomment)
        [names,values] = regexp(l_nocomment.expr,varnames_expr,'match','split');
        for j=1:length(names)
            k = k+1;
            output_cell{k,1} = names{j};
            %fprintf(1,[names{j} '\n\t' values{j+1} '\n']);
            numval = regexp(values{j+1},numeric_arg,'match');
            if ~isempty(numval)
                output_cell{k,2} = sscanf(numval{1},'%f')';
            else
                cellstr_out = regexp(values{j+1},split_arg,'match');
                if length(cellstr_out) == 1
                    output_cell{k,2} = cellstr_out{1};
                else
                    output_cell{k,2} = cellstr_out;
                end
            end
        end
    end
end
fclose(fid);
end