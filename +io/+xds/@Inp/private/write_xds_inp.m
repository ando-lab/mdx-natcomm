function write_xds_inp(input_struct,fn,header_comment_lines)

    if nargin<3
        header_comment_lines = {['XDS input file created by write_xds_inp.m on ' date]};
    end

    field_names = fieldnames(input_struct);
    nfields = length(field_names);
    
    input_cell = cell(nfields,2);
    xds_commands = xds_inp_translate();
    for j=1:nfields
        input_cell{j,1} = xds_commands.(field_names{j});
        input_cell{j,2} = input_struct.(field_names{j});
    end
    
    if nargin==1
        write_xds_inp_cell(input_cell,1);
    else
        fid = fopen(fn,'w');
        for j=1:length(header_comment_lines)
            fprintf(fid,'!%s\n',header_comment_lines{j});
        end
        write_xds_inp_cell(input_cell,fid);
    	fclose(fid);
    end
end

function write_xds_inp_cell(input_cell,fid)

for j=1:length(input_cell)
    cmd = input_cell{j,1};
    args = input_cell{j,2};
    
    if isempty(args) && isnumeric(args)
        continue;
    elseif isnumeric(args)
        for k=1:size(args,1)
            fprintf(fid,'%s=',cmd);
            for l=1:size(args,2)
                fprintf(fid,' %s',num2str(args(k,l)));
            end
            fprintf(fid,'\n');
        end
    elseif ischar(args)
        fprintf(fid,'%s= %s\n',cmd,args);
    elseif isempty(args) % empty cell array <-- PRINT!
        fprintf(fid,'%s= !\n',cmd);
    else % must be cell array of strings
        for k=1:size(args,1)
            fprintf(fid,'%s=',cmd);
            for l=1:size(args,2)
                fprintf(fid,' %s',args{k,l});
            end
            fprintf(fid,'\n');
        end
    end
end

end