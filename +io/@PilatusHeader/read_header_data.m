function value = read_header_data(obj)
fid = fopen(obj.filepath,'r');

if fid==-1
    % did not read file. throw error
    error('no file found at %s\n',obj.filepath);
end

value = fread(fid,[1,obj.non_binary_length],'*char');
fclose(fid);
%header_end_location = regexp(value,'\x{0C}\x{1A}\x{04}\x{d5}');
%value = value(1:(header_end_location-1));
end