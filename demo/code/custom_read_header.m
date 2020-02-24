function h = custom_read_header(fn)
info = imfinfo(fn);
header_text = info.ImageDescription;
parsed_header = regexp(header_text,'(?<name>\S+)\s*:\s*(?<value>\S*)','names');

h = struct();
for j=1:numel(parsed_header)
    h.(parsed_header(j).name) = str2double(parsed_header(j).value);
end

end