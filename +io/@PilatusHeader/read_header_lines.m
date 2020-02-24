function value = read_header_lines(obj)
contents_pattern = '_array_data.header_contents\s+;(.*?);';
contents_match = regexp(obj.rawheader,contents_pattern,'tokens');
assert(~isempty(contents_match),'did not find header_contents');
value = regexp(contents_match{1}{1},'\r\n','split');
end