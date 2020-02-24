function value = read_cif_lines(obj)
contents_pattern = ['\n--CIF-BINARY-FORMAT-SECTION--',...
    '(.*?)\x{0C}\x{1A}\x{04}\x{d5}'];
contents_match = regexp(obj.rawheader,contents_pattern,'tokens');
assert(~isempty(contents_match),'did not find header_contents');
value = regexp(contents_match{1}{1},'\r\n','split');
end