function value = seek_binary_start(obj)
nPatternBytes = 4;
%pattern = '\x{0C}\x{1A}\x{04}\x{d5}';
pattern = '\x{0C}\x{1A}\x{04}'; % last byte sometimes missing?
pattern_match = regexp(obj.rawheader,pattern);
value = pattern_match + nPatternBytes - 1;
end