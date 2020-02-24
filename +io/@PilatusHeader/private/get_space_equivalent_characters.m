function out = space_equivalent_characters()
% the following characters are ignored (treated as spaces) when parsing the
% header comments in a pilatus CBF image file
    out = '#:=,()';
end