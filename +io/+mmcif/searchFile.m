function [m,ind] = searchFile(fileName,sequence)
if ischar(sequence)
    sequence = {sequence};
end

sequence = upper(sequence);

aaToMatch = unique(sequence);
searchExpr = ['data_(?<name>' strjoin(aaToMatch,'|') ')\s(?<record>.*?)(?=($|data_))'];
m = regexpi(fileread(fileName),searchExpr,'names');

[ia,ind] = ismember(sequence,upper({m.name}));
if ~all(ia)
    warning('%d records were not found in %s:',sum(~ia),fileName);
    disp(sequence(~ia));
end

end