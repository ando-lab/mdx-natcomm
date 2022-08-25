function AA = sequence2AminoAcids(seq)
% sequence2AminoAcids - search for amino acid parameters in aa-variants-v1.cif
%
% returns an array of model.chem.AminoAcid objects.
%
% for example:
%   seq = {'Ala','Val','Tyr','Ile','Ala','Thr'};
%   AA = model.chem.sequence2AminoAcids(seq);
[mmcifData,ind] = io.mmcif.searchFile('+model/+chem/private/aa-variants-v1.cif',seq);
mmcifParsed = cellfun(@io.mmcif.parseRecord,{mmcifData.record},'UniformOutput',false);
M = cellfun(@io.mmcif.convert2Molecule,mmcifParsed,'UniformOutput',false);
M = [M{:}];
M = M(ind);
AA = model.chem.AminoAcid.empty();
for j=1:length(M)
    AA(j) = model.chem.AminoAcid(M(j));
end
end