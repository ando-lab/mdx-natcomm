function m = hkl2mult(h,k,l,Crystal)
% multiplicity of a given reflection
lg = Crystal.SpaceGroup.PointGroup.LaueGroup;
Ops = lg.generalPositionsMat;
m = zeros(size(h));
for j=1:numel(h)
    hklSym = cell2mat(cellfun(@(v) ([h(j),k(j),l(j)]*v)',Ops,...
        'UniformOutput',false))';
    m(j) = size(unique(hklSym,'rows'),1);
end
end