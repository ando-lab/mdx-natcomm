function [Atoms,Basis,SpaceGroup] = pdbImport(fileName,varargin)

if isempty(varargin)
    chainIDs = {};
else
    chainIDs = varargin;
end

[Cryst,Atoms,HetAtoms,AnisoTemps]  = io.pdb.read(fileName);

isHet = [false(size(Atoms,1),1);...
    true(size(HetAtoms,1),1)];

Atoms = [Atoms;HetAtoms]; 

if ~isempty(chainIDs)
    isHet = isHet | ~cellfun(@(ch) any(strcmp(ch,chainIDs)),Atoms.chainID);
end

%clear HetAtoms

[fatom,Elements,elindex] = latt.GaussianAtom.pdbImport(Atoms,AnisoTemps);

[~,~,alt] = unique(Atoms.altLoc); % convert alternate conformer to numeric
[~,~,chain] = unique(Atoms.chainID); % convert chain to numeric

occ = Atoms.occupancy;
x = Atoms.x;
y = Atoms.y;
z = Atoms.z;

masses = [Elements.atomicMass];
mass = masses(elindex)';

r = [Elements(elindex).vdwRadius]'; % vdw Radius

% assign resSeq of heteroatoms to that of nearest protein atom
resNum = Atoms.resSeq;
ix = find(~isHet);
A = Atoms(ix,:);
B = Atoms(isHet,:);
ind = knnsearch([A.x,A.y,A.z],[B.x,B.y,B.z]);
resNum(isHet) = resNum(ix(ind));
chain(isHet) = chain(ix(ind));
serial = Atoms.serial; 

% set atom output table
name = Atoms.name;
element = {Elements(elindex).atomicSymbol}';

Atoms = table(name,x,y,z,occ,chain,resNum,alt,isHet,element,mass,r,fatom,serial);

Basis = latt.OrientedBasis(Cryst.a,Cryst.b,Cryst.c,Cryst.alpha,Cryst.beta,Cryst.gamma);
SpaceGroup = symm.SpaceGroup(strrep(Cryst.sGroup,' ','')); % P1

end