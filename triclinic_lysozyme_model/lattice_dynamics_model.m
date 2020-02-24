%% Lattice Dynamics Model for Triclinic Lysozyme

%% get c-alpha coordinates

load model/atomic_model.mat Atoms SpaceGroup Basis

% select the c-alpha atoms for the network model
isCA = strcmpi(Atoms.name,'CA') & strcmpi(Atoms.element,'C');
A = Atoms(isCA,:);

% do a weighted average of alpha carbon positions for alternate conformers
numRes = max(A.resNum);% 129;
ind = accumarray(A.resNum,(1:size(A,1)),[numRes,1],@(v) {v});
xyz = cell2mat(cellfun(@(j) nm.Group(A.x(j),A.y(j),A.z(j),A.occ(j)).com,ind(:)','UniformOutput',false));

% create a unit cell
G = nm.Group(xyz(1,:),xyz(2,:),xyz(3,:),1);
UC = nm.Cell(Basis,SpaceGroup,G); % unit cell
clear xyz ind numRes is isCA G

% create a second unit cell with all of the atoms (for calculating the mass
% tensor)
G = nm.Group(Atoms.x,Atoms.y,Atoms.z,Atoms.mass.*Atoms.occ);

UCfull = nm.Cell(Basis,SpaceGroup,G);

clear Basis SpaceGroup G A numRes Atoms

%% calculate mass tensor
M = UCfull.tlmass;

%% find c-alpha contacts based on elastic network model
enm = load('model/elastic_network_model.mat');

T0 = enm.N.Edges;

asu = enm.N.Cell.AsymmetricUnit;
ix = cell2mat(arrayfun(@(g,n) n*ones(1,numel(g.x)),asu,1:numel(asu),'UniformOutput',false));
clear asu
a1 = ix(enm.N.Edges.a1);
a2 = ix(enm.N.Edges.a2);
clear ix
T0.a1 = a1(:);
T0.a2 = a2(:); % with reference to residue number rather than atom number

[~,ix] = unique(T0(:,{'c2','a1','a2'}),'rows');
T0 = T0(ix,:);
T0 = T0(T0.interface ~= 0,:);

N = nm.Network(T0,UC);

%% save result
save model/lattice_dynamics_model.mat M N