%% Internal Vibrations Model for Triclinic Lysozyme

%% load coordinates

load model/atomic_model.mat Atoms SpaceGroup Basis

% remove hetero atoms
Atoms = Atoms(~Atoms.isHet,:);

% for network model, remove hydrogens and average alt conformers
A = Atoms(~ismember(upper(Atoms.element),{'H'}),:);

% average alt conformers.
[~,ib,ix] = unique(A(:,{'resNum','name'}),'rows','stable');
ind = accumarray(ix,(1:size(A,1)),[numel(ib),1],@(v) {v});
xyz = cell2mat(cellfun(...
    @(ind1) nm.Group(A.x(ind1),A.y(ind1),A.z(ind1),A.occ(ind1)).com',...
    ind,'UniformOutput',false));
A = A(ib,{'name','x','y','z','resNum','element','mass'});
A.x = xyz(:,1);
A.y = xyz(:,2);
A.z = xyz(:,3);

clear ind xyz ib ix

%% assign atoms to groups

% sort by resSeq and make an array of groups (one for each residue)
[~,ib,ix] = unique(A(:,'resNum'),'rows','stable');
ind = accumarray(ix,(1:size(A,1)),[numel(ib),1],@(v) {v});
G = cellfun(@(ix) nm.Group(A.x(ix),A.y(ix),A.z(ix),1),ind');

% full group including occupancy*mass information
[~,ib,ix] = unique(Atoms(:,'resNum'),'rows','stable');
ind =  accumarray(ix,(1:size(Atoms,1)),[numel(ib),1],@(v) {v});
Gfull = cellfun(@(ix) nm.Group(Atoms.x(ix),Atoms.y(ix),Atoms.z(ix),Atoms.occ(ix).*Atoms.mass(ix)),ind');

%% create unit cell model
UC = nm.Cell(Basis,SpaceGroup,G); % unit cell

UCfull = nm.Cell(Basis,SpaceGroup,Gfull);

%% find bonds

contactDistance = 4; 
T = UC.contactSearch(contactDistance);

N = nm.Network(T,UC);

M = UCfull.tlmass;

%% save
save model/elastic_network_model.mat M N