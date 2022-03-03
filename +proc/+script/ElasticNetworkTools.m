classdef ElasticNetworkTools < util.propertyValueConstructor
    %ElasticNetworkTools - generate elastic network models for lattice
    %dynamics and diffuse scattering simulation
    
    properties
        Atoms table
        tlsorigin (1,3) double = [0,0,0]
        cutoffdistance (1,1) double = 4
        Basis (1,1) latt.Basis = latt.Basis(100,100,100,90,90,90)
        UnitCellOperators (1,:) symm.SymmetryOperator = symm.SymmetryOperator();
        groupAttributes (1,:) cell = {'mdxGroupByResidue'};
    end
    properties(Dependent)
        com % compute the center of mass
        M % by default, return the
    end
    
    methods
        function obj = ElasticNetworkTools(varargin)
            %ElasticNetworkTools
            obj@util.propertyValueConstructor(varargin{:});
        end
        
        function val = get.com(obj)
            G = obj.table2group(obj.Atoms);
            val = G.com;
        end
        
        function M = calc_M(obj,mode1,mode2)
            if nargin < 2 || isempty(mode1)
                mode1 = 'grouped';
            end
            if nargin < 3 || isempty(mode2)
                mode2 = 'cell';
            end
            switch lower(mode1)
                case 'grouped'
                    G = obj.group_atoms();
                case 'ungrouped'
                    G = obj.table2group(obj.Atoms);
                otherwise
                    error('mode1 not recognized');
            end
            M = G.tlmass;
            switch lower(mode2)
                case 'cell'
                    ncellops = numel(obj.UnitCellOperators);
                    M = kron(speye(ncellops),M);
                case 'asu'
                    % do nothing
                otherwise
                    error('mode2 not recognized');
            end
        end
        
        function val = get.M(obj)
            % for compatibility with earlier scripts
            val = obj.calc_M('ungrouped','cell');
        end
        
        function CellOps = map_operators_to_cell(obj,Ops,asucenter)
            % add unit cell translations to the cell operators so that
            % the point asucenter is always mapped into the region [0,1]
            %
            if nargin < 3 || isempty(asucenter)
                asucenter = obj.tlsorigin;
            end
            B = obj.Basis.orthogonalizationMatrix;
            
            x0 = inv(B)*asucenter(:); % in fractional coordinates
            
            % deal with the possibility that asucenter is outside the
            % defined unit cell
            [~,ind0] = ismember(symm.SymmetryOperator(),Ops);
            x0 = Ops(ind0)*x0;
            cellshift0 = mod(x0,[1;1;1]) - x0;
            
            CellOps = Ops;
            
            %make sure the new com is within the unit cell (fractional coords [0,1])
            for j=1:length(Ops)
                xj = Ops(j)*x0;
                cellshift = mod(xj,[1;1;1]) - xj;
                CellOps(j).t = Ops(j).t + cellshift - cellshift0;
            end
            
        end
        
        function T = internalContactSearch(obj)
            
            nlist = proc.script.CoordinateTools.asu_neighbor_search(...
                obj.Atoms.x,obj.Atoms.y,obj.Atoms.z,...
                obj.cutoffdistance,obj.Atoms.mdxAltGroup);
            
            % filter internal neighbors by group, generate table
            group_assignment = obj.find_atom_groups();
            [a1,a2] = nlist2indices(nlist);
            
            % remove contacts between groups
            isSameGroup = group_assignment(a1)==group_assignment(a2);
            a1 = a1(~isSameGroup);
            a2 = a2(~isSameGroup);
            
            % remove duplicates
            isunique = a1 < a2;
            a1 = a1(isunique);
            a2 = a2(isunique);
            
            T = table(a1,a2);
            
        end
        
        function T = externalContactSearch(obj)
            
            B = obj.Basis.orthogonalizationMatrix;
            Ops = obj.UnitCellOperators;
            
            [nlist,Aout] = proc.script.CoordinateTools.symmetry_neighbor_search(...
                obj.Atoms.x,obj.Atoms.y,obj.Atoms.z,...
                B,Ops,obj.cutoffdistance);
            
            [o1,o2,c2,interface,a1,a2] = expandnlist2cell(nlist,Aout,Ops);
            
            T = table(o1,o2,c2,interface,a1,a2);
            
        end
        
        function [T,G,index,Tg] = externalModel(obj)
            
            T = obj.externalContactSearch();
            
            [T,G,index,Tg] = obj.contacts2nodes(T);
        end
        
        function [Tasu,G,index,Tg] = internalModel(obj)
            
            Tasu = obj.internalContactSearch();
            
            nr = size(Tasu,1);
            c2 = zeros(nr,3);
            interface = zeros(nr,1);
            o1 = ones(nr,1);
            o2 = ones(nr,1);
            
            Tasu = [table(o1,o2,c2,interface),Tasu];
            
            [Tasu,G,index,Tg] = obj.contacts2nodes(Tasu);
            
        end
        
        function [T,G,index,Tg] = fullModel(obj)
            
            Tasu = obj.internalContactSearch();
            
            nr = size(Tasu,1);
            c2 = zeros(nr,3);
            interface = zeros(nr,1);
            o1 = ones(nr,1);
            o2 = ones(nr,1);
            
            Tasu = [table(o1,o2,c2,interface),Tasu];
            
            T = repmat({Tasu},numel(obj.UnitCellOperators),1);
            for j=1:size(T,1)
                T{j}.o1(:) = j;
                T{j}.o2(:) = j;
            end
            T = cat(1,T{:});
            
            Tcell = obj.externalContactSearch();
            
            T = cat(1,T,Tcell);
            
            [T,G,index,Tg] = obj.contacts2nodes(T);
        end
        
        function ENM = exportModel(obj,T,G)
            ENM = nm.ElasticNetwork_v2(...
                'Cell', nm.Cell_v2(...
                'AsymmetricUnit',G,...
                'Basis',obj.Basis,...
                'UnitCellOperators',obj.UnitCellOperators),...
                'Edges',T);
        end
        
        function [ENMasu,ENMcell] = exportHybridModel(obj,Tasu,Gasu,Tcell,Gcell)
            ENMasu = nm.ElasticNetwork_v2(...
                'Cell', nm.Cell_v2(...
                'AsymmetricUnit',Gasu,...
                'Basis',obj.Basis),...
                'Edges',Tasu);
            ENMcell = obj.exportModel(Tcell,Gcell);
        end
        
        function [T,G,index,Tg] = contacts2nodes(obj,T)
            nodes = unique([T.a1;T.a2]);
            
            isNode = ismember((1:size(obj.Atoms,1))',nodes);
            
            [G,index,Tg] = obj.group_atoms(isNode);
            
            ind = cat(1,index{:});
            if ~isempty(ind) 
                nodemap = accumarray(ind,1:numel(ind),[size(obj.Atoms,1),1],@max);
            else % special case where T is empty
                nodemap = [];
            end
            
            T.a1 = nodemap(T.a1);
            T.a2 = nodemap(T.a2);
            
        end
        
        function [group_assignments,Tg] = find_atom_groups(obj)
            
            if isempty(obj.groupAttributes)
                group_assignments = ones(size(obj.Atoms,1),1);
                Tg = table(1);
                return;
            end
            
            A = obj.Atoms(:,obj.groupAttributes);
            [A,sortorder] = sortrows(A);
            [~,unsortorder] = sort(sortorder);
            [group_assignments,Tg] = findgroups(A);
            group_assignments = group_assignments(unsortorder);
        end
        
        function [G,index,Tg] = group_atoms(obj,isIncl)
            
            [grp,Tg] = obj.find_atom_groups();
            
            ngroups = size(Tg,1);
            
            A = obj.Atoms;
            ind0 = 1:size(A,1);
            
            if nargin > 1 && ~isempty(isIncl)
                grp = grp(isIncl);
                ind0 = ind0(isIncl);
            end
            
            if isempty(grp)
                G = nm.Group([],[],[],[]);
                index = {};
                return;
            end
            
            assert(all(~isnan(grp)))
            
            %G = splitapply(@obj.array2group,A,A.tmp_group(:))';
            index = accumarray(grp,ind0,[ngroups,1],@(v) {v});
            G = cellfun(@(g) obj.table2group(A(g,:)),index);
            
        end
        
        function G = table2group(obj,A)
            G = obj.array2group(A.x,A.y,A.z,A.mdxAtomicMass,A.occupancy);
        end
        
        function G = array2group(obj,x,y,z,mass,occupancy)
            if nargin == 4
                G = nm.Group(x,y,z);
            elseif nargin == 5
                G = nm.Group(x,y,z,mass);
            elseif nargin == 6
                G = nm.Group(x,y,z,mass.*occupancy);
            else
                error('incorrect number of arguments');
            end
            G.ori = obj.tlsorigin;
        end
        
    end
    
    methods(Static)
        function ENT = initialize(Atoms,Basis,SpaceGroup)
            Atoms = Atoms(~Atoms.isHet & Atoms.mdxAtomicSymbol~="H",:); % remove waters and hydrogens
            ENT = proc.script.ElasticNetworkTools('Atoms',Atoms,'Basis',Basis);
            ENT.tlsorigin = ENT.com;
            ENT.UnitCellOperators = ENT.map_operators_to_cell(SpaceGroup.generalPositions);
        end
    end
end




function [a1,a2] = nlist2indices(nlist)

nlist = cellfun(@(n) n(:),nlist,'Uni',0);

n = cellfun(@numel,nlist);
a1 = arrayfun(@(nn,jj) repmat(jj,nn,1),n,(1:numel(nlist))','Uni',0);

a1 = cat(1,a1{:});
a2 = cat(1,nlist{:});

end



function [o1,o2,c2,interface,a1,a2] = expandnlist2cell(nlist,Aout,Ops)

verbose = true; % for debugging, its good to leave this on.

if verbose, fprintf(1,'Expanding ASU contact model to the unit cell\n'); end

[AsuNeighbors,~,ic] = unique(Aout(:,{'isymop','xshift','yshift','zshift'}),'rows');
[Tnmap,CellNeighborOps] = nhoodmapping(Ops,AsuNeighbors);

[a1,r2] = nlist2indices(nlist);
op2 = ic(r2);
a2 = Aout.iatom(r2);

if verbose, fprintf(1,'  The ASU model has %d edges\n',numel(a1)); end

nOps = numel(Ops);

T = cell(nOps,1);
nmap = accumarray([Tnmap.icellop,Tnmap.iasunop],Tnmap.icellnop); % matrix form for fast lookup
for j=1:nOps
    isymop1 = repmat(j,size(a1));
    inop2 = nmap(j,op2)';
    asuindex = (1:numel(a1))';
    T{j} = [table(isymop1,a1,a2,asuindex),CellNeighborOps(inop2,:)];
    T{j}.interface = CellNeighborOps.interface(nmap(1,op2)); % <-- HACK ALERT!
end

%Tasu = T{1};

T = cat(1,T{:});

% find redundant entries
isInternal = T.xshift==0 & T.yshift==0 & T.zshift==0;
isRedundant = isInternal & (T.isymop < T.isymop1);

if verbose
    fprintf(1,'  After symmetry expansion there are %d edges ',size(T,1));
    fprintf(1,'    (%d internal + %d external)\n',nnz(isInternal),nnz(~isInternal));
    fprintf(1,'  Of the %d internal edges, %d are redundant\n',nnz(isInternal),nnz(isRedundant));
    fprintf(1,'  Deleting the redundant entries.\n');
end

T = T(~isRedundant,:);

% assign outputs
c2 = table2array(T(:,{'xshift','yshift','zshift'}));
o1 = T.isymop1;
o2 = T.isymop;
a1 = T.a1;
a2 = T.a2;
interface = T.interface;

if verbose, fprintf(1,'  Now there are %d edges.\nDone!\n\n',size(T,1)); end

end



function [T,CellNeighborOps] = nhoodmapping(CellOps,AsuNeighbors)

AsuNeighborOps = symmtable2array(AsuNeighbors,CellOps);

interface = assigninterface(AsuNeighborOps);

nCellOps = numel(CellOps);

T = cell(nCellOps,1);

for j=1:nCellOps
    o = arrayfun(@(op) CellOps(j)*op,AsuNeighborOps);
    T{j} = symmarray2table(o,CellOps);
    T{j}.icellop = repmat(j,size(T{j},1),1);
    T{j}.iasunop = (1:size(T{j},1))';
    T{j}.interface = interface;
end

T = cat(1,T{:});

[~,ind,T.icellnop] = unique(T(:,{'isymop','xshift','yshift','zshift'}),'rows','stable');

CellNeighborOps = T(ind,{'isymop','xshift','yshift','zshift','interface'});

T = T(:,{'icellop','iasunop','icellnop'});

end


function interface = assigninterface(o)

[ia,locb] = ismember(arrayfun(@inv,o),o);
assert(all(ia));

interface = (1:numel(o))';
isflipped = locb < interface; % choose these to flip
interface(isflipped) = interface(locb(isflipped)); % assign equivalent interface
[~,~,interface] = unique(interface,'stable'); % re-number them
interface(isflipped) = -interface(isflipped); % add minus sign if flipped

end


function A = symmtable2array(T,Ops)

A = rowfun(@(op,x,y,z) isymm2op(op,x,y,z,Ops),T,'OutputFormat','uniform');

end

function T = symmarray2table(A,Ops)

[isymop,xshift,yshift,zshift] = arrayfun(@(op) op2isymm(op,Ops),A);

T = table(isymop,xshift,yshift,zshift);

end


function Op = isymm2op(isymop,xshift,yshift,zshift,Ops)

Op = symm.SymmetryOperator(eye(3),[xshift,yshift,zshift])*Ops(isymop);

end

function [isymop,xshift,yshift,zshift] = op2isymm(Op,Ops)

denom = 48;
opfun = @(op) round([op.r(:)',op.t(:)']*denom);

Op = opfun(Op);
Ops = cell2mat(arrayfun(opfun,Ops(:),'Uni',0));

a = Op - Ops;
b = a;
b(:,10:12) = mod(b(:,10:12),denom);

isymop = find(all(b==0,2),1,'first');

xshift = a(isymop,10)/denom;
yshift = a(isymop,11)/denom;
zshift = a(isymop,12)/denom;

end