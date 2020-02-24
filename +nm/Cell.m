classdef Cell
    %CELL - model of unit cell contents
    
    properties
        AsymmetricUnit = nm.Group.empty(); % array of Group objects
        Basis = latt.Basis.empty(); % basis vectors of the unit cell
        SpaceGroup = symm.SpaceGroup.empty(); % space group
        mapToCell = true; % if true, makes sure that COM is mapped to unit cell. otherwise, UnitCellOperators returns general positions
    end
    properties (Dependent = true)
        UnitCellOperators
        B
    end
    
    methods
        
        function obj = Cell(B,SG,ASU)
            %CELL
            if nargin >= 1 && ~isempty(B)
                obj.Basis = B;
            end
            if nargin >= 2 && ~isempty(SG)
                obj.SpaceGroup = SG;
            end
            if nargin >= 3 && ~isempty(ASU)
                obj.AsymmetricUnit = ASU;
            end
        end
        
        function Bval = get.B(obj)
            if isa(obj.Basis,'latt.OrientedBasis')
                Om = obj.Basis.orientationMatrix*obj.Basis.orthogonalizationMatrix;
            else
                Om = obj.Basis.orthogonalizationMatrix;
            end
            Bval = symm.AffineTransformation(Om,[0;0;0]);
            
        end
        
        function P = tl2uxyz(obj)

                    Ops = obj.UnitCellOperators;
                    B = obj.B;
                    Binv = inv(B);
                    numAtoms = numel([obj.AsymmetricUnit.x]);
                    numOps = numel(Ops);
                    P0 = tl2uxyz(obj.AsymmetricUnit);
                    Pn = cell(1,numOps);
                    
                    for n = 1:numel(Ops)
                        Op = B*Ops(n)*Binv;
                        R = sparse(Op.r);
                        shiftVec = sparse(n,1,1,numOps,1);
                        shiftMat = kron(shiftVec,speye(numAtoms));
                        Pn{1,n} = kron(shiftMat,R)*P0;
                    end
                    P = cell2mat(Pn);
        end
        
        function M = tlmass(obj)
            M = kron(speye(numel(obj.UnitCellOperators)),obj.AsymmetricUnit.tlmass);
        end
        
        function Ops = get.UnitCellOperators(obj)
            
            Ops = obj.SpaceGroup.generalPositions; % initialize unit cell generators
                        
            if obj.mapToCell
            
            B = obj.B;
            
            r0com = com(obj.AsymmetricUnit);
            r0comfrac = inv(B)*r0com; % in fractional coordinates
            
            
            [~,ind0] = ismember(symm.SymmetryOperator(eye(3),[0;0;0]),Ops);
            r1comfrac = Ops(ind0)*r0comfrac;
            cellshift0 = mod(r1comfrac,[1;1;1]) - r1comfrac;
            
            % make sure the new com is within the unit cell (fractional coords [0,1])
            for j=1:length(Ops)
                r1comfrac = Ops(j)*r0comfrac;
                cellshift = mod(r1comfrac,[1;1;1]) - r1comfrac;
                Ops(j).t = Ops(j).t + cellshift - cellshift0;
            end
            
            end
        end
        
        function UC = UnitCell(obj)
            
            Ops = obj.UnitCellOperators;
            B = obj.B;
            
            numOps = numel(Ops);
            numGroups = numel(obj.AsymmetricUnit);
            
            allGroups = cell(1,numOps); % column vector
            
            for n=1:numOps
                allGroups{n} = obj.AsymmetricUnit(:)'; % row vector
                for g=1:numGroups
                    r = [allGroups{n}(g).x;allGroups{n}(g).y;allGroups{n}(g).z];
                    r1 = B*Ops(n)*inv(B)*r;
                    allGroups{n}(g).x = r1(1,:);
                    allGroups{n}(g).y = r1(2,:);
                    allGroups{n}(g).z = r1(3,:);
                end
            end
            
            UC = allGroups; % row is operator, column is group
            
        end
        
        function [opindex,cellShift,interface,d] = neighborSearch(obj,distanceCutoff)
            % neighborSearch - returns the unit cell operator indices and
            % fractional cell shift vectors for each ASU that neighbors the
            % current one. A neighbor is defined according to the center of
            % mass and the distanceCutoff: two ASUs with coms within d <=
            % distance Cutoff are considered to neighbors.
            
            Ops = obj.UnitCellOperators;
            
            % make the symmetry operators
            [n1,n2,n3] = ndgrid(-1:1,-1:1,-1:1); % 3x3x3 supercell
            n123 = [n1(:),n2(:),n3(:)];
            
            SuperCellOps = arrayfun(...
                @(a1,a2,a3) symm.SymmetryOperator(eye(3),[a1;a2;a3]),...
                n123(:,1),n123(:,2),n123(:,3));
            
            % compose super cell and unit cell operators
            AllOps = symm.SymmetryOperator.empty();
            for s=1:numel(SuperCellOps)
                for u=1:numel(Ops)
                    AllOps(s,u) = SuperCellOps(s)*Ops(u);
                end
            end
            
            % calculate distance between ASU center of mass and those of
            % neighbors
            B = obj.B;
            r0com = com(obj.AsymmetricUnit);
            r1com = arrayfun(@(op) B*op*inv(B)*r0com,AllOps,'UniformOutput',false);
            d = cellfun(@(v) sqrt(sum((v-r0com).^2,1)),r1com);
            
            [slist,ulist] = find(d <= distanceCutoff); % potential contacts
            d = d(sub2ind(size(d),slist,ulist));
            
            [d,dsortorder] = sort(d,'ascend'); % sort in ascending order
            
            % temporarily remove the first one, which has d~0 and is the
            % self-operator
            d = d(2:end);
            ulist = ulist(dsortorder(2:end));
            slist = slist(dsortorder(2:end));
            
            % now, label the remaining interfaces
            neighborOps = AllOps(sub2ind(size(AllOps),slist,ulist));
            [ia,locb] = ismember(arrayfun(@inv,neighborOps),neighborOps);
            
            
            %assert(all(ia)) % just in case... this should always be true
            if ~all(ia)
                
                % the above assertion sometimes fails if one of the neighbor
                % relations does not have its inverse within the confines of a
                % 3x3x3 supercell. In this case, it is fair to assume that
                % there is not a lattice contact with this neighbor, and it can
                % be eliminated.
                
                warning(['at least one symmetry neighbor does not have its '...
                    'inverse within the confines of a 3x3x3 supercell, ',...
                    'and will be removed']);
                
                locb = locb(ia);
                neighborOps = neighborOps(ia);
                slist = slist(ia);
                ulist = ulist(ia);
                d = d(ia);
            end
            
            
            
            interface = (1:numel(neighborOps))'; % initialize list of interfaces
            isflipped = locb < interface; % choose these to flip
            interface(isflipped) = interface(locb(isflipped)); % assign equivalent interface
            [~,~,interface] = unique(interface); % re-number them
            interface(isflipped) = -interface(isflipped); % add minus sign if flipped
            
            % assign outputs
            opindex = ulist;
            cellShift = n123(slist,:);
            
            % add the self operator back in at the beginning, labeled as
            % interface 0
            opindex = [1;opindex];
            cellShift = [0,0,0;cellShift];
            interface = [0;interface];
            d = [0;d];
        end
        
        function [i,j,d,vij] = latticeContactSearch(obj,contactDistance,opindex,cellShift)
            % latticeContactSearch - find all atoms within contactDistance of an
            % atom in the ASU that are related by a symmetry operator
            
            % get a copy of the unit cell operator
            Op = obj.UnitCellOperators(opindex);
            
            % make the symmetry operator
            CellOp = symm.SymmetryOperator(eye(3),cellShift');
            
            % operator for unit cell basis transformation
            B = obj.B;
            
            % copy all the coordinates of atoms in the ASU
            r0 = [[obj.AsymmetricUnit.x]; [obj.AsymmetricUnit.y]; [obj.AsymmetricUnit.z]];
            
            r1 = (B*CellOp*Op*inv(B))*r0;
            [i,j,d] = findContactsKD(r0,r1,contactDistance);
            vij = (r1(:,j)' - r0(:,i)').*repmat(1./d,1,3);
            
        end
        
        function [i,j,d,vij] = internalContactSearch(obj,contactDistance)
            % internalContactSearch - find atoms belonging to different groups
            % within the ASU that are within contactDistance of each other
            i = [];
            j = [];
            d = [];
            vij = [];
            nGroups = numel(obj.AsymmetricUnit);
            
            r = [obj.AsymmetricUnit.x; obj.AsymmetricUnit.y; obj.AsymmetricUnit.z];
            
            nAtoms = arrayfun(@(v) numel(v.x),obj.AsymmetricUnit);
            atomIndex = arrayfun(@(n,n0) (1:n) + n0,nAtoms,cumsum([0,nAtoms(1:(end-1))]),'UniformOutput',false);
            for ind1 = 1:(nGroups-1)
                index1 = atomIndex{ind1};
                index2 = [atomIndex{(ind1+1):nGroups}];
                %disp([numel(index1),numel(index2)]);
                [i1,j1,d1] = findContactsKD(r(:,index1),r(:,index2),contactDistance);
                index1 = index1(i1);
                index2 = index2(j1);
                
                v = (r(:,index2)' - r(:,index1)').*repmat(1./d1,1,3);
                
                i = [i; index1(:)];
                j = [j; index2(:)];
                d = [d; d1];
                vij = [vij; v];
            end
        end
        
        function T = contactSearch(obj,contactDistance)
            
            % 1) find neighbors
            
            distanceCutoff = 2*rmax(obj.AsymmetricUnit)+contactDistance;
            [opindex,cellShift,interface] = neighborSearch(obj,distanceCutoff);
            
            % 2) search for contacts over all the unique interfaces
            contactList = cell(numel(interface),1);
            for ind=1:numel(interface)
                if interface(ind) < 0
                    continue; % deal with this later
                end
                
                if interface(ind) == 0
                    [i,j,d] = internalContactSearch(obj,contactDistance);
                else % interface(ind) > 0
                    % neighbor interface
                    [i,j,d] = latticeContactSearch(obj,contactDistance,opindex(ind),cellShift(ind,:));
                end
                contactList{ind} = [i,j,d];
            end
            
            % 3) assign the non-unique interfaces (no search required)
            for ind=1:numel(interface)
                if interface(ind) >= 0
                    continue;
                end
                symind = find(-1*interface(ind) == interface,1,'first');
                
                % assign contacts, but swap i,j
                contactList{ind} = contactList{symind}(:,[2,1,3]);
                
            end
            
            % 4) get rid of neighbors that don't have contacts
            hasContact = ~cellfun(@isempty,contactList);
            contactList = contactList(hasContact);
            opindex = opindex(hasContact);
            cellShift = cellShift(hasContact,:);
            interface = interface(hasContact);
            
            % re-name interfaces
            isInternal = interface == 0;
            [~,~,irow] = unique(abs(interface(~isInternal)));
            interface(~isInternal) = sign(interface(~isInternal)).*irow;
            
            % 5) expand contacts to the unit cell
            T0 = expandUniqueInterfaces(obj,opindex,cellShift,interface);
            
            % 6) combine the results in an output table
            
            % o = operator, g = group, a = atom, c = cell
            % d = distance from 1 to 2
            T = table([],[],[],[],[],[],[],...
                'VariableNames',{'o1','o2','c2','interface','a1','a2','d'});
            
            for ind=1:size(T0,1) % loop over all interfaces
                t0 = T0(ind,:);
                [~,ia] = ismember(t0.interface,interface);
                a1 = contactList{ia}(:,1);
                a2 = contactList{ia}(:,2);
                d = contactList{ia}(:,3);
                t0 = repmat(t0,numel(a1),1);
                t0 = [t0, table(a1,a2,d)];
                T = [T; t0];
            end
            
        end
        
        function T = expandUniqueInterfaces(obj,opindex,cellShift,interfaceIn)
            % make a copy of the unit cell operators
            Ops = obj.UnitCellOperators;
            
            % generate nearest neighbor super cell operators
            [n1,n2,n3] = ndgrid(-1:1,-1:1,-1:1);
            n123 = [n1(:),n2(:),n3(:)];
            
            SuperCellOps = arrayfun(...
                @(a1,a2,a3) symm.SymmetryOperator(eye(3),[a1;a2;a3]),...
                n123(:,1),n123(:,2),n123(:,3));
            
            % compose super cell and unit cell operators
            AllOps = symm.SymmetryOperator.empty();
            for s=1:numel(SuperCellOps)
                for u=1:numel(Ops)
                    AllOps(s,u) = SuperCellOps(s)*Ops(u);
                end
            end
            
            [~,cellindex] = ismember(cellShift,n123,'rows');
            
            % function that converts an operator into a unique list of
            % integers for fast lookup (the factor of 48 should work
            % for all space groups)
            opfun = @(op) round([op.r(:)',op.t(:)']*48);
            
            % list of all operators in integer form for fast lookup
            A = cell2mat(arrayfun(opfun,AllOps(:),'UniformOutput',false));
            
            % initialize the output table
            T= table([],[],[],[],'VariableNames',{'o1','o2','c2','interface'});
            
            for ind=1:numel(Ops)
                thisop = Ops(ind);
                neighborOps = AllOps(sub2ind(size(AllOps),cellindex,opindex));
                newops = arrayfun(@(o) thisop*o,neighborOps);
                
                a = cell2mat(arrayfun(opfun,newops(:),'UniformOutput',false));
                [ia,locb] = ismember(a,A,'rows');
                
                % remove ones that are outside
                locb = locb(ia);
                newops = newops(ia);
                
                [cellindex2,opindex2] = ind2sub(size(AllOps),locb);
                o2 = opindex2;
                c2 = n123(cellindex2,:);
                o1 = ind*ones(size(o2));
                interface = interfaceIn(ia);
                T = [T; table(o1,o2,c2,interface)];
                
            end
            
            % get rid of duplicates
            isDuplicate = ismember(T.c2,[0,0,0],'rows') & T.o1 > T.o2;
            T = T(~isDuplicate,:);
        end
        
    end
end


function [i,j,d] = findContactsKD(X,Y,rcutoff)
% KD-tree-based search to find all pairs of points in X,Y that are within a
% distance rcutoff of each other. Requires stats toolbox to use rangesearch

% search using KD tree method
[idx, dist] = rangesearch(X',Y',rcutoff);

% put values into the ijList array
jList = find(~cellfun(@isempty,idx));
nNeighbors = sum(cellfun(@numel,idx));

i = zeros(nNeighbors,1);
j = zeros(nNeighbors,1);
d = zeros(nNeighbors,1);

e = 0;
for ind1= 1:length(jList)
    jIndex = jList(ind1);
    iList = idx{jIndex};
    dList = dist{jIndex};
    for ind2 = 1:length(iList)
        e = e + 1;
        i(e) = iList(ind2);
        j(e) = jIndex;
        d(e) = dList(ind2);
    end
end

end