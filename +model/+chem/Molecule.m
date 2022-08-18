classdef Molecule
    %MOLECULE represent chemical bonding in a molecule
    
    properties
        name = '' % long name for display purposes
        id = ''   % short name to serve as unique identifier
    end
    properties(Dependent = true)
        nAtoms % number of atoms
        nBonds % number of bonds
        nElectrons % total number of electrons
        formulaWeight
        totalCharge
        atomInventory
    end
    properties(SetAccess = protected)
        G % Graph object storing atom/bond information as nodes/edges
    end
    properties(Dependent = true, Hidden) % wrappers for accessing common elements of G
        x
        y
        z
        charge
        isDouble
        isCovalent
        isAromatic
    end
    
    methods
        function obj = Molecule(AtomTable,BondTable,varargin)
            %MOLECULE
            %
            % AtomTable contains information about the atoms in the
            % molecule. The first column is the name of the atom, and
            % should be a unique identifier. The other columns are
            % optional, however for member functions to work, the following
            % should be defined:
            %
            % AtomTable.x
            % AtomTable.y
            % AtomTable.z
            % AtomTable.type_symbol
            % AtomTable.charge
            %
            % BondTable contains information about bonding within the
            % molecule. The first two columns are either the indices of the
            % (undirected) bond between atoms, or the unique identifiers
            % (first column of AtomTable). The other columns are optional,
            % but the following should be defined for all member functions
            % to work properly:
            %
            % BondTable.isDouble % true/false for double/single bonds
            %
            % The following are also used by plot():
            %
            % BondTable.isCovalent % true/false for covalent/non-covanent
            
            if nargin==0
                return;
            end
            
            for j=2:2:length(varargin)
                obj.(varargin{j-1}) = varargin{j};
            end
            
            EdgeTable = table(table2cell(BondTable(:,1:2)),'VariableNames',{'EndNodes'});
            EdgeTable = [EdgeTable,BondTable(:,3:end)];
            NodeTable = AtomTable;
            NodeTable.Properties.VariableNames{1} = 'Name';
            
            if ~iscell(NodeTable.Name)
                NodeTable.Name = {NodeTable.Name}; % hack to fix tables with only 1 row
                NodeTable.type_symbol = {NodeTable.type_symbol}; % hack to fix tables with only 1 row
            end
            
            obj.G = graph(EdgeTable,NodeTable);
            
        end
        
        function val = get.nAtoms(obj)
            val = obj.G.numnodes;
        end
        function val = get.nBonds(obj)
            val = obj.G.numedges;
        end
        function obj = setCharge(obj,id,val)
            obj.G.Nodes.charge(obj.G.findnode(id)) = val;
        end
        function [AtomTable,atomIndex] = findAtom(obj,atomID)
            atomIndex = obj.G.findnode(atomID);
            AtomTable = obj.G.Nodes(atomIndex,:);
        end
        
        function M = removeAtom(obj,atomID)
            M = obj;
            M.G = obj.G.rmnode(atomID);
        end
        
        function [BondTable,bondIndex] = findBond(obj,atomID,atomID2)
            if nargin<3 || isempty(atomID2)
                atomID2 = obj.G.neighbors(atomID);
            end
            bondIndex = obj.G.findedge(atomID,atomID2);
            BondTable = obj.G.Edges(bondIndex,:);
        end
        
        function M = removeBond(obj,atomID,atomID2)
            M = obj;
            M.G = obj.G.rmedge(atomID,atomID2);
        end
        
        function Mnew = addAtom(obj,AtomTable)
            AtomTable.Properties.VariableNames{1} = 'Name';
            Mnew = obj; % copy of current molecule
            Mnew.G = obj.G.addnode(AtomTable);
        end
        
        function Mnew = addBond(obj,BondTable)
            EdgeTable = table(table2cell(BondTable(:,1:2)),'VariableNames',{'EndNodes'});
            EdgeTable = [EdgeTable,BondTable(:,3:end)];
            Mnew = obj; % copy of current molecule
            Mnew.G = obj.G.addedge(EdgeTable);
        end
        
        function val = get.x(obj)
            val = obj.G.Nodes.x;
        end
        function val = get.y(obj)
            val = obj.G.Nodes.y;
        end
        function val = get.z(obj)
            val = obj.G.Nodes.z;
        end
        function val = get.isDouble(obj)
            val = obj.G.Edges.isDouble;
        end
        function val = get.isCovalent(obj)
            val = obj.G.Edges.isCovalent;
        end
        function val = get.isAromatic(obj)
            val = obj.G.Edges.isAromatic;
        end
        function val = get.charge(obj)
            val = obj.G.Nodes.charge;
        end
        
        function Mnew = addMolecule(obj,obj2)
            Mnew = obj;
            if isempty(Mnew.G)
                Mnew.G = graph(obj2.G.Edges,obj2.G.Nodes);
            else
                Mnew.G = Mnew.G.addnode(obj2.G.Nodes);
                Mnew.G = Mnew.G.addedge(obj2.G.Edges);
            end
        end
        
        function Mnew = prependToName(obj,str)
            Mnew = obj;
            Mnew.G.Nodes.Name = ...
                cellfun(@(v) [str v],Mnew.G.Nodes.Name,...
                'UniformOutput',false);
        end
        function Mnew = appendToName(obj,str)
            Mnew = obj;
            Mnew.G.Nodes.Name = ...
                cellfun(@(v) [v str],Mnew.G.Nodes.Name,...
                'UniformOutput',false);
        end
        
        function Mnew = rotate(obj,Rmat)
            coords = table2array(obj.G.Nodes(:,{'x','y','z'}))';
            coords = Rmat*coords;
            Mnew = obj;
            Mnew.G.Nodes.x = coords(1,:)';
            Mnew.G.Nodes.y = coords(2,:)';
            Mnew.G.Nodes.z = coords(3,:)';
        end
        
        function Mnew = translate(obj,vec)
            Mnew = obj;
            Mnew.G.Nodes.x = obj.G.Nodes.x + vec(1);
            Mnew.G.Nodes.y = obj.G.Nodes.y + vec(2);
            Mnew.G.Nodes.z = obj.G.Nodes.z + vec(3);
        end
        
        function val = get.atomInventory(obj)
            [Type,~,ib] = unique(obj.G.Nodes.type_symbol);
            Number = accumarray(ib,1);
            val = table(Type,Number);
        end
        
        function val = get.formulaWeight(obj)
            atoms = obj.atomInventory;
            val = 0;
            for j=1:size(atoms,1)
                E = model.atom.Element(atoms.Type{j});
                val = val + E.atomicMass*atoms.Number(j);
            end
        end
        
        function val = get.totalCharge(obj)
            val = sum(obj.charge);
        end
        
        function val = get.nElectrons(obj)
            atoms = obj.atomInventory;
            val = 0;
            for j=1:size(atoms,1)
                E = model.atom.Element(atoms.Type{j});
                val = val + E.atomicNumber*atoms.Number(j);
            end
            val = val - obj.totalCharge;
        end
        
        function plot(obj)
            elements = obj.G.Nodes.type_symbol;
            xcoord = obj.x;
            ycoord = obj.y;
            zcoord = obj.z;
            for j=1:obj.nAtoms
                element = elements{j};
                switch upper(element)
                    case 'H'
                        opts = {'MarkerFaceColor',[1,1,1]*.7,'MarkerEdgeColor','none'};
                    case 'C'
                        opts = {'MarkerFaceColor',[1,.8,.6],'MarkerEdgeColor','none'};
                    case 'O'
                        opts = {'MarkerFaceColor',[.5,.5,1],'MarkerEdgeColor','none'};
                    case 'N'
                        opts = {'MarkerFaceColor',[1,.5,.5],'MarkerEdgeColor','none'};
                    case 'S'
                        opts = {'MarkerFaceColor',[1,.6,0],'MarkerEdgeColor','none'};
                    case 'P'
                        opts = {'MarkerFaceColor',[.2,.8,.2],'MarkerEdgeColor','none'};
                    otherwise
                        opts = {'MarkerEdgeColor',[0,0,0]};
                        warning('did not recognize the %dth element: %s',j,element);
                end
                plot3(xcoord(j),ycoord(j),zcoord(j),'o','MarkerSize',10,opts{:});hold on;
            end
            [s,t] = obj.G.findedge();
            isD = obj.isDouble;
            isC = obj.isCovalent;
            for j=1:obj.nBonds
                ind = [s(j),t(j)];
                if isD(j)
                    plot3(xcoord(ind),ycoord(ind),zcoord(ind),'k-','lineWidth',4);
                elseif ~isC(j)
                    plot3(xcoord(ind),ycoord(ind),zcoord(ind),'k--');
                else
                    plot3(xcoord(ind),ycoord(ind),zcoord(ind),'k-');
                end
            end
            
            xlabel('X');ylabel('Y');zlabel('Z');
            title(obj.name);
            axis equal
            set(gca,'FontSize',14);
            drawnow;
        end
        
        function [tt,ss,dd] = distanceRestraints(obj)
            [tt,ss,dd] = findDistanceRestraints(obj.G);
        end
        
        function [Icoh,Iincoh,Ibond] = scatteringFactors(obj,s,mode)
            
            if nargin<3 || isempty(mode)
                mode = 'direct';
            end
            
            switch lower(mode)
                case 'interp'
                    scalc = linspace(0,max(s(:))*1.01,100);
                    
                    [Icoh,Iincoh] = calculateScatteringFactors(scalc,...
                        obj.atomInventory.Type,obj.atomInventory.Number);
                    
                    % interpolate
                    Icoh = interp1(scalc,Icoh,s,'spline');
                    Iincoh = interp1(scalc,Iincoh,s,'spline');
                    
                    if nargout == 3
                        % calculate bonding correction
                        Ibond = calculateBondingInterferenceTerm(scalc,...
                            obj.G,obj.atomInventory.Type);
                        
                        % interpolate
                        Ibond = interp1(scalc,Ibond,s,'spline');
                    end
                case 'direct'
                    % default behavior
                    [Icoh,Iincoh] = calculateScatteringFactors(s,...
                        obj.atomInventory.Type,obj.atomInventory.Number);
                    
                    if nargout == 3
                        Ibond = calculateBondingInterferenceTerm(s,...
                            obj.G,obj.atomInventory.Type);
                    end
                otherwise
                    error('did not recognize mode');
            end
        end
    end
end


function [Icoh,Iincoh] = calculateScatteringFactors(s,atomType,atomNumber)

% 2) calculate atomic scattering factors
Elements = cellfun(@model.atom.Element,atomType,'UniformOutput',false);
Elements = [Elements{:}];
sf = [Elements.ScatteringFactor];
allf = arrayfun(@(a) f_coh(a,s),sf,'UniformOutput',false);
allincoh = arrayfun(@(a) I_incoh(a,s),sf,'UniformOutput',false);

% calculate Icoh, Iincoh
Icoh = zeros(size(s));
Iincoh = zeros(size(s));
for k=1:length(Elements)
    Icoh = Icoh + atomNumber(k)*allf{k}.^2;
    Iincoh = Iincoh + atomNumber(k)*allincoh{k};
end

end


function Ibond = calculateBondingInterferenceTerm(s,G0,atomType)

if size(G0.Edges,1) == 0 % for molecules with no bonds
    Ibond = zeros(size(s));
    return;
end

[tt,ss,dd] = findDistanceRestraints(G0);

% 2) calculate atomic scattering factors
Elements = cellfun(@model.atom.Element,atomType,'UniformOutput',false);
Elements = [Elements{:}];
sf = [Elements.ScatteringFactor];
allf = arrayfun(@(a) f_coh(a,s),sf,'UniformOutput',false);
[~,eInd] = ismember(G0.Nodes.type_symbol,atomType);

% 3) calculate Ibond
Ibond = zeros(size(s));

fprintf(1,'calculating Ibond for %d atom pairs\n',length(tt));

for k=1:length(tt)
    f1 = allf{eInd(tt(k))};
    f2 = allf{eInd(ss(k))};
    Ibond = Ibond + 2*f1.*f2.*mysinc(2*s*dd(k)); % q*dd(k)/pi = 2*pi*s*dd(k)/pi = 2*s*dd(k)
end

end

function y = mysinc(x)
% recreate the MATLAB sinc function 
%
% sinc = sin(pi*x)/(pi*x)
x = x * pi;
y = sin(x)./x;
y(x==0) = 1;

end


function [tt,ss,dd] = findDistanceRestraints(G0)

initialWeights = double(~(G0.Edges.isDouble | G0.Edges.isAromatic));

xx = G0.Nodes.x;
yy = G0.Nodes.y;
zz = G0.Nodes.z;
tt = [];
ss = [];
dd = [];
for k=1:G0.numnodes
    
    G0.Edges.Weight = initialWeights; % reset weights
    
    % find all edges joining node k with its neighbors.
    neighborEdges = findedge(G0,k,neighbors(G0,k));
    
    % set weights for these edges to zero
    G0.Edges.Weight(neighborEdges) = 0;
    
    % find all nodes within an edge distance of 1
    nn = nearest(G0,k,1);
    
    % ignore nodes with index less than k (avoids double-counting)
    nn = nn(nn>k);
    
    dx = xx(nn)-xx(k);
    dy = yy(nn)-yy(k);
    dz = zz(nn)-zz(k);
    tt = [tt; k*ones(size(nn))];
    ss = [ss; nn];
    dd = [dd; sqrt(dx.*dx + dy.*dy + dz.*dz)];
end

end