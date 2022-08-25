classdef ElasticNetwork_v2 < util.propertyValueConstructor
    %ElasticNetwork - calculate the Hessian matrix for TLS lattice ENM
    
    properties
        Cell % model of the unit cell (nm.Cell_v2 type)
        Edges % Table of cell edges, as produced by Cell.contactSearch()
    end
    properties(Constant)
        G_n = latt.PeriodicGrid([3,3,3],[-1,-1,-1],[3,3,3]);
    end
    
    methods
        function obj = ElasticNetwork_v2(varargin)
            %ElasticNetwork 
            obj@util.propertyValueConstructor(varargin{:});
        end
        
        function EdgeInfo = get_edge_info(obj)
            % information about edges to use in grouping
            
            E = obj.Edges;
            EdgeInfo = E(:,{'a1','a2','interface'});
            EdgeInfo.interface = abs(EdgeInfo.interface);
            
            ifswitch = EdgeInfo.a1 > EdgeInfo.a2;
            EdgeInfo.a1(ifswitch) = E.a2(ifswitch);
            EdgeInfo.a2(ifswitch) = E.a1(ifswitch);
            
            asu = obj.Cell.AsymmetricUnit;
            ix = cell2mat(arrayfun(@(g,n) n*ones(numel(g.x),1),asu(:),(1:numel(asu))','Uni',0));
            EdgeInfo.g1 = ix(EdgeInfo.a1);
            EdgeInfo.g2 = ix(EdgeInfo.a2);
        end
        
        % TODO: functions to parameterize the Hessian
        
        function [param2k,p0] = parameterize(obj,groupType,k0,userGroup)
            if nargin < 3 || isempty(k0)
                k0 = 1;
            end
            if size(k0,1)==1
                k0 = repmat(k0,size(obj.Edges,1),1);
            end
            
            useUserGroup = nargin > 3 && ~isempty(userGroup);
            
            EdgeInfo = obj.get_edge_info();
            
            if useUserGroup
                EdgeInfo.g1 = userGroup(EdgeInfo.a1);
                EdgeInfo.g2 = userGroup(EdgeInfo.a2);
            end
            
            switch lower(groupType)
                case 'global'
                    param2k = @(p) p;
                    p0 = mean(k0,1);
                case 'uniqueuser'
                    [~,ind0,ix] = unique(userGroup,'rows');
                    param2k = @(p) p(ix,:);
                    p0 = k0(ind0,:);
                case 'unique'
                    [~,ind0,ix] = unique(EdgeInfo(:,{'interface','a1','a2'}),'rows');
                    param2k = @(p) p(ix,:);
                    p0 = k0(ind0,:);
                case 'uniquegrouped'
                    [~,ind0,ix] = unique(EdgeInfo(:,{'interface','g1','g2'}),'rows');
                    param2k = @(p) p(ix,:);
                    p0 = k0(ind0,:);
                case 'interface'
                    [~,ind0,ix] = unique(EdgeInfo(:,{'interface'}),'rows');
                    param2k = @(p) p(ix,:);
                    p0 = k0(ind0,:);
                case 'group'
                    g1 = EdgeInfo.g1;
                    g2 = EdgeInfo.g2;
                    ngroups = numel(obj.Cell.AsymmetricUnit);

                    param2k = @(p) sqrt(p(g1,:).*p(g2,:)); % y is the coupling strength
                    p0 = repmat(sqrt(mean(k0.^2,1)),ngroups,1);
                case 'groupmin' % same as group, non-bonded groups are removed
                    g1 = EdgeInfo.g1;
                    g2 = EdgeInfo.g2;
                    g0 = unique([g1;g2]);
                    [~,g1] = ismember(g1,g0);
                    [~,g2] = ismember(g2,g0);
                    
                    ngroups = numel(g0);

                    param2k = @(p) sqrt(p(g1,:).*p(g2,:)); % y is the coupling strength
                    p0 = repmat(sqrt(mean(k0.^2,1)),ngroups,1);
                otherwise
                    error('did not recognize grouping type');
            end
            
        end
        
        function V = Hessian(obj,springType,k)
            switch lower(springType)
                case {'gauss','gaussian'}
                    V = obj.Hessian_Gauss(k);
                case {'par','parallel'}
                    V = obj.Hessian_parallel(k);
                case {'hybrid'}
                    assert(size(k,2)==2)
                    V = obj.Hessian_hybrid(k(:,1),k(:,2));
                otherwise
                    error('springType not recognized');
            end
        end
        
        function V = Hessian_hybrid(obj,kperp,kpar)
            assert(nargin==3 & ~isempty(kperp) & all(size(kperp)==size(kpar)));
            if numel(kperp)==1
                kperp = repmat(kperp,size(obj.Edges,1),1);
                kpar = repmat(kpar,size(obj.Edges,1),1);
            end
            E = obj.E();
            E_par = obj.E2Eparallel(E);
            E_perp = E_par;
            
            gamma = kron(sqrt(kperp(:)),[1;1;1]);
            gamma_par = sqrt(kpar(:));
            gamma_perp = sqrt(kperp(:));
            
            for n=1:numel(E)
                E{n} = gamma.*E{n};
                E_par{n} = gamma_par.*E_par{n};
                E_perp{n} = gamma_perp.*E_perp{n};
            end
            
            V = cell(size(E));
            for n = 1:numel(E)
                V{n} = E{2,2,2}'*E{n} + E_par{2,2,2}'*E_par{n} - E_perp{2,2,2}'*E_perp{n};
            end
            %V = obj.symavg(V);
        end
        
        function V = Hessian_parallel(obj,k)
            if nargin < 2 || isempty(k)
                k = 1;
            end
            if numel(k)==1
                k = repmat(k,size(obj.Edges,1),1);
            end
            E = obj.E2Eparallel(obj.E);
            gamma = sqrt(k(:));
            for n=1:numel(E)
                E{n} = gamma.*E{n};
            end
            V = cell(size(E));
            for n = 1:numel(E)
                V{n} = E{2,2,2}'*E{n};
            end
            %V = obj.symavg(V);
        end
        
        function V = Hessian_Gauss(obj,k)
            if nargin < 2 || isempty(k)
                k = 1;
            end
            if numel(k)==1
                k = repmat(k,size(obj.Edges,1),1);
            end
            E = obj.E();
            gamma = kron(sqrt(k(:)),[1;1;1]);
            for n=1:numel(E)
                E{n} = gamma.*E{n};
            end
            V = cell(size(E));
            for n = 1:numel(E)
                V{n} = E{2,2,2}'*E{n};
            end
            %V = obj.symavg(V);
        end

        
        function E = E(obj)
            P = obj.Cell.tl2uxyz();
            edgeMat = obj.edgeMatrix();
            
            E = cell(obj.G_n.N);
            for n = 1:numel(E)
                E{n} = kron(edgeMat{n},speye(3,3))*P;
            end
        end
        
        function E = E2Eparallel(obj,E)
            
            v12 = obj.bondVectors();
            
            [r,c] = ndgrid(1:size(v12,1),1:3);
            c2 = c + (r-1)*3;
            u = sparse(r(:),c2(:),v12(:));
            
            for n = 1:numel(E)
                E{n} = u*E{n};
            end
        end
        
        function v12 = bondVectors(obj)
            % calculate unit vectors along each bond
            [r1,r2] = obj.nodes();
            v12 = r2-r1;
            d12 = sqrt(sum(v12.*v12,2));
            v12 = v12.*repmat(1./d12,1,3);
            assert(size(v12,2)==3);
        end
        
        function [r1,r2] = nodes(obj)
            
            r = obj.Cell.unitCellCoordinates();
            B = obj.Cell.B;
            Binv = inv(B);
            
            r = cellfun(@(rr) Binv*rr,r,'Uni',0); % convert to fractional coordinates
            
            numEdges = size(obj.Edges,1);
            numOps = numel(r);
            
            r1 = zeros(3,numEdges);
            r2 = zeros(3,numEdges);
            
            ind1 = accumarray(obj.Edges.o1,1:numEdges,[numOps,1],@(v) {v});
            ind2 = accumarray(obj.Edges.o2,1:numEdges,[numOps,1],@(v) {v});
            
            for j=1:numOps
                r1(:,ind1{j}) = r{j}(:,obj.Edges.a1(ind1{j}));
                r2(:,ind2{j}) = r{j}(:,obj.Edges.a2(ind2{j}));
            end
            
            r2 = r2 + obj.Edges.c2'; % add shift
            
            r1 = B*r1; % convert back to cartesian coordinates
            r2 = B*r2; 
            
            r1 = r1';
            r2 = r2';
            
        end
        
        function edgeMat = edgeMatrix(obj)
            [n1,n2,n3] = obj.G_n.grid();
            n123 = [n1(:),n2(:),n3(:)];
            [~,cellindex] = ismember(obj.Edges.c2,n123,'rows');

            numOps = numel(obj.Cell.UnitCellOperators);
            numAtoms = numel([obj.Cell.AsymmetricUnit.x]);
            numEdges = size(obj.Edges,1);

            % Order: all atoms of ASU1, then all atoms of ASU2, etc
            ao1 = sub2ind([numAtoms,numOps],obj.Edges.a1, obj.Edges.o1);
            ao2 = sub2ind([numAtoms,numOps],obj.Edges.a2, obj.Edges.o2);

            % group indices by cell index
            ao2Cell = cell(obj.G_n.N);
            edgeID = cell(obj.G_n.N);
            for n=1:prod(obj.G_n.N)
                isCell = cellindex==n;
                ao2Cell{n} = ao2(isCell);
                edgeID{n} = find(isCell);
            end

            % define a function to calculate the edge matrix, for edges
            % connecting nodes i and j
            edgeMatFun = @(e,j) sparse(e,j(:)',ones(1,numel(j)),... % values
                numEdges,numOps*numAtoms);
            
            edgeMat = cellfun(@(e,j) -1*edgeMatFun(e,j),edgeID,ao2Cell,'Uni',0);
            edgeMat{2,2,2} = edgeMat{2,2,2} + edgeMatFun(1:numEdges,ao1);
        end
    end
end


