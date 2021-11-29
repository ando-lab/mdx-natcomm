classdef ElasticNetwork < util.propertyValueConstructor
    %ElasticNetwork - calculate the Hessian matrix for TLS lattice ENM
    
    properties
        Cell  % model of the unit cell (nm.Cell type)
        Edges % Table of cell edges, as produced by Cell.contactSearch()
    end
    properties(Constant)
        G_n = latt.PeriodicGrid([3,3,3],[-1,-1,-1],[3,3,3]);
    end
    
    methods
        function obj = ElasticNetwork(varargin)
            %ElasticNetwork 
            obj@util.propertyValueConstructor(varargin{:});
        end
        
        % TODO: functions to parameterize the Hessian
        
        function [param2k,p0] = parameterize(obj,groupType,k0)
            if size(k0,1)==1
                k0 = repmat(k0,size(obj.Edges,1),1);
            end
            switch lower(groupType)
                case 'global'
                    param2k = @(p) p;
                    p0 = mean(k0,1);
                case 'unique'
                    T = obj.Edges;
                    ifswitch = T.interface < 1;
                    SpringInfo = T;
                    SpringInfo.a1(ifswitch) = T.a2(ifswitch);
                    SpringInfo.a2(ifswitch) = T.a1(ifswitch);
                    SpringInfo.interface(ifswitch) = -T.interface(ifswitch);

                    [~,ind0,ix] = unique(SpringInfo(:,{'interface','a1','a2'}),'rows');

                    param2k = @(p) p(ix,:);
                    p0 = k0(ind0,:);
                case 'interface'
                    
                    [~,ind0,ix] = unique(abs(obj.Edges.interface));
                    param2k = @(p) p(ix,:);
                    p0 = k0(ind0,:);
                case 'group'
                    asu = obj.Cell.AsymmetricUnit;
                    ix = cell2mat(arrayfun(@(g,n) n*ones(1,numel(g.x)),asu,1:numel(asu),'Uni',0));
                    g1 = ix(obj.Edges.a1);
                    g2 = ix(obj.Edges.a2);

                    param2k = @(p) sqrt(p(g1,:).*p(g2,:)); % y is the coupling strength
                    p0 = sqrt(mean(k0.^2,1));
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
            V = obj.symavg(V);
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
            V = obj.symavg(V);
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
            V = obj.symavg(V);
        end
        
        function V = symavg(obj,V)
            % compute symmetry average of V to account for small numerical errors
            [f1,f2,f3] = obj.G_n.grid();
            
            for j=1:13 % the first 13 terms is a non-redundant set
                [n1,n2,n3] = obj.G_n.frac2ind(f1(j),f2(j),f3(j));
                [p1,p2,p3] = obj.G_n.frac2ind(-f1(j),-f2(j),-f3(j));
                V{n1,n2,n3} = 0.5*(V{n1,n2,n3} + V{p1,p2,p3}');
                V{p1,p2,p3} = V{n1,n2,n3}';
            end
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
            Ops = obj.Cell.UnitCellOperators;
            B = symm.AffineTransformation(obj.Cell.Basis.orthogonalizationMatrix,[0;0;0]);
            
            r = [obj.Cell.AsymmetricUnit.r];
             
            r = inv(B)*r; % convert to fractional coordinates
            
            [n1,n2,n3] = obj.G_n.grid();
            n123 = [n1(:),n2(:),n3(:)];
            SuperCellOps = arrayfun(...
                @(a1,a2,a3) symm.SymmetryOperator(eye(3),[a1;a2;a3]),...
                n123(:,1),n123(:,2),n123(:,3));
            [~,cellindex] = ismember(obj.Edges.c2,n123,'rows');

            numEdges = size(obj.Edges,1);
                       
            Op1 = Ops(obj.Edges.o1);
            Op2 = Ops(obj.Edges.o2);
            S2 = SuperCellOps(cellindex);
            r1 = r(:,obj.Edges.a1);
            r2 = r(:,obj.Edges.a2);
            
            for j=1:numEdges
                r1(:,j) = Op1(j)*r1(:,j);
                r2(:,j) = S2(j)*Op2(j)*r2(:,j);
            end
            
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


