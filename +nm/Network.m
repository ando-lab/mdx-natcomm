classdef Network
    % NETWORK - crystalline elastic network model
    
    properties % I will remove these eventually. They are no longer needed, since they can be passed as arguments to Hessian()
        springType = 'Gaussian'
        springConstants = 1
    end
    properties(SetAccess = private)
        Edges % Table of cell edges, as produced by Cell.contactSearch()
        Cell  % model of the unit cell
        Nodes
    end
    
    methods
        function obj = Network(Edges,Cell,springType,springConstants)
            % NETWORK 
            assert(nargin >= 2)
            obj.Edges = Edges;
            obj.Cell = Cell;
            obj.Nodes = obj.getNodes;
            if nargin >= 3
                warning('assignment of spring type in constructor will be deprecated');
                obj.springType = springType;
            end
            if nargin >= 4
                warning('assignment of spring constant in constructor will be deprecated');
                obj.springConstants = springConstants;
            end
        end
        
        function T = getNodes(obj)
            Ops = obj.Cell.UnitCellOperators;
            B = symm.AffineTransformation(obj.Cell.Basis.orthogonalizationMatrix,[0;0;0]);
            
            r = [obj.Cell.AsymmetricUnit.r];
             
            r = inv(B)*r; % convert to fractional coordinates
            
            [n1,n2,n3] = ndgrid(-1:1,-1:1,-1:1);
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
            
            T = table(r1,r2);
        end
        
        function [Kab,Kab_perp,Kab_parallel] = Hessian(obj,springTypeOverride,springConstantsOverride)
            
            if nargin >= 2
                springMode = springTypeOverride;
            else
                springMode = obj.springType;
            end
            
            if nargin >= 3
                k = springConstantsOverride;
            else
                k = obj.springConstants;
            end
            
            % use ASU coordinate convention
            P = obj.Cell.tl2uxyz();
            
            if isempty(k)
                k = 1;
            end
            if size(k,1)==1
                k = repmat(k,size(obj.Edges,1),1);
            end
            
            switch lower(springMode)
                case {'gaussian','gauss'}
                    assert(size(k,2)==1);
                    gamma_perp = k;
                    gamma_parallel = k;
                case {'parallel'}
                    assert(size(k,2)==1);
                    gamma_perp = 0*k;
                    gamma_parallel = k;
                case {'hybrid'}
                    assert(size(k,2)==2);
                    gamma_perp = k(:,1);
                    gamma_parallel = k(:,2);
                otherwise
                    error('unrecognized spring type');
            end
           
            edgeMat = edgeMatrix(obj);
            
            numCells = numel(edgeMat);
            [n1,n2,n3] = ndgrid(-1:1,-1:1,-1:1);
            n123 = [n1(:),n2(:),n3(:)];
            [~,ixself] = ismember([0,0,0],n123,'rows'); % (= 14)
            assert(size(n123,1)==numCells); % just to make sure (= 27)
            
            E = cell(numCells,1);
            for n = 1:numCells
                E{n} = kron(edgeMat{n},speye(3,3))*P;
            end
            
            % calculate unit vectors along each bond
            N = obj.Nodes;
            v12 = N.r2 - N.r1;
            d12 = sqrt(sum(v12.*v12,2));
            v12 = v12.*repmat(1./d12,1,3);

            assert(size(v12,2)==3);
            [r,c] = ndgrid(1:size(v12,1),1:3);
            c2 = c + (r-1)*3;
            u = sparse(r(:),c2(:),v12(:));

            Kab = cell(numCells,1);
            Kab_perp = cell(numCells,1);
            Kab_parallel = cell(numCells,1);

            kmat_perp = (diag(sparse(kron(gamma_perp(:)',[1,1,1]))) - u'*diag(sparse(gamma_perp))*u);
            kmat_parallel = (u'*diag(sparse(gamma_parallel))*u);

            for n = 1:numCells
                
                Kab_perp{n} = E{ixself}'*kmat_perp*E{n};
                Kab_parallel{n} = E{ixself}'*kmat_parallel*E{n};
                Kab{n} = Kab_perp{n} + Kab_parallel{n};
            end
            
        end
        
        
        function edgeMat = edgeMatrix(obj)
            [n1,n2,n3] = ndgrid(-1:1,-1:1,-1:1);
            n123 = [n1(:),n2(:),n3(:)];
            [~,cellindex] = ismember(obj.Edges.c2,n123,'rows');

            numCells = size(n123,1); % 27
            numOps = numel(obj.Cell.UnitCellOperators);
            numAtoms = numel([obj.Cell.AsymmetricUnit.x]);
            numEdges = size(obj.Edges,1);

            % Order: all atoms of ASU1, then all atoms of ASU2, etc
            ao1 = sub2ind([numAtoms,numOps],obj.Edges.a1, obj.Edges.o1);
            ao2 = sub2ind([numAtoms,numOps],obj.Edges.a2, obj.Edges.o2);

            % group indices by cell index
            ao2Cell = cell(numCells,1);
            edgeID = cell(numCells,1);
            for n=1:numCells
                isCell = cellindex==n;
                ao2Cell{n} = ao2(isCell);
                edgeID{n} = find(isCell);
            end

            % define a function to calculate the edge matrix, for edges
            % connecting nodes i and j
            edgeMatFun = @(e,j) sparse(e,j(:)',ones(1,numel(j)),... % values
                numEdges,numOps*numAtoms);
            
            [~,ixself] = ismember([0,0,0],n123,'rows');
            edgeMat = cellfun(@(e,j) -1*edgeMatFun(e,j),edgeID,ao2Cell,'UniformOutput',false);
            edgeMat{ixself} = edgeMat{ixself} + edgeMatFun(1:numEdges,ao1);
        end

    end
end

