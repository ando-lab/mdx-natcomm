classdef Group
    %GROUP - set of cartesian coordinates and (optional) masses for an atomic model.
    
    properties
        x
        y
        z
        mass = 1 % mass (default value = 1)
    end
    properties(Dependent = true)
        massvec % a vector of masses
        r
    end
    
    methods
        function obj = Group(x,y,z,mass)
            if nargin<3
                error('not enough arguments');
            end
            
            natoms = numel(x);
            assert(numel(y) == natoms);
            assert(numel(z) == natoms);
            
            obj.x = x(:)'; % row vector
            obj.y = y(:)';
            obj.z = z(:)';
            
            if nargin >= 4 % mass is defined
                assert(numel(mass)==natoms | numel(mass)==1);
                obj.mass = mass(:)';
            end

        end
        
        function val = get.r(obj)
            val = [obj.x;obj.y;obj.z];
        end
        
        function val = rmax(obj)
            % maximum distance from the center of mass
            r0 = com(obj);
            dx = [obj.x] - r0(1);
            dy = [obj.y] - r0(2);
            dz = [obj.z] - r0(3);
            
            r2 = dx.*dx + dy.*dy + dz.*dz;
            val = sqrt(max(r2));
        end
       
        function val = com(obj)
            % COM - calculate the center of mass of a Group, or an array of
            % Groups
            if numel(obj) == 1 % normal condition
                m = obj.massvec;
                val = sum(repmat(m,3,1).*[obj.x;obj.y;obj.z],2)/...
                        sum(m);
            else % run on an array of objects
                Mr = [0;0;0];
                M = 0;
                for j=1:numel(obj)
                    % 1) calculate total mass
                    Mj = sum(obj(j).massvec);
                    M = M + Mj;
                    Mr = Mr + Mj*obj(j).com;
                end
                val = Mr/M;
            end
        end
        
                
        function m = get.massvec(obj)
            % MASSVEC - return a vector of masses
            % (replicates a single mass for all atoms, if necessary)

            nAtoms = numel(obj.x);
            m = ones(1,nAtoms);
            if isempty(obj.mass)
                % do nothing. Mass = 1
            elseif numel(numel(obj.mass))==1
                m(:) = obj.mass;
            else % size(obj.x) == size(obj.m)
                m(:) = obj(n).mass(:);
            end
        end
        
        function P = tl2uxyz(obj)
            % TL2UXYZ - returns an operator (sparse matrix) that converts from generalized
            % rigid-body vibrational coordinates (t1,t2,t3,l1,l2,l3) to
            % cartesian coordinates (ux,uy,uz)
            
            cmat = @(r1,r2,r3) sparse([3,1,2,2,3,1],[2,3,1,3,1,2],[r1,r2,r3,-r1,-r2,-r3],3,3);
            
            numGroups = numel(obj);

            Pn = {};
            for n=1:numGroups
                thisobj = obj(n);
                Pn = [Pn; arrayfun(@(r1,r2,r3) kron(sparse(1,n,1,1,numGroups),[speye(3),-1*cmat(r1,r2,r3)]),...
                    thisobj.x(:),thisobj.y(:),thisobj.z(:),'UniformOutput',false)];
            end
            P = cell2mat(Pn);
        end
        
        function M = tlmass(obj)
            % TLMASS - mass tensor (sparse matrix) in rigid-body vibrational
            % coordinates (t1,t2,t3,l1,l2,l3)
            %
            % (also works on arrays of groups)
            
            m = [obj.massvec];
            
            P = tl2uxyz(obj);
            mdiag = diag(sparse(kron(m,[1,1,1])));
            M = P'*mdiag*P;
        end
    end
end

