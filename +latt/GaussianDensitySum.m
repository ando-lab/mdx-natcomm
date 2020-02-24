classdef GaussianDensitySum
    %GAUSSIANDENSITYSUM
    %
    % Electron density and scattering factor calculations for latt.GaussianAtom
    
    properties(SetAccess = immutable)
        a
        V
    end
    properties(SetAccess = immutable, Hidden = true)
        isIsotropic % T/F
        Vinv % = inv(V)
        ap % = obj.a(j)/sqrt((2*pi*V)^3);
    end
    properties(Dependent = true)
        minEig
    end
    
    methods
        function obj = GaussianDensitySum(a,V)
            %GAUSSIANDENSITYSUM
            
            assert(nargin==2) % a and V are required
            
            N = numel(a);
            a = a(:)'; % make sure it's 1 x N
            
            if numel(V)==N
                obj.isIsotropic = true;
            else
                obj.isIsotropic = false;
                assert(numel(V) == 9*N && size(V,1)==3 && size(V,2)==3 && size(V,3)==N);
            end
            
            if obj.isIsotropic % V is a vector
                V = V(:)'; % make sure it's 1 x N
                Vinv = 1./V;
                ap = a./sqrt((2*pi*V).^3);
            else % V is a 3x3xN matrix
                Vinv = zeros(3,3,N);
                ap = zeros(1,N);
                for j=1:N
                    ap(j) = a(j)/sqrt((2*pi)^3*det(V(:,:,j)));
                    % make sure V is symmetric
                    V(:,:,j) = 0.5*(V(:,:,j) + V(:,:,j)');
                    Vinv(:,:,j) = inv(V(:,:,j));
                end
            end
            
            obj.a = a;
            obj.V = V;
            obj.Vinv = Vinv;
            obj.ap = ap;
            
        end
        
        function GA = addU(obj,U)
            V0 = obj.V;
            if obj.isIsotropic && numel(U)==1
                % simple case
                Vnew = V0 + U;
            else % handle more complex cases
                
                if numel(U)==1
                    U = eye(3)*U;
                end
                if obj.isIsotropic
                    Vnew = zeros(3,3,numel(V0));
                    for j=1:numel(V0)
                        Vnew(:,:,j) = V0(j)*eye(3) + U;
                    end
                else
                    Vnew = zeros(3,3,size(V0,3));
                    for j=1:size(V0,3)
                        Vnew(:,:,j) = V0(:,:,j) + U;
                    end
                end
            end
            
            GA = latt.GaussianDensitySum(obj.a,Vnew);
        end
        
        function rho = electronDensity(obj,x,y,z)
            
            rho = zeros(size(x));
            r = [x(:),y(:),z(:)]';
            
            N = numel(obj.a); % = numel(obj.b)
            
            if obj.isIsotropic % V is a vector
                r2 = sum(r.*r,1);
                for j=1:N
                    thisrho = obj.ap(j)*exp(-0.5*r2*obj.Vinv(j));
                    rho(:) = rho(:) + thisrho';
                end
            else % V is a 3x3xN matrix
                for j=1:N
                    thisrho = obj.ap(j)*exp(-0.5*dot(r,obj.Vinv(:,:,j)*r));
                    rho(:) = rho(:) + thisrho';
                end
            end
            
        end
        
        function F = scatteringAmplitude(obj,sx,sy,sz)
            F = zeros(size(sx));
            s = [sx(:),sy(:),sz(:)]';
            
            N = numel(obj.a);
            
            if obj.isIsotropic % V is a vector
                s2 = sum(s.*s,1);
                for j=1:N
                    thisF = obj.a(j)*exp(-2*pi^2*obj.V(j)*s2);
                    F(:) = F(:) + thisF';
                end
            else % V is a 3x3xN matrix
                for j=1:N
                    thisF = obj.a(j)*exp(-2*pi^2*dot(s,obj.V(:,:,j)*s));
                    F(:) = F(:) + thisF';
                end
            end
        end
        
        function d = get.minEig(obj)
            % return the minimum eigenvalue for all matrices V
            % (or, if is isotropic, returns the minimum value of V)
            if obj.isIsotropic
                d = min(obj.V);
            else
                d = Inf;
                for j=1:numel(obj.a)
                    v = eig(obj.V(:,:,j),'vector');
                    d = min([d,v(:)']);
                end
            end
        end
        
        function Vadd = fix(obj,dmin)
            Vadd = zeros(size(obj.V));
            
            if obj.minEig > dmin
                % don't need to fix anything
                return;
            end
            
            if obj.isIsotropic
                for j = 1:numel(obj.a)
                    if obj.V(j) > dmin
                        Vadd(j) = dmin - obj.V(j);
                    end
                end
            else
                for j=1:numel(obj.a)
                    [v,d] = eig(obj.V(:,:,j),'vector');
                    if any(d < dmin)
                        dp = lsqnonlin(@(x) [100*sum(x),x],...
                            [0,0,0],dmin-d,[Inf,Inf,Inf],...
                            optimset('Display','off'));
                        Vadd(:,:,j) = v*diag(dp)*v';
                    end
                end
            end
        end
        
        function newobj = transform(obj,A)
            if isa(A,'symm.AffineTransformation')
                A = A.r; % get the rotation part
            end
            assert(all(size(A)==[3,3]));
            
            if obj.isIsotropic && (abs(abs(det(A))-1) < 10*eps)
                newobj = obj;
                return;
            end
            
            Vnew = zeros(3,3,numel(obj.a));
            
            if obj.isIsotropic
                AA = (A*A');
                for j=1:numel(obj.a)
                    Vnew(:,:,j) = AA*obj.V(j);
                end
            else
                for j=1:numel(obj.a)
                    Vnew(:,:,j) = A*obj.V(:,:,j)*A';
                end
            end
            newobj = latt.GaussianDensitySum(obj.a,Vnew);
        end
    end
    methods(Static)
        
        function G = approximateEllipsoid(U,density)

            % first, make a sphere with radius 1 and unit mass.
            % then, rescale the coefficients so that the equivalent "U" is
            % the identity matrix.
            
            sphereVolume = (4*pi/3);
            ellipsoidVolume = sphereVolume*sqrt(det(U)/.2^3);
            
            G0 = latt.GaussianDensitySum.approximateSphere(1,1/sphereVolume);
            
            u0 = sum(G0.V.*G0.a)/sum(G0.a); % diagonal element of U0
            
            v = G0.V/u0; % scale B so that diagonal elements of U are 1
            
            N = numel(v);
            V = zeros(3,3,N);
            for j=1:N
                V(:,:,j) = v(j)*U;
            end
            
            a = G0.a*density*ellipsoidVolume;
            
            G = latt.GaussianDensitySum(a,V);
        end
        
        function G = approximateSphere(radius,density)
            a0 = [-0.572213, 3.488039, 77.859464, -79.775291];
            b0 = [2.521880, 4.003510, 7.394241, 7.165021];
            
            % the above are coefficients to approximate a sphere with radius 1 and
            % electron density 3/(4*pi) = 0.2387 (one electron)
            
            b = b0*radius.^2;
            a = a0*radius.^3*density*4*pi/3;
            
            V = b/(8*pi^2);
            
            G = latt.GaussianDensitySum(a,V);
            
        end
        
    end
end

