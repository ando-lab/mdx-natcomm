classdef PeriodicGrid
    properties
        N % vector [N1, N2, N3] with the number of divisions along each direction
        ori = [0,0,0] % origin in fractional coordinates of element(1,1,1)
        P = [1,1,1] % vector [P1, P2, P3] with the period (in fractional coordinates) along each direction
    end
    properties(Dependent)
        delta % the grid spacing in fractional coordinates in each direction
    end
    methods
        function obj = PeriodicGrid(N,ori,P)
            obj.N = N;
            if nargin>=2
                obj.ori = ori;
            end
            if nargin>=3
                obj.P = P;
            end
        end
        function val = get.delta(obj)
            val = obj.P./obj.N;
        end
        function [ix,iy,iz] = frac2ind(obj,x,y,z,isPeriodic)
            % frac2ind - map fractional coordinates to indices in N1 x N2 x N3 array.
            % periodic boundary conditions apply in default case.
            if nargin==4 || isempty(isPeriodic)
                isPeriodic = true;
            end
            
            ix = 1 + round((x - obj.ori(1))*obj.N(1)/obj.P(1));
            iy = 1 + round((y - obj.ori(2))*obj.N(2)/obj.P(2));
            iz = 1 + round((z - obj.ori(3))*obj.N(3)/obj.P(3));
            
            if isPeriodic
                ix = 1 + mod(ix - 1,obj.N(1));
                iy = 1 + mod(iy - 1,obj.N(2));
                iz = 1 + mod(iz - 1,obj.N(3));
            end
        end
        function [x,y,z] = grid(obj)
            [n1,n2,n3] = ndgrid(1:obj.N(1),1:obj.N(2),1:obj.N(3));
            [x,y,z] = obj.ind2frac(n1,n2,n3);
        end
        function [x,y,z] = ind2frac(obj,ix,iy,iz)
            x = obj.ori(1) + (ix-1)*obj.P(1)/obj.N(1);
            y = obj.ori(2) + (iy-1)*obj.P(2)/obj.N(2);
            z = obj.ori(3) + (iz-1)*obj.P(3)/obj.N(3);
        end
        function M = invert(obj)
            M = obj;
            M.N = obj.N;
            M.P = 1./(obj.delta);
            M.ori = -floor(obj.N/2)./obj.P;
        end
        function [A,I] = splat(obj,k1,k2,k3,distfun,addfun,initVal,x,y,z,varargin)
            % SPLAT accumulates a function on the periodic grid within a
            % predefined kernel.
            %
            % Usage examples:
            %
            % Example 1 - Map of distance to closest atom within a neighborhood of
            % 27 grid points. Atom coordinates are given by x, y, z. r is
            % the van der waals radius of each atom.
            %
            %    [k1,k2,k3] = ndgrid(-1:1,-1:1,-1:1);
            %    distfun = @(x,y,z,r) sqrt(x.*x + y.*y + z.*z) - r
            %    addfun = @min
            %    initVal = Inf;
            %    A = obj.splat(k1,k2,k3,distfun,addfun,initVal,x,y,z,r)
            %
            % NOTE (new feature as of sept 18, 2021)
            % If addfun is @le, @leq, @ge, or @geq, splat returns a second variable I
            % which is a map of integer indices for each voxel.
            %
            % Example 2 - Map of "electron density" with a 1-Gaussian form
            % factor and a cutoff of nearest neighbor grid points. b is a
            % parameter controlling the gaussian width which may be
            % different for each atom.
            %
            %    [k1,k2,k3] = ndgrid(-1:1,-1:1,-1:1);
            %    ffun = @(x,y,z,b) exp(-b.*(x.*x + y.*y + z.*z))
            %    addfun = @plus
            %    initVal = 0;
            %    A = obj.splat(k1,k2,k3,ffun,addfun,initVal,x,y,z,b)
            
            % initialize output matrix
            A = zeros(obj.N);
            A(:) = initVal;
            
            if any(strcmp(functions(addfun).function,{'le','ge','leq','geq'}))
                outputMode = 'index';
                % there are two outputs. Initialize the second one.
                I = zeros(obj.N);
            else
                outputMode = 'simple';
            end

            npts = numel(x);
            
            [ix,iy,iz] = obj.frac2ind(x,y,z,false);
            [x0,y0,z0] = obj.ind2frac(ix,iy,iz);
            dx = x - x0;
            dy = y - y0;
            dz = z - z0;
            kx = k1*obj.P(1)/obj.N(1);
            ky = k2*obj.P(2)/obj.N(2);
            kz = k3*obj.P(3)/obj.N(3);
            
            for j=1:npts
                k1shift = mod(k1 + ix(j) - 1,obj.N(1)) + 1;
                k2shift = mod(k2 + iy(j) - 1,obj.N(2)) + 1;
                k3shift = mod(k3 + iz(j) - 1,obj.N(3)) + 1;
                kshiftind = sub2ind(obj.N,k1shift,k2shift,k3shift);
                args = cellfun(@(c) c(j,:),varargin,'uniformOutput',false);
                kval = distfun(kx - dx(j),ky - dy(j),kz - dz(j),args{:});
                
                addFunOut = addfun(kval(:),A(kshiftind(:)));
                
                switch outputMode
                    case 'index'
                        I(kshiftind(addFunOut)) = j;
                        A(kshiftind(addFunOut)) = kval(addFunOut);
                    case 'simple'
                        A(kshiftind(:)) = addFunOut;
                end
            end
            
        end
        function [F,M] = fft(obj,f)
            % FFT fast Fourier transform
            %
            %  [F,M] = fft(obj,f)
            
            assert(size(f,1)==obj.N(1) & ...
                size(f,2)==obj.N(2) & ...
                size(f,3)==obj.N(3),...
                'size of input array does not match Map specification');
            
            M = obj.invert;
            
            if ~isempty(obj.ori) && any(obj.ori)
            
                [x,y,z] = M.ind2frac(1:M.N(1),1:M.N(2),1:M.N(3));
                [ph,pk,pl] = ndgrid(...
                    exp(-2i*pi*obj.ori(1)*x),...
                    exp(-2i*pi*obj.ori(2)*y),...
                    exp(-2i*pi*obj.ori(3)*z));
            
                phaseFactor = ph.*pk.*pl;
                
            else
                phaseFactor = 1;
            end
            
            f = f*prod(obj.delta); % normalize so that mapfft, mapifft are a transform pair.
            
            F = fftn(f); % calculate Fourier transform
            F = fftshift(F);
            F = F.*phaseFactor;
        end
        
        function [f,M] = ifft(obj,F)
            % IFFT fast inverse Fourier transform
            %
            %   [f,M] = ifft(obj,F)
            
            assert(size(F,1)==obj.N(1) & ...
                size(F,2)==obj.N(2) & ...
                size(F,3)==obj.N(3),...
                'size of input array does not match Map specification');
            
            M = obj.invert;
            
            if ~isempty(obj.ori) && any(obj.ori)
                [x,y,z] = M.ind2frac(1:M.N(1),1:M.N(2),1:M.N(3));
                [ph,pk,pl] = ndgrid(...
                    exp(2i*pi*obj.ori(1)*x),...
                    exp(2i*pi*obj.ori(2)*y),...
                    exp(2i*pi*obj.ori(3)*z));
            
                phaseFactor = ph.*pk.*pl;
                
            else
                phaseFactor = 1;
            end
            
            f = ifftn(F); % calculate Fourier transform
            f = fftshift(f);
            f = f.*phaseFactor;
            f = f/prod(M.delta); % normalize so (integral of P) = sum(p(:))*prod(delta(1)) = F(0)
        end
    end
end