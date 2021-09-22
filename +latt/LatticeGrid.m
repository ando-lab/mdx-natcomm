classdef LatticeGrid
    properties
        PeriodicGrid
        Basis
    end
    properties (Dependent)
        N
        ori
        delta
    end
    methods
        function obj = LatticeGrid(PeriodicGrid,Basis)
            obj.PeriodicGrid = PeriodicGrid;
            obj.Basis = Basis;
        end
        
        function [r1,r2,r3] = ind2lab(obj,i1,i2,i3)
            [x1,x2,x3] = obj.PeriodicGrid.ind2frac(i1,i2,i3);
            [r1,r2,r3] = obj.Basis.frac2lab(x1,x2,x3);
        end
        
        function [r1,r2,r3] = grid(obj)
            [x1,x2,x3] = obj.PeriodicGrid.grid();
            [r1,r2,r3] = obj.Basis.frac2lab(x1,x2,x3);
        end
        
        function [ix,iy,iz] = lab2ind(obj,r1,r2,r3,varargin)
            [x1,x2,x3] = obj.Basis.lab2frac(r1,r2,r3);
            [ix,iy,iz] = obj.PeriodicGrid.frac2ind(x1,x2,x3,varargin{:});
        end
        
        function val = get.N(obj)
            val = obj.PeriodicGrid.N;
        end
        function val = get.ori(obj)
            oriFrac = obj.PeriodicGrid.ori;
            [X,Y,Z] = obj.Basis.frac2lab(oriFrac(1),oriFrac(2),oriFrac(3));
            val = [X,Y,Z];
        end
        function val = get.delta(obj)
            deltaFrac = diag(obj.PeriodicGrid.delta);
            [X,Y,Z] = obj.Basis.frac2lab(deltaFrac(:,1),deltaFrac(:,2),deltaFrac(:,3));
            val = [X,Y,Z];
        end
        function M = invert(obj)
            M = obj;
            M.PeriodicGrid = obj.PeriodicGrid.invert;
            M.Basis = obj.Basis.invert;
        end
        function [F,L] = fft(obj,varargin)
            % extends PeriodicGrid.fft to normalize by cell volume

            % The volume integral of F is equal to f(x=0,y=0,z=0)
            [F,~] = obj.PeriodicGrid.fft(varargin{:});
            L = obj.invert;
            F = F/L.Basis.volume; % normalize
        end
        function [f,L] = ifft(obj,varargin)
            % extends PeriodicGrid.fft to normalize by cell volume

            % The volume integral of f is equal to F(h=0,k=0,l=0)
            [f,~] = obj.PeriodicGrid.ifft(varargin{:});
            L = obj.invert;
            f = f/L.Basis.volume; % normalize
        end
        function [k1,k2,k3] = sphkernel(obj,rmax)

            % brute force method: make a spherical grid of points
            vecedge = linspace(-rmax,rmax,20);
            [V1,V2,V3] = ndgrid(vecedge,vecedge,vecedge);
            isIncl = V1.*V1 + V2.*V2 + V3.*V3 <= rmax*rmax;
            V1 = V1(isIncl);
            V2 = V2(isIncl);
            V3 = V3(isIncl);
            
            % transform to fractional coordinates
            [v1,v2,v3] = obj.Basis.lab2frac(V1,V2,V3);
            
            % then transform to indices and take the maximum distance
            gridDelta = obj.PeriodicGrid.delta; % = P/N
            i1max = max(ceil(v1(:)/gridDelta(1)));
            i2max = max(ceil(v2(:)/gridDelta(2)));
            i3max = max(ceil(v3(:)/gridDelta(3)));
            
            % construct a cube of relative indices
            [k1,k2,k3] = ndgrid(-i1max:i1max,-i2max:i2max,-i3max:i3max);
            % transform to fractional delta coordinates
            c1 = k1*gridDelta(1);
            c2 = k2*gridDelta(2);
            c3 = k3*gridDelta(3);
            % transform to lab coordinates
            [C1,C2,C3] = obj.Basis.frac2lab(c1,c2,c3);
            % include only those points with R<rmax
            isIncl = C1.*C1 + C2.*C2 + C3.*C3 < rmax.*rmax;
            k1 = k1(isIncl);
            k2 = k2(isIncl);
            k3 = k3(isIncl);
        end
        
        function [varargout] = splat(obj,k1,k2,k3,distfun,addfun,initVal,X,Y,Z,varargin)
            
            % convert to fractional coordinates
            [x,y,z] = obj.Basis.lab2frac(X,Y,Z);
            
            % make a wrapper for distfun that converts from fractional to
            % lab coordinates
            function d = distfunwrapper(xx,yy,zz,varargin)
                [XX,YY,ZZ] = obj.Basis.frac2lab(xx,yy,zz);
                d = distfun(XX,YY,ZZ,varargin{:});
            end
            
            varargout = cell([1,nargout]);
            
            %for j=1:nargout % do you have to initialize varargout?
            %    varargout{j} = [];
            %end
            
            [varargout{:}] = obj.PeriodicGrid.splat(k1,k2,k3,@distfunwrapper,addfun,initVal,x,y,z,varargin{:});
        end
    end
end