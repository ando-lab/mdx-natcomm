classdef InterpLin3
    % InterpLin3 - 3D linear interpolation model
    properties
        Nwx % = 5          % default number of samples
        Nwy % = 5          % default number of samples
        Nwz % = 5
        x % = repmat(linspace(0,1,11)',1,11)% points to interpolate at
        y % = repmat(linspace(0,1,11),11,1);% points to interpolate at
        z % = 0.5*ones(11,11);
        xmin % = 0          % range of interpolating parameters
        xmax % = 1 + eps
        ymin % = 0
        ymax % = 1 + eps
        zmin % = 0
        zmax % = 1 + eps
    end
    properties(Dependent = true)
        A  % interpolation operator: z(x) = A*w(x)
        wx  % vector of sample points
        wy
        wz
        dwx
        dwy
        dwz
        Bxy % 2D laplacian operator: laplacian(w(x,y,z),x,y) = Bxy*w
        Bx  % second derivative operator in x: D(w,x) = Bx*w
        By  % second derivative operator in y: D(w,y) = By*w
        Bz  % second derivative operator in z: D(w,z) = Bz*w
        B   % 3D laplacian operator: laplacian(w(x,y,z),x,y,z) = B*w
        N  % numel(x) ( == numel(y) == numel(z))
    end
    methods
        function obj = InterpLin3(varargin)
            if ~isempty(varargin)
                for j=1:2:length(varargin)
                    obj.(varargin{j}) = varargin{j+1};
                end
            end
        end
        function val = get.A(obj)
            [ind,ix,iy,iz,t,s,r] = obj.map();
            
            val = sparse(repmat(ind,1,8),...
                mysub2ind([obj.Nwx,obj.Nwy,obj.Nwz],...
                [ix,(ix+1),ix,(ix+1),ix,(ix+1),ix,(ix+1)],...
                [iy,iy,(iy+1),(iy+1),iy,iy,(iy+1),(iy+1)],...
                [iz,iz,iz,iz,(iz+1),(iz+1),(iz+1),(iz+1)]),...
                [(1-r).*(1-t).*(1-s),(1-r).*t.*(1-s),(1-r).*(1-t).*s,(1-r).*t.*s,...
                r.*(1-t).*(1-s),r.*t.*(1-s),r.*(1-t).*s,r.*t.*s],...
                obj.N,obj.Nwx*obj.Nwy*obj.Nwz);
        end
        function val = interp(obj,p)
            % faster calculation of A*p(:)
            [ind,ix,iy,iz,t,s,r] = obj.map();
            sz = [obj.Nwx,obj.Nwy,obj.Nwz];
            ind1 = mysub2ind(sz,ix,iy,iz);
            ind2 = mysub2ind(sz,ix+1,iy,iz);
            ind3 = mysub2ind(sz,ix,iy+1,iz);
            ind4 = mysub2ind(sz,ix+1,iy+1,iz);
            ind5 = mysub2ind(sz,ix,iy,iz+1);
            ind6 = mysub2ind(sz,ix+1,iy,iz+1);
            ind7 = mysub2ind(sz,ix,iy+1,iz+1);
            ind8 = mysub2ind(sz,ix+1,iy+1,iz+1);
            val = zeros(obj.N,1);
            val(ind) = p(ind1).*(1-r).*(1-t).*(1-s) + ...
                p(ind2).*(1-r).*t.*(1-s) + ...
                p(ind3).*(1-r).*(1-t).*s + ...
                p(ind4).*(1-r).*t.*s + ...
                p(ind5).*r.*(1-t).*(1-s) + ...
                p(ind6).*r.*t.*(1-s) + ...
                p(ind7).*r.*(1-t).*s + ...
                p(ind8).*r.*t.*s;
        end
        
        function [ind,ix,iy,iz,t,s,r] = map(obj)
            thisx = obj.x(:);
            [~,ix] = histc(thisx,[-Inf;obj.wx;Inf]);
            ix = ix-1;
            
            thisy = obj.y(:);
            [~,iy] = histc(thisy,[-Inf;obj.wy;Inf]);
            iy = iy-1;
            
            thisz = obj.z(:);
            [~,iz] = histc(thisz,[-Inf;obj.wz;Inf]);
            iz = iz-1;
            
            ind = (1:obj.N)';
            isInside = iy > 0 & iy < obj.Nwy & ix > 0 & ix < obj.Nwx & ...
                iz > 0 & iz < obj.Nwz;
            
            ix = ix(isInside);
            iy = iy(isInside);
            iz = iz(isInside);
            ind = ind(isInside);
            
            t = (thisx(isInside) - obj.wx(ix))*(1/obj.dwx);
            s = (thisy(isInside) - obj.wy(iy))*(1/obj.dwy);
            r = (thisz(isInside) - obj.wz(iz))*(1/obj.dwz);

        end
        
        function val = get.N(obj)
            assert(numel(obj.x)==numel(obj.y) & numel(obj.x)==numel(obj.z),...
                'x, y, and z must have the same number of elements');
            val = numel(obj.x);
        end
        function val = get.dwx(obj)
            val = (obj.xmax - obj.xmin)/(obj.Nwx-1);
        end
        function val = get.wx(obj)
            val = linspace(obj.xmin,obj.xmax,obj.Nwx)';
        end
        function val = get.dwy(obj)
            val = (obj.ymax - obj.ymin)/(obj.Nwy-1);
        end
        function val = get.wy(obj)
            val = linspace(obj.ymin,obj.ymax,obj.Nwy)';
        end
        function val = get.dwz(obj)
            val = (obj.zmax - obj.zmin)/(obj.Nwz-1);
        end
        function val = get.wz(obj)
            val = linspace(obj.zmin,obj.zmax,obj.Nwz)';
        end
        function val = get.Bxy(obj)
            [k,l,m] = ndgrid(1:obj.Nwx,1:obj.Nwy,1:obj.Nwz);
            
            right_neighbors = mysub2ind([obj.Nwx,obj.Nwy,obj.Nwz],...
                k(2:end,:,:),l(2:end,:,:),m(2:end,:,:));
            left_neighbors = mysub2ind([obj.Nwx,obj.Nwy,obj.Nwz],...
                k(1:(end-1),:,:),l(1:(end-1),:,:),m(1:(end-1),:,:));
            top_neighbors = mysub2ind([obj.Nwx,obj.Nwy,obj.Nwz],...
                k(:,1:(end-1),:),l(:,1:(end-1),:),m(:,1:(end-1),:));
            bottom_neighbors = mysub2ind([obj.Nwx,obj.Nwy,obj.Nwz],...
                k(:,2:end,:),l(:,2:end,:),m(:,2:end,:));
            
            L = sparse(...
                [left_neighbors(:)',right_neighbors(:)',...
                top_neighbors(:)',bottom_neighbors(:)'],...
                [right_neighbors(:)',left_neighbors(:)',...
                bottom_neighbors(:)',top_neighbors(:)'],...
                1,obj.Nwx*obj.Nwy*obj.Nwz,obj.Nwx*obj.Nwy*obj.Nwz);
            
            % weight according to number of neighbors
            L = diag(sparse(1./full(sum(L,2))))*L; 
            
            % center minus average of neighbors
            val = speye(obj.Nwx*obj.Nwy*obj.Nwz) - L; 
        end
        function val = get.B(obj) % 3D laplacian
            
            [k,l,m] = ndgrid(1:obj.Nwx,1:obj.Nwy,1:obj.Nwz);
            
            right_neighbors = mysub2ind([obj.Nwx,obj.Nwy,obj.Nwz],...
                k(2:end,:,:),l(2:end,:,:),m(2:end,:,:));
            left_neighbors = mysub2ind([obj.Nwx,obj.Nwy,obj.Nwz],...
                k(1:(end-1),:,:),l(1:(end-1),:,:),m(1:(end-1),:,:));
            top_neighbors = mysub2ind([obj.Nwx,obj.Nwy,obj.Nwz],...
                k(:,1:(end-1),:),l(:,1:(end-1),:),m(:,1:(end-1),:));
            bottom_neighbors = mysub2ind([obj.Nwx,obj.Nwy,obj.Nwz],...
                k(:,2:end,:),l(:,2:end,:),m(:,2:end,:));
            front_neighbors = mysub2ind([obj.Nwx,obj.Nwy,obj.Nwz],...
                k(:,:,1:(end-1)),l(:,:,1:(end-1)),m(:,:,1:(end-1)));
            back_neighbors = mysub2ind([obj.Nwx,obj.Nwy,obj.Nwz],...
                k(:,:,2:end),l(:,:,2:end),m(:,:,2:end));

            L = sparse(...
                [left_neighbors(:)',right_neighbors(:)',...
                top_neighbors(:)',bottom_neighbors(:)',...
                front_neighbors(:)',back_neighbors(:)'],...
                [right_neighbors(:)',left_neighbors(:)',...
                bottom_neighbors(:)',top_neighbors(:)',...
                back_neighbors(:)',front_neighbors(:)'],...
                1,obj.Nwx*obj.Nwy*obj.Nwz,obj.Nwx*obj.Nwy*obj.Nwz);
            % weight according to number of neighbors
            L = diag(sparse(1./full(sum(L,2))))*L; 
            
            % center minus average of neighbors
            val = speye(obj.Nwx*obj.Nwy*obj.Nwz) - L;
        end
        function val = get.Bx(obj)
            [k,l,m] = ndgrid(1:obj.Nwx,1:obj.Nwy,1:obj.Nwz);
            sz = [obj.Nwx,obj.Nwy,obj.Nwz];
            left  = mysub2ind(sz,k(1:(end-2),:,:),l(1:(end-2),:,:),m(1:(end-2),:,:));
            self  = mysub2ind(sz,k(2:(end-1),:,:),l(2:(end-1),:,:),m(2:(end-1),:,:));
            right = mysub2ind(sz,k(3:end,:,:),l(3:end,:,:),m(3:end,:,:));
            
            rows  = mysub2ind([obj.Nwx - 2,obj.Nwy,obj.Nwz],...
                k(1:(end-2),:,:),l(1:(end-2),:,:),m(1:(end-2),:,:));
            
            nrows = (obj.Nwx-2)*obj.Nwy*obj.Nwz;

            val = sparse(repmat(rows(:)',1,3),[left(:)',self(:)',right(:)'],...
                [-.5*ones(1,nrows),ones(1,nrows),-.5*ones(1,nrows)],...
                nrows,prod(sz));
        end
        function val = get.By(obj)
            [k,l,m] = ndgrid(1:obj.Nwx,1:obj.Nwy,1:obj.Nwz);
            sz = [obj.Nwx,obj.Nwy,obj.Nwz];
            left  = mysub2ind(sz,k(:,1:(end-2),:),l(:,1:(end-2),:),m(:,1:(end-2),:));
            self  = mysub2ind(sz,k(:,2:(end-1),:),l(:,2:(end-1),:),m(:,2:(end-1),:));
            right = mysub2ind(sz,k(:,3:end,:),l(:,3:end,:),m(:,3:end,:));
            
            rows  = mysub2ind([obj.Nwx,obj.Nwy-2,obj.Nwz],...
                k(:,1:(end-2),:),l(:,1:(end-2),:),m(:,1:(end-2),:));
            
            nrows = obj.Nwx*(obj.Nwy-2)*obj.Nwz;

            val = sparse(repmat(rows(:)',1,3),[left(:)',self(:)',right(:)'],...
                [-.5*ones(1,nrows),ones(1,nrows),-.5*ones(1,nrows)],...
                nrows,prod(sz));
        end
        function val = get.Bz(obj)
            [k,l,m] = ndgrid(1:obj.Nwx,1:obj.Nwy,1:obj.Nwz);
            sz = [obj.Nwx,obj.Nwy,obj.Nwz];
            left  = mysub2ind(sz,k(:,:,1:(end-2)),l(:,:,1:(end-2)),m(:,:,1:(end-2)));
            self  = mysub2ind(sz,k(:,:,2:(end-1)),l(:,:,2:(end-1)),m(:,:,2:(end-1)));
            right = mysub2ind(sz,k(:,:,3:end),l(:,:,3:end),m(:,:,3:end));
            
            rows  = mysub2ind([obj.Nwx,obj.Nwy,obj.Nwz-2],...
                k(:,:,1:(end-2)),l(:,:,1:(end-2)),m(:,:,1:(end-2)));
            
            nrows = obj.Nwx*obj.Nwy*(obj.Nwz-2);

            val = sparse(repmat(rows(:)',1,3),[left(:)',self(:)',right(:)'],...
                [-.5*ones(1,nrows),ones(1,nrows),-.5*ones(1,nrows)],...
                nrows,prod(sz));
        end
    end
end

function ind = mysub2ind(sz,x,y,z)
% faster than the built-in method (no input error checking)
%   note: ONLY WORKS ON 3D ARRAYS
    ind = x + sz(1)*(y-1) + sz(2)*sz(1)*(z-1);
end