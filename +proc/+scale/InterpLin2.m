classdef InterpLin2
    % InterpLin2 - 2D linear interpolation model
    properties
        Nwx % = 5          % default number of samples
        Nwy % = 5          % default number of samples
        x % = repmat(linspace(0,1,11)',1,11)% points to interpolate at
        y % = repmat(linspace(0,1,11),11,1);% points to interpolate at
        xmin % = 0          % range of interpolating parameters
        xmax % = 1.000001
        ymin % = 0
        ymax % = 1.000001
    end
    properties(Dependent = true)
        A  % interpolation operator: z(x) = A*w(x)
        wx  % vector of points sampling the peak
        wy
        dwx
        dwy
        B  % second derivative operator: laplacian(w(x,y)) = B*w(x,y)
        Bx % second derivative with respect to x: D(w,x) = Bx*w
        By
        N  % length(x) ( == length(y) )
    end
    methods
        function obj = InterpLin2(varargin)
            if ~isempty(varargin)
                for j=1:2:length(varargin)
                    obj.(varargin{j}) = varargin{j+1};
                end
            end
        end
        function val = get.A(obj)
            [ind,ix,iy,t,s] = obj.map();
            
            val = sparse([ind,ind,ind,ind],...
                mysub2ind([obj.Nwx,obj.Nwy],...
                [ix,(ix+1),ix,(ix+1)],[iy,iy,(iy+1),(iy+1)]),...
                [(1-t).*(1-s),t.*(1-s),(1-t).*s,t.*s],...
                obj.N,obj.Nwx*obj.Nwy);
        end
        function val = interp(obj,p)
            % faster calculation of A*p(:)
            [ind,ix,iy,t,s] = obj.map();
            val = zeros(obj.N,1);
            sz = [obj.Nwx,obj.Nwy];
            ind1 = mysub2ind(sz,ix,iy);
            ind2 = mysub2ind(sz,ix+1,iy);
            ind3 = mysub2ind(sz,ix,iy+1);
            ind4 = mysub2ind(sz,ix+1,iy+1);
            
            val(ind) =  p(ind1).*(1-t).*(1-s) + ...
                        p(ind2).*t.*(1-s) + ...
                        p(ind3).*(1-t).*s + ...
                        p(ind4).*t.*s;
        end
        function [ind,ix,iy,t,s] = map(obj)
            thisx = obj.x(:);
            [~,ix] = histc(thisx,[-Inf;obj.wx;Inf]);
            ix = ix-1;
            
            thisy = obj.y(:);
            [~,iy] = histc(thisy,[-Inf;obj.wy;Inf]);
            iy = iy-1;
            
            ind = (1:obj.N)';
            isInside = iy > 0 & iy < obj.Nwy & ix > 0 & ix < obj.Nwx;
            
            ix = ix(isInside);
            iy = iy(isInside);
            ind = ind(isInside);
            
            t = (thisx(isInside) - obj.wx(ix))*(1/obj.dwx);
            s = (thisy(isInside) - obj.wy(iy))*(1/obj.dwy);
        end
        
        function val = get.N(obj)
            assert(numel(obj.x)==numel(obj.y),...
                'x and y must be the same length');
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
        function val = get.B(obj)
            [k,l] = ndgrid(1:obj.Nwx,1:obj.Nwy);
            
            right_neighbors = mysub2ind([obj.Nwx,obj.Nwy],k(2:end,:),l(2:end,:));
            left_neighbors = mysub2ind([obj.Nwx,obj.Nwy],k(1:(end-1),:),l(1:(end-1),:));
            top_neighbors = mysub2ind([obj.Nwx,obj.Nwy],k(:,1:(end-1)),l(:,1:(end-1)));
            bottom_neighbors = mysub2ind([obj.Nwx,obj.Nwy],k(:,2:end),l(:,2:end));
            
            L = sparse(...
                [left_neighbors(:)',right_neighbors(:)',...
                top_neighbors(:)',bottom_neighbors(:)'],...
                [right_neighbors(:)',left_neighbors(:)',...
                bottom_neighbors(:)',top_neighbors(:)'],...
                1,obj.Nwx*obj.Nwy,obj.Nwx*obj.Nwy);
            L = diag(sparse(1./full(sum(L,2))))*L; % weight according to number of neighbors
            val = diag(sparse(ones([obj.Nwx*obj.Nwy,1]))) - L; % center minus average of neighbors
            
        end
        function val = get.Bx(obj)
            [k,l] = ndgrid(1:obj.Nwx,1:obj.Nwy);
            sz = [obj.Nwx,obj.Nwy];
            left  = mysub2ind(sz,k(1:(end-2),:),l(1:(end-2),:));
            self  = mysub2ind(sz,k(2:(end-1),:),l(2:(end-1),:));
            right = mysub2ind(sz,k(3:end,:),l(3:end,:));
            
            rows  = mysub2ind([obj.Nwx - 2,obj.Nwy],k(1:(end-2),:),l(1:(end-2),:));
            
            nrows = (obj.Nwx-2)*obj.Nwy;

            val = sparse(repmat(rows(:)',1,3),[left(:)',self(:)',right(:)'],...
                [-.5*ones(1,nrows),ones(1,nrows),-.5*ones(1,nrows)],...
                nrows,prod(sz));
        end
        function val = get.By(obj)
            [k,l] = ndgrid(1:obj.Nwx,1:obj.Nwy);
            sz = [obj.Nwx,obj.Nwy];
            left  = mysub2ind(sz,k(:,1:(end-2)),l(:,1:(end-2)));
            self  = mysub2ind(sz,k(:,2:(end-1)),l(:,2:(end-1)));
            right = mysub2ind(sz,k(:,3:end),l(:,3:end));
            
            rows  = mysub2ind([obj.Nwx,obj.Nwy-2],k(:,1:(end-2)),l(:,1:(end-2)));
            
            nrows = (obj.Nwx)*(obj.Nwy-2);

            val = sparse(repmat(rows(:)',1,3),[left(:)',self(:)',right(:)'],...
                [-.5*ones(1,nrows),ones(1,nrows),-.5*ones(1,nrows)],...
                nrows,prod(sz));
        end
    end
    
end

function ind = mysub2ind(sz,x,y)
% faster than the built-in method (no input error checking)
%   note: ONLY WORKS ON 2D ARRAYS
    ind = x + sz(1)*(y-1);
end