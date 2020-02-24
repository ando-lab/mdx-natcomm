classdef InterpLin1
    %InterpLin1 - 1D linear interpolation model
    properties
        Nw % = 50           % default number of samples
        x % = linspace(0,1)'% points to interpolate at
        xmin % = 0          % range of interpolating parameters
        xmax % = 1
    end
    properties(Dependent = true)
        A  % interpolation operator: y(x) = A*w(x)
        w  % vector of samples
        dw
        B  % second derivative operator: w''(x) = B*w(x)
        Nx
    end
    methods
        function obj = InterpLin1(varargin)
            if ~isempty(varargin)
                for j=1:2:length(varargin)
                    obj.(varargin{j}) = varargin{j+1};
                end
            end
        end
        function val = get.A(obj)
            [v,ix,u] = obj.map();
            val = sparse([v;v],[ix;ix+1],[1-u;u],obj.Nx,obj.Nw);
        end
        function val = interp(obj,p)
            % fast calculation of val = A*p(:);
            [v,ix,u] = obj.map();
            val = zeros(obj.Nx,1);
            val(v) = p(ix).*(1-u) + p(ix+1).*u;
        end
        function [v,ix,u] = map(obj)
            thisx = obj.x(:);
            [~,ix] = histc(thisx,[-Inf;obj.w;Inf]);
            ix = ix-1;
            isInside = ix > 0 & ix < obj.Nw;
            v = (1:obj.Nx)';
            ix = ix(isInside);
            v = v(isInside);
            u = (thisx(isInside) - obj.w(ix))*(1/obj.dw);
        end
        function val = get.dw(obj)
            val = (obj.xmax - obj.xmin)/(obj.Nw-1);
        end
        function val = get.w(obj)
            val = linspace(obj.xmin,obj.xmax,obj.Nw)';
        end
        function val = get.Nx(obj)
            val = length(obj.x);
        end
        function val = get.B(obj)
            val = sparse(1:(obj.Nw-2),1:(obj.Nw-2),-.5,obj.Nw-2,obj.Nw) + ...
                sparse(1:(obj.Nw-2),2:(obj.Nw-1),1,obj.Nw-2,obj.Nw) + ...
                sparse(1:(obj.Nw-2),3:(obj.Nw),-.5,obj.Nw-2,obj.Nw);
        end
    end
    
end

