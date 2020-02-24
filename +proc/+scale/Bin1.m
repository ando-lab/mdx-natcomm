classdef Bin1
    %Bin1 - 1D binning model
    properties
        Nw % = 50           % nubmer of bins
        x % = linspace(0,1)'% points to bin
        xmin % = 0          % range of bin centers
        xmax % = 1
    end
    properties(Dependent = true)
        A  % binning operator: y(x) = A*w(x)
        w  % vector of bin centers
        dw % bin width
        B  % second derivative operator: w''(x) = B*w(x)
        Nx
    end
    properties(Dependent = true, Access = private)
        we % vector of bin edges
    end
    methods
        function obj = Bin1(varargin)
            if ~isempty(varargin)
                for j=1:2:length(varargin)
                    obj.(varargin{j}) = varargin{j+1};
                end
            end
        end
        function val = get.A(obj)
            [v,ix] = obj.map();
            val = sparse(v,ix,1,obj.Nx,obj.Nw);
        end
        function val = interp(obj,p)
            % fast calculation of A*p(:)
            [v,ix] = obj.map();
            val = zeros(obj.Nx,1);
            val(v) = p(ix);
        end
        function [v,ix] = map(obj)
            thisx = obj.x(:);
            [~,ix] = histc(thisx,[-Inf;obj.we;Inf]);
            ix = ix-1;
            isInside = ix > 0 & ix <= obj.Nw;
            v = (1:obj.Nx)';
            ix = ix(isInside);
            v = v(isInside);
        end
        function val = get.dw(obj)
            if obj.Nw==1
                val = 1;
            else
                val = (obj.xmax - obj.xmin)/(obj.Nw-1);
            end
        end
        function val = get.w(obj)
            val = linspace(obj.xmin,obj.xmax,obj.Nw)';
        end
        function val = get.we(obj)
            val = obj.w;
            val = [val - 0.5*obj.dw; val(end) + 0.5*obj.dw];
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

