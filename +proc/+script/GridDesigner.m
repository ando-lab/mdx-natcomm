classdef GridDesigner < util.propertyValueConstructor
    %GRIDDESIGNER
    
    properties
        dgrid = 0.3
        SpaceGroup = symm.SpaceGroup.empty()
        Basis = latt.Basis(1,1,1,90,90,90)
        maxPrimeFactor = 7
    end
    
    methods
        function obj = GridDesigner(varargin)
            %GRIDDESIGNER
            obj@util.propertyValueConstructor(varargin{:});
        end
        
        function [LG] = run(obj)
            
            if isempty(obj.SpaceGroup)
                m = [1;1;1];
            else
                m = obj.find_grid_factors();
            end
            
            nx = next_good_dimension(obj,obj.Basis.a/obj.dgrid,m(1));
            ny = next_good_dimension(obj,obj.Basis.b/obj.dgrid,m(2));
            nz = next_good_dimension(obj,obj.Basis.c/obj.dgrid,m(3));
            
            % make the lattice grid object
            LG = latt.LatticeGrid(latt.PeriodicGrid([nx,ny,nz],[0,0,0],[1,1,1]),obj.Basis);
            
        end
        
        
        
        function m = find_grid_factors(obj)
            
            ops = obj.SpaceGroup.generalPositions;
            
            r = [1,1,1]';
            
            for j = 1:numel(ops)
                t = ops(j).t;
                isLess = t~=0 & t < r;
                r(isLess) = t(isLess);
            end
            
            m = round(1./r); % round just in case
            
        end
        
        
        function n = next_good_dimension(obj,nMin,requiredFactor)
            
            if isempty(obj.maxPrimeFactor)
                f = Inf;
            else
                assert(numel(factor(obj.maxPrimeFactor))==1,'maxPrimeFactor must be prime!');
                f = obj.maxPrimeFactor;
            end
            
            n = floor(nMin);
            while mod(n,requiredFactor) ~= 0 || max(factor(n)) > f
                n = n + 1;
            end
            
        end
        

    end
end

