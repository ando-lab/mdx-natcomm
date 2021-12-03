classdef MapTools < util.propertyValueConstructor
    %MAPTOOLS 
    
    properties
        Grid (1,1) latt.PeriodicGrid = latt.PeriodicGrid([1,1,1],[0,0,0],[1,1,1])
        Basis (1,1) latt.Basis = latt.Basis(1,1,1,90,90,90)
        SpaceGroup (1,1) symm.SpaceGroup = symm.SpaceGroup(1)
        type {mustBeMember(type,{'density','structurefactor','patterson','intensity'})} = 'density'
        isPeriodic (1,1) logical = false
    end
    properties(Dependent)
        Symmetry
    end
    
    methods
        function obj = MapTools(varargin)
            %MAPTOOLS 
            obj@util.propertyValueConstructor(varargin{:});
        end
        
        function S = get.Symmetry(obj)
            switch obj.type
                case 'density'
                    S = obj.SpaceGroup;
                case 'structurefactor'
                    warning('Laue group symmetry applies only to amplitudes... proceed with caution');
                    S = obj.SpaceGroup.LaueGroup;
                case 'patterson'
                    S = obj.SpaceGroup.PointGroup;
                case 'intensity'
                    S = obj.SpaceGroup.LaueGroup;
            end
        end
        
        function newobj = invert(obj)
            newobj = obj;
            newobj.Grid = obj.Grid.invert();
            newobj.Basis = obj.Basis.invert();
            switch obj.type
                case 'density'
                    newobj.type = 'structurefactor';
                case 'structurefactor'
                    newobj.type = 'density';
                case 'patterson'
                    newobj.type = 'intensity';
                case 'intensity'
                    newobj.type = 'patterson';
            end
        end
        
        function asumask = isASU(obj)
            S = obj.Symmetry;
            assert(isa(S,'symm.LaueGroup'));
            [h,k,l] = obj.Grid.grid();
            ndiv = obj.Grid.invert.P;
            h = round(ndiv(1)*h);
            k = round(ndiv(2)*k);
            l = round(ndiv(3)*l);
            asumask = S.testASU(h,k,l);
        end
        
        
        
        function varargout = table2array(obj,T)

            ncols = size(T,2) - 3; % number of columns to convert
            assert(ncols > 0);
            varargout = cell(1,ncols);
            
            [n1,n2,n3] = obj.Grid.frac2ind(T.(1),T.(2),T.(3),obj.isPeriodic);
            isincl = n1 >= 1 & n2 >= 1 & n3 >= 1 & ...
                     n1 <= obj.Grid.N(1) & n2 <= obj.Grid.N(2) & n3 <= obj.Grid.N(3);
            ind = sub2ind(obj.Grid.N,n1(isincl),n2(isincl),n3(isincl));
            T = T(isIncl,:);
            
            for n=1:ncols
                varargout{n} = NaN*ones(obj.Grid.N);
                varargout{n}(ind) = T.(3 + n);
            end
        end
        
        function T = array2table(obj,A,mask)
            if nargin < 3 || isempty(mask)
                mask = isnan(A);
            end
            [G0,m,cropfun] = obj.shrink2mask(obj.Grid,mask);
            A = cropfun(A);
            [f1,f2,f3] = G0.grid();
            f1 = f1(m);
            f2 = f2(m);
            f3 = f3(m);
            A = A(m);
            
            switch obj.type
                case 'density'
                    colnames = {'x','y','z','rho'};
                case 'structurefactor'
                    colnames = {'h','k','l','F'};
                case 'patterson'
                    colnames = {'x','y','z','P'};
                case 'intensity'
                    colnames = {'h','k','l','I'};
            end
            
            T = table(f1,f2,f3,A,'VariableNames',colnames);
        end
        
    end
    methods(Static)
        
        function [NewGrid,resizefun] = resize2bounds(G,f1min,f1max,f2min,f2max,f3min,f3max)
            [n1min,n2min,n3min] = G.frac2ind(f1min,f2min,f3min,false);
            [n1max,n2max,n3max] = G.frac2ind(f1max,f2max,f3max,false);
            
            N = 1 + [n1max,n2max,n3max] - [n1min,n2min,n3min];
            P = G.delta.*N;
            [o1,o2,o3] = G.ind2frac(n1min,n2min,n3min);
            
            NewGrid = latt.PeriodicGrid(N,[o1,o2,o3],P);
            
            ncmin = max([n1min,n2min,n3min],[1,1,1]);
            ncmax = min([n1max,n2max,n3max],G.N);
            
            npadpre = ncmin - [n1min,n2min,n3min];
            npadpost = [n1max,n2max,n3max] - ncmax;
            padval = NaN;
            
            resizefun = @(m) ...
                padarray(...
                  padarray(...
                    m(ncmin(1):ncmax(1),ncmin(2):ncmax(2),ncmin(3):ncmax(3)),...
                  npadpre,padval,'pre'),...
                npadpost,padval,'post');
            
        end
        
        function [NewGrid,newmask,cropfun] = shrink2mask(G,mask)
            
            mask3 = squeeze(any(mask,[1,2]));
            mask2 = squeeze(any(mask,[1,3]));
            mask1 = squeeze(any(mask,[2,3]));
            n1min = find(mask1,1,'first');
            n1max = find(mask1,1,'last');
            n2min = find(mask2,1,'first');
            n2max = find(mask2,1,'last');
            n3min = find(mask3,1,'first');
            n3max = find(mask3,1,'last');
            
            N = 1 + [n1max,n2max,n3max] - [n1min,n2min,n3min];
            P = G.delta.*N;
            [o1,o2,o3] = G.ind2frac(n1min,n2min,n3min);
            
            NewGrid = latt.PeriodicGrid(N,[o1,o2,o3],P);

            cropfun = @(m) m(n1min:n1max,n2min:n2max,n3min:n3max);
            
            newmask = cropfun(mask);
        end
        
    end
end

