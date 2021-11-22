classdef ExpandTableToArray < util.propertyValueConstructor
    % ExpandTableToArray -
    properties
        hklcols = {'h','k','l','I'}; % first 3 must be miller indices (h,k,l)
        SpaceGroup = symm.SpaceGroup(1); % P1 is default
        symexpand = true;
        ndiv = [1,1,1]; % used if auto-generating P (default behavior)
        fid = 1; % where to print stats. 1 is command window. can also print to log file.
    end
    
    methods
        function obj = ExpandTableToArray(varargin)
            obj@util.propertyValueConstructor(varargin{:});
        end
        function [P,varargout] = run(obj,T,P)
            T = T(:,obj.hklcols);
            
            if nargin <3 || isempty(P)
                % auto-generate P
                hmin = min(T.(1));  hmax = max(T.(1));
                kmin = min(T.(2));  kmax = max(T.(2));
                lmin = min(T.(3));  lmax = max(T.(3));
                P = obj.hkl2grid(hmin,hmax,kmin,kmax,lmin,lmax);
            end
            
            ncols = size(T,2) - 3; % number of columns to convert
            varargout = cell(1,ncols);
            
            for n=1:ncols
                if obj.symexpand
                    varargout{n} = latt.hkl2grid(T.(1),T.(2),T.(3),T.(n+3),P,obj.SpaceGroup);
                else
                    varargout{n} = latt.hkl2grid(T.(1),T.(2),T.(3),T.(n+3),P);
                end
            end
            
            
        end
        
        function P = hkl2grid(obj,hmin,hmax,kmin,kmax,lmin,lmax)
            if obj.symexpand
                [hh,kk,ll] = ndgrid([hmin,hmax],[kmin,kmax],[lmin,lmax]);
                hkl0 = [hh(:),kk(:),ll(:)];
                Ops = obj.SpaceGroup.LaueGroup.generalPositions;
                hkl = cell(numel(Ops),1);
                for j = 1:numel(Ops)
                    hkl{j} = hkl0*Ops(j).r;
                end
                hkl = unique(cell2mat(hkl),'rows');
                
                hmin = min(hkl(:,1));  hmax = max(hkl(:,1));
                kmin = min(hkl(:,2));  kmax = max(hkl(:,2));
                lmin = min(hkl(:,3));  lmax = max(hkl(:,3));
            end
            N = 1 + obj.ndiv.*[hmax-hmin,kmax-kmin,lmax-lmin];
            N = round(N); % must be integers
            P = latt.PeriodicGrid(N,[hmin,kmin,lmin],N./obj.ndiv);
        end
        
        function R = design_grid(obj,Basis,smax)
            G0 = proc.script.GridDesigner(...
                'dgrid',0.5/smax,...
                'SpaceGroup',obj.SpaceGroup,...
                'Basis',Basis.invert).run();
            P0 = G0.PeriodicGrid;
            P0.N = P0.N.*obj.ndiv;
            P0.P = obj.ndiv;
            R = P0.invert;
        end
        
    end
end