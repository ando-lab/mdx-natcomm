classdef ScalingModel < util.propertyValueConstructor
    % ScalingModel - scaling model for diffuse scattering
    %
    % interpolated model with detector (d), absorption (a), overall (b),
    % and offset (c) factors
    properties
        a % absorption factor - detector location (ix by iy) by frame number (iz)
        b % overall factor - frame number (iz)
        c % constant offset - scattering vector |s| by frame number (iz)
        d % detector panel - index (ip)
        ixLim % detector pixel range (first-dim)
        iyLim % detector pixel range (second-dim)
        izLim % frame range
        ipLim % detector panel index range
        sLim  % range of momentum transfer
        sza % size of a array
        szb % size of b array
        szc % size of c array
        szd % size of d array
    end
    properties(Dependent = true)
        Ia % interpolator for a
        Ib % interpolator for b
        Ic % interpolator for c
        Id % interpolator for d
    end
    methods
        function obj = ScalingModel(varargin)
            obj@util.propertyValueConstructor(varargin{:});
        end
        function val = get.Ia(obj)
            sz = obj.sza;
            val = proc.scale.InterpLin3(...
                'xmin',obj.ixLim(1),'xmax',obj.ixLim(2),...
                'ymin',obj.iyLim(1),'ymax',obj.iyLim(2),...
                'zmin',obj.izLim(1),'zmax',obj.izLim(2),...
                'Nwx',sz(1),'Nwy',sz(2),'Nwz',sz(3));
        end
        function val = get.Ib(obj)
            sz = obj.szb;
            val = proc.scale.InterpLin1(...
                'xmin',obj.izLim(1),'xmax',obj.izLim(2),...
                'Nw',prod(sz));
        end
        function val = get.Ic(obj)
            sz = obj.szc;
            val = proc.scale.InterpLin2(...
                'xmin',obj.sLim(1),'xmax',obj.sLim(2),...
                'ymin',obj.izLim(1),'ymax',obj.izLim(2),...
                'Nwx',sz(1),'Nwy',sz(2));
        end
        function val = get.Id(obj)
            sz = obj.szd;
            val = proc.scale.Bin1(...
                'xmin',obj.ipLim(1),'xmax',obj.ipLim(2),...
                'Nw',prod(sz));
        end
        function val = aVal(obj,ix,iy,iz)
            if isempty(obj.a)
                val = ones(size(ix));
            else
                I = obj.Ia;
                I.x = ix(:);
                I.y = iy(:);
                I.z = iz(:);
                val = reshape(I.interp(obj.a(:)),size(ix));
            end
        end
        function val = bVal(obj,iz)
            if isempty(obj.b)
                val = ones(size(iz));
            else
                I = obj.Ib;
                I.x = iz(:);
                val = reshape(I.interp(obj.b(:)),size(iz));
            end
        end
        function val = cVal(obj,s,iz)
            if isempty(obj.c)
                val = zeros(size(s));
            else
                I = obj.Ic;
                I.x = s(:);
                I.y = iz(:);
                val = reshape(I.interp(obj.c(:)),size(s));
            end
        end
        function val = dVal(obj,ip)
            if isempty(obj.d)
                val = ones(size(ip));
            else
                I = obj.Id;
                I.x = ip(:);
                val = reshape(I.interp(obj.d(:)),size(ip));
            end
        end
        
        function [a,b,c,d] = getScales(obj,ix,iy,iz,ip,s,n)
            assert((length(obj)==1 & nargin==6) | nargin==7,...
                'if more than one scaling model, n argument needed');
            nBatch = length(obj);
            a = NaN*zeros(size(ix));
            b = NaN*zeros(size(ix));
            c = NaN*zeros(size(ix));
            d = NaN*zeros(size(ix));
            for j=1:nBatch
                if nBatch > 1
                    isBatch = n==j;
                else
                    isBatch = true(size(iz));
                end
                a(isBatch) = obj(j).aVal(ix(isBatch),iy(isBatch),iz(isBatch));
                b(isBatch) = obj(j).bVal(iz(isBatch));
                c(isBatch) = obj(j).cVal(s(isBatch),iz(isBatch));
                d(isBatch) = obj(j).dVal(ip(isBatch));
            end
        end
        
        function [A,Bxy,Bz] = aCalc(obj,ix,iy,iz)
            Ia = obj.Ia;
            Ia.x = ix;
            Ia.y = iy;
            Ia.z = iz;
            A = Ia.A;
            Bxy = Ia.Bxy;
            Bz = Ia.Bz;
        end
        function [A,B] = bCalc(obj,iz)
            Ib = obj.Ib;
            Ib.x = iz;
            A = Ib.A;
            B = Ib.B;
        end
        function [A,Bx,By] = cCalc(obj,s,iz)
            Ic = obj.Ic;
            Ic.x = s;
            Ic.y = iz;
            A = Ic.A;
            Bx = Ic.Bx;
            By = Ic.By;
        end
        function [A] = dCalc(obj,ip)
            Id = obj.Id;
            Id.x = ip;
            A = Id.A;
        end
    end
end
