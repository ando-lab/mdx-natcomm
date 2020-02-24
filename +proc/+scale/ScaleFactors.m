classdef ScaleFactors < handle
    %SCALEFACTORS is a wrapper for proc.scale.ScalingModel that makes the
    %fitting routine more efficient (i.e. it does not duplicate
    %calculations)
    
    properties(SetAccess = private)
        ScalingModel = proc.scale.ScalingModel.empty();
    end
    properties(SetAccess = private, GetAccess = private)
        aHasChanged = true;
        bHasChanged = true;
        cHasChanged = true;
        dHasChanged = true;
        aValCached = [];
        bValCached = [];
        cValCached = [];
        dValCached = [];
    end
    properties(SetAccess = immutable) % must be set in constructor
        sigma = [] % used by cCalc
        ix = []
        iy = []
        iz = []
        ip = []
        s  = []
    end
    properties(SetAccess = immutable, GetAccess = private)
        % pre-computed values of ScalingModel.aCalc(...) etc
        aCalcCached = {};
        bCalcCached = {};
        cCalcCached = {};
        dCalcCached = {};
    end
    properties(Dependent = true) 
        % wrappers for ScalingModel.a, b, c, d, sza, ...
        a
        b
        c
        d
        sza
        szb
        szc
        szd
        % wrappers for functions ScalingModel.aVal, b, c, ...
        aVal
        bVal
        cVal
        dVal
    end
    
    methods
        function obj = ScaleFactors(varargin)
            % give name,value pairs
            if ~isempty(varargin)
                for j=1:2:length(varargin)
                    obj.(varargin{j}) = varargin{j+1};
                end
            end
            obj.aCalcCached = cell(3,1); % A0, Bxy, Bz
            [obj.aCalcCached{:}] = obj.ScalingModel.aCalc(obj.ix,obj.iy,obj.iz);
            
            obj.bCalcCached = cell(2,1); % A, B
            [obj.bCalcCached{:}] = obj.ScalingModel.bCalc(obj.iz);
            
            obj.cCalcCached = cell(4,1); 
            [A0,Bx,By] = obj.ScalingModel.cCalc(obj.s,obj.iz);
            A = diag(sparse(1./obj.sigma))*A0;
            AA = A'*A; %<-- a slow step: moved here from ScaleFactors.cCalc
            obj.cCalcCached{1} = A;
            obj.cCalcCached{2} = Bx;
            obj.cCalcCached{3} = By;
            obj.cCalcCached{4} = AA;
            
            
            obj.dCalcCached = cell(1,1); % A0 = obj.Model.dCalc(obj.ip);
            [obj.dCalcCached{:}] = obj.ScalingModel.dCalc(obj.ip);
        end
        function set.a(obj,val)
            assert(isempty(val) || testSize(val,obj.sza),...
                'attempted to assign something of the wrong size');
            % test if val is different from obj.a
            if ~testEquality(val,obj.ScalingModel.a)
                obj.ScalingModel.a = val;
                obj.aHasChanged = true;
            end
        end
        function set.b(obj,val)
            assert(isempty(val) || testSize(val,obj.szb),...
                'attempted to assign something of the wrong size');
            % test if val is different from obj.b
            if ~testEquality(val,obj.ScalingModel.b)
                obj.ScalingModel.b = val;
                obj.bHasChanged = true;
            end
        end
        function set.c(obj,val)
            assert(isempty(val) || testSize(val,obj.szc),...
                'attempted to assign something of the wrong size');
            % test if val is different from obj.c
            if ~testEquality(val,obj.ScalingModel.c)
                obj.ScalingModel.c = val;
                obj.cHasChanged = true;
            end
        end
        function set.d(obj,val)
            assert(isempty(val) || testSize(val,obj.szd),...
                'attempted to assign something of the wrong size');
            % test if val is different from obj.d
            if ~testEquality(val,obj.ScalingModel.d)
                obj.ScalingModel.d = val;
                obj.dHasChanged = true;
            end
        end
        function val = get.aVal(obj)
            if obj.aHasChanged % update stored value
                obj.aValCached = obj.ScalingModel.aVal(obj.ix,obj.iy,obj.iz);
                obj.aHasChanged = false;
            end
            val = obj.aValCached;
        end
        function val = get.bVal(obj)
            if obj.bHasChanged % update stored value
                obj.bValCached = obj.ScalingModel.bVal(obj.iz);
                obj.bHasChanged = false;
            end
            val = obj.bValCached;
        end
        function val = get.cVal(obj)
            if obj.cHasChanged % update stored value
                obj.cValCached = obj.ScalingModel.cVal(obj.s,obj.iz);
                obj.cHasChanged = false;
            end
            val = obj.cValCached;
        end
        function val = get.dVal(obj)
            if obj.dHasChanged % update stored value
                obj.dValCached = obj.ScalingModel.dVal(obj.ip);
                obj.dHasChanged = false;
            end
            val = obj.dValCached;
        end
        function val = get.a(obj)
            val = obj.ScalingModel.a;
        end
        function val = get.b(obj)
            val = obj.ScalingModel.b;
        end
        function val = get.c(obj)
            val = obj.ScalingModel.c;
        end
        function val = get.d(obj)
            val = obj.ScalingModel.d;
        end
        function val = get.sza(obj)
            val = obj.ScalingModel.sza;
        end
        function val = get.szb(obj)
            val = obj.ScalingModel.szb;
        end
        function val = get.szc(obj)
            val = obj.ScalingModel.szc;
        end
        function val = get.szd(obj)
            val = obj.ScalingModel.szd;
        end
        function varargout = aCalc(obj)
            varargout = obj.aCalcCached;
        end
        function varargout = bCalc(obj)
            varargout = obj.bCalcCached;
        end
        function varargout = cCalc(obj)
            varargout = obj.cCalcCached;
        end
        function varargout = dCalc(obj)
            varargout = obj.dCalcCached;
        end
    end
    
end

function tf = testSize(val,sz)
    if length(sz)==1
        sz = [sz,1];
    end
    tf = ndims(val) == length(sz) && all(size(val)==sz);
end

function tf = testEquality(val1,val2)
tf = isa(val1,class(val2)) && ...
    numel(val1) == numel(val2) && ...
    all(val1(:)==val2(:));
end