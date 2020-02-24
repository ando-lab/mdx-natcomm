classdef Integrate < util.propertyValueConstructor
    properties
        WedgeGeometry
        ImageSeries
        verbose = true
    end
    methods
        function obj = Integrate(varargin)
            obj@util.propertyValueConstructor(varargin{:});
        end
        
        [h,k,l] = hklPredict(obj,smax)
        
        [countHist,countOverflow] = countHistogram(obj,Grid,maxCount)
        
        [m,mExcl] = pixelMultiplicity(obj,Grid,binMode,binExcludedVoxels)
        
        [count,pixel,n,countExcl,pixelExcl,nExcl] = bin(obj,...
            Grid,binMode,arrayType,binExcludedVoxels)     
    end
end