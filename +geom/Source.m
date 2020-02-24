classdef Source < util.propertyValueConstructor
    properties
        wavelength = [] % in Angstroms
        direction = []  % 3-vector
        fractionOfPolarization = [] % (0.5)1 = (un)polarized
        polarizationPlaneNormal = [] % 3-vector
    end
    properties(Dependent)
        wavevector
    end
    methods
        function obj = Source(varargin)
            obj@util.propertyValueConstructor(varargin{:});
        end

        function value = get.wavevector(obj)
            if ~isempty(obj.wavelength) && length(obj.direction)==3
                d = obj.direction;
                d = d/sqrt(dot(d,d));
                value = d/obj.wavelength;
            else
                value = [];
            end
        end
    end
end
