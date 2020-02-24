classdef DiffractionExperiment < util.propertyValueConstructor
    properties
        Source
        Spindle
        Crystal
        Detector
    end
    methods
        function obj = DiffractionExperiment(varargin)
            obj@util.propertyValueConstructor(varargin{:});
        end
        function [h,k,l] = frame2hkl(obj,frame)
            [h,k,l] = obj.phi2hkl(obj.Spindle.phiMid(frame));
        end
        function [h,k,l] = phi2hkl(obj,phi)
            % wrapper for geom.Scattering.hkl
            [h,k,l] = geom.Scattering.hkl(obj.Source,obj.Detector,...
                obj.Crystal,obj.Spindle,phi);
        end
        function [sx,sy,sz] = hkl2s(obj,h,k,l)
            % wrapper for geom.Scattering.hkl2s0
            [sx,sy,sz] = geom.Scattering.hkl2s0(h,k,l,obj.Crystal);
        end
        function [sx,sy,sz] = s(obj)
            % wrapper for geom.Scattering.s
            [sx,sy,sz] = geom.Scattering.s(obj.Source,obj.Detector);
        end
        function [sx,sy,sz] = q(obj)
            % wrapper for geom.Scattering.s
            [sx,sy,sz] = geom.Scattering.q(obj.Source,obj.Detector);
        end
        function d3s = d3s(obj)
            % wrapper for geom.Corrections.d3s
            d3s = geom.Corrections.d3s(...
                obj.Source,...
                obj.Detector,...
                obj.Spindle);
        end
        function [totalCorr,solidAngle,polarization,efficiency,...
                attenuation] = intensityCorrections(obj)
            % wrapper for geom.Corrections.total
            [totalCorr,solidAngle,polarization,efficiency,attenuation] = ...
                geom.Corrections.total(obj.Source,obj.Detector);
        end
    end
end
