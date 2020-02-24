classdef Corrections
    %CORRECTIONS 

    methods(Static)
        function [x,y,z,d,cosw] = calcGeometry(Detector)
            [x,y,z] = Detector.xyz;
            d = sqrt(x.*x + y.*y + z.*z);
            cosw = (Detector.ed(1,3)*x + Detector.ed(2,3)*y + Detector.ed(3,3)*z)./d;
        end

        function value = solidAngle(qx,qy,d,cosw)
            % The solid angle in units of mrad^2
            % qx, qy, and d are in the same units (nominally mm)
            % so I multiply by 10^6 to go from rad^2 to mrad^2
            value = 1E6*qx*qy*cosw./(d.*d);
        end

        function value = polarization(s0,lambda,p,pn,x,y,z,d)
            pv1 = cross(s0,pn)*lambda;
            pv2 = cross(pv1,s0)*lambda;
            cos_phi1 = (pv1(1)*x + pv1(2)*y + pv1(3)*z)./d;
            cos_phi2 = (pv2(1)*x + pv2(2)*y + pv2(3)*z)./d;

            value = p*(1-cos_phi1.*cos_phi1) + ...
                (1-p)*(1-cos_phi2.*cos_phi2);
        end

        function value = efficiency(material,delta,lambda,cosw)
            kappa = geom.MaterialProperties.calcMu(material,lambda);
            value = 1 - exp(-kappa*delta./cosw);
        end

        function value = attenuation(material,lambda,d)
            mu = geom.MaterialProperties.calcMu(material,lambda);
            value = exp(-mu*d);
        end

        function [dx,dy] = parallax(XraySource,Detector)
            [x,y,~,~,cosw] = geom.Corrections.calcGeometry(Detector);
            qx = Detector.qx;
            qy = Detector.qy;
            lambda = XraySource.wavelength;
            material = Detector.sensorMaterial;
            kappa = geom.MaterialProperties.calcMu(material,lambda);
            delta = Detector.sensorThickness;

            sinw = sqrt(1-cosw.*cosw);
            rho = sqrt(x.*x + y.*y);

            v = -kappa*delta./cosw;
            dr = -(1/kappa)*sinw.*(1+v.*exp(v)./(1-exp(v)));

            dx = (x./rho).*dr/qx;
            dy = (y./rho).*dr/qy;
        end

        function [totalCorr,solidAngle,polarization,efficiency,...
                attenuation] = total(XraySource,Detector)
            [x,y,z,d,cosw] = geom.Corrections.calcGeometry(Detector);
            qx = Detector.qx;
            qy = Detector.qy;
            s0 = XraySource.wavevector;
            pn = XraySource.polarizationPlaneNormal;
            p = XraySource.fractionOfPolarization;
            lambda = XraySource.wavelength;

            solidAngle = geom.Corrections.solidAngle(qx,qy,d,cosw);
            polarization = geom.Corrections.polarization(s0,lambda,p,pn,x,y,z,d);

            efficiency = geom.Corrections.efficiency(...
                Detector.sensorMaterial,...
                Detector.sensorThickness,lambda,cosw);

            % assume air is the intervening material
            attenuation = geom.Corrections.attenuation('air',lambda,d);

            totalCorr = solidAngle.*polarization.*efficiency.*attenuation;
        end

        function value = d3s(XraySource,Detector,Spindle) 
            % swept reciprocal space volume in ang^-3
            
            qx = Detector.qx;
            qy = Detector.qy;
            s0 = XraySource.wavevector;
            lambda = XraySource.wavelength;
            m2 = Spindle.rotationAxis;
            dPhi = Spindle.oscillationRange*pi/180; % rotation range in radians

            [x,y,z,d,cosw] = geom.Corrections.calcGeometry(Detector);
            [sx,sy,sz] = geom.Scattering.xyz2s(x,y,z,XraySource);
            % solid angle in radians^2
            dOmega = 1E-6*geom.Corrections.solidAngle(qx,qy,d,cosw);

            % calculate inverse lorentz factor: Linv = |m2 . (s0 x s1)|/|s0||s1|
            s1x = sx+s0(1);
            s1y = sy+s0(2);
            s1z = sz+s0(3);
            % components of (s0 x s1)
            e1 = s1y*s0(3) - s1z*s0(2);
            e2 = s1z*s0(1) - s1x*s0(3);
            e3 = s1x*s0(2) - s1y*s0(1);
            % calculate Linv (using |s0|=|s1|=1/lambda)
            Linv = abs(m2(1)*e1 + m2(2)*e2 + m2(3)*e3)*lambda^2;
            value = Linv.*dOmega*dPhi/lambda^3;
        end


    end

end
