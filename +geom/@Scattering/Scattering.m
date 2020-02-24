classdef Scattering
    %SCATTERING a collection of methods for scattering geometry calculations
    %  
    
    methods(Static) 
        
        [sx,sy,sz] = xyz2s(x,y,z,Source)

        [sx,sy,sz] = hkl2s0(h,k,l,Crystal)
        
        [sx,sy,sz] = hkl2s(h,k,l,Crystal,Spindle,phi)
        
        [h,k,l] = s2hkl(sx,sy,sz,Crystal,Spindle,phi)
        
        [h,k,l] = ewald2hkl(Source,Crystal,Spindle,phi,smax)
        
        m = hkl2mult(h,k,l,Crystal)
        
        [hnew,knew,lnew] = hkl2asu(h,k,l,Crystal)
        
        [hnew,knew,lnew,ind] = asu2hkl(h,k,l,Crystal)
        
        function [h,k,l] = hkl(Source,Detector,Crystal,Spindle,phi)
            [sx,sy,sz] = geom.Scattering.s(Source,Detector);
            [h,k,l] = geom.Scattering.s2hkl(sx,sy,sz,Crystal,Spindle,phi);
        end
        
        function [sx,sy,sz] = s(Source,Detector)
            % s = (s1-s0)
            [x,y,z] = Detector.xyz;
            [sx,sy,sz] = geom.Scattering.xyz2s(x,y,z,Source);
        end
        
        function [qx,qy,qz] = q(Source,Detector)
            % q = 2*pi/s
            [qx,qy,qz] = geom.Scattering.s(Source,Detector);
            qx = 2*pi*qx;
            qy = 2*pi*qy;
            qz = 2*pi*qz;
        end
        
        function value = sMag(Source,Detector)
            [sx,sy,sz] = geom.Scattering.s(Source,Detector);
            value = sqrt(sx.*sx+sy.*sy+sz.*sz);
        end
        
        function value = qMag(Source,Detector)
            sMag = geom.Scattering.sMag(Source,Detector);
            value = 2*pi*sMag;
        end
        
        function value = dSpacing(Source,Detector)
            % d = 1/|s| = 2*pi/|q|
            sMag = geom.Scattering.sMag(Source,Detector);
            value = 1./sMag;
        end
    end
    
end