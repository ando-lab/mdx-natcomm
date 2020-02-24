classdef Spindle < util.propertyValueConstructor
    properties
        rotationAxis = [] % 3-vector
        startingFrame = [] % integer
        startingAngle = [] % scalar, in degrees
        oscillationRange = [] % scalar, in degrees
    end
    
    methods
        function obj = Spindle(varargin)
            obj@util.propertyValueConstructor(varargin{:});
        end

        function value = rotationMatrix(obj,phi)
            % ROTATES RELATIVE TO A STARTING ANGLE OF ZERO

            phi0 = 0;
            dphi_radians = (phi-phi0)*pi/180;
            m2 = obj.rotationAxis;
            m2 = m2/sqrt(dot(m2,m2)); % normalize just in case

            % rotation matrix (rodrigues formula)
            value = cos(dphi_radians)*eye(3) + sin(dphi_radians)*...
                [0, -m2(3), m2(2); m2(3),0,-m2(1); -m2(2), m2(1),0] + ...
                (1-cos(dphi_radians))*(m2'*m2);
        end

        function value = frame(obj,phi)
           % return the frame number containing the angle phi
           x = (phi - obj.startingAngle)/obj.oscillationRange;
           y = floor(x);
           if (y-x+1) < (1E-6) % this is needed to fix roundoff errors
               y = y + 1;
           end
           value = y + obj.startingFrame;
        end

        function value = frameRange(obj,phiRange)
            % range of frames where at least half of each oscillation falls
            % within the phi range given
            x = (phiRange - obj.startingAngle)/obj.oscillationRange;
            value = [ceil(x(1) - 0.5 - 1E-6),floor(x(2) - 0.5 - 1E-6)] + ...
                obj.startingFrame;
            if value(2)<value(1)
                % there are no frames for which >half of one of them falls
                % in the range given
                value = [];
            end

        end

        function value = phiRange(obj,frameRange)
            if isempty(frameRange)
                value = []; return;
            elseif length(frameRange)==1
                nFrames = 1;
            elseif length(frameRange)==2
                nFrames = 1 + frameRange(2) - frameRange(1);
            end
            phiStart = (frameRange(1) - obj.startingFrame)*...
                obj.oscillationRange + obj.startingAngle;

            value = phiStart + [0,nFrames*obj.oscillationRange];
        end

        function value = phiMid(obj,frameRange)
            value = mean(obj.phiRange(frameRange));
        end

        function value = frame2phi(obj,frameList)
            phiStart = (frameList - obj.startingFrame)*...
                obj.oscillationRange + obj.startingAngle;

            value = phiStart + 0.5*obj.oscillationRange;
        end

    end
end
