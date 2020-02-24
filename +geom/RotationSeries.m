classdef RotationSeries < geom.Spindle
    properties
        seriesFrameRange
    end
    
    methods
        
        function obj = RotationSeries(varargin)
            obj@geom.Spindle(varargin{:});
        end
        
        function frames = getMatchingFrames(obj,TargetSeries)
            
            angleRange = obj.phiRange(obj.seriesFrameRange);

            targetFrameRange = TargetSeries.seriesFrameRange;
            
            targetAngleRange = TargetSeries.phiRange(targetFrameRange);
            targetSpan = targetAngleRange(2)-targetAngleRange(1);
            
            if angleRange(1)>=0 && angleRange(end)<=360
                
                % map target so that the start angle is in the range [0,360]
                targetAngleRange = mod(targetAngleRange(1),360) + [0,targetSpan];
                
                if targetAngleRange(2) > 360
                    
                    targetAngleRange1 = [targetAngleRange(1),360];
                    targetAngleRange2 = [0,targetAngleRange(2) - 360];
                    
                    frameRange1 = obj.frameRange(targetAngleRange1);
                    frameRange2 = obj.frameRange(targetAngleRange2);
                    
                    frames = [frameRange1(1):frameRange1(2),...
                        frameRange2(1):frameRange2(2)];
                    
                else
                    frameRange1 = obj.frameRange(targetAngleRange);
                    frames = frameRange1(1):frameRange1(2);
                end
            elseif targetAngleRange(1) >= angleRange(1) && ...
                    targetAngleRange(2) <= angleRange(2)
                frameRange1 = obj.frameRange(targetAngleRange);
                frames = frameRange1(1):frameRange1(2);
            else
                % TODO: implement special case (or make more general)
                error('have not implemented this case yet');
            end
            
            assert(all(frames>=obj.seriesFrameRange(1) & ...
                frames <= obj.seriesFrameRange(2)),'frames out of range?');
            
            
        end
        
    end
end
