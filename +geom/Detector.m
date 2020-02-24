classdef Detector < util.propertyValueConstructor
    properties
        name % 'Pilatus6m'
        nx   % number of pixels in 'fast' direction
        ny   % number of pixels in 'slow' direction
        qx   % pixel pitch in mm, 'fast' direction
        qy   % pixel pitch in mm, 'slow' direction
        ed = eye(3) % 3x3 matrix: identity = default value for detector orientation
        f    % detector distance in mm
        orgx % detector origin, 'fast' direction, pixel units
        orgy % detector origin, 'slow' direction, pixel units
        xShift % position correction, 'fast' direction, pixel units
        yShift % position correction, 'slow' direction, pixel units
        sensorMaterial % char string (such as 'Silicon')
        sensorThickness % thickness in mm
        mask % nx-by-ny binary mask. not used by any methods at the moment
    end
    
    methods
        function obj = Detector(varargin)
            obj@util.propertyValueConstructor(varargin{:});
        end

        function [x,y] = xy(obj)
            % pixel coordinate in the detector plane
            [ix,iy] = ndgrid(1:obj.nx,1:obj.ny);
            if ~isempty(obj.mask) % apply mask if supplied
                ix(~obj.mask) = NaN;
                iy(~obj.mask) = NaN;
            end
            x = obj.qx*(obj.xShift + ix - obj.orgx);
            y = obj.qy*(obj.yShift + iy - obj.orgy);
        end

        function [x,y] = ij2xy(obj,ix,iy)
            x = NaN*zeros(size(ix));
            y = NaN*zeros(size(iy));

            isok = ix>=1 & iy>=1 & ix<=obj.nx & iy<=obj.ny;

            % note: this assumes that both xShift and yShift
            % have the same dimensions, which should be true but is not
            % checked
            if numel(obj.xShift)>1 && numel(obj.yShift)>1
                ixy = sub2ind([obj.nx,obj.ny],ix(isok),iy(isok));
            elseif numel(obj.xShift)==1 && numel(obj.yShift)==1
                ixy = 1;
            else
                ixy = [];
            end
            x(isok) = obj.qx*(obj.xShift(ixy) + ix(isok) - obj.orgx);
            y(isok) = obj.qy*(obj.yShift(ixy) + iy(isok) - obj.orgy);
        end

        function [x,y,z] = ij2xyz(obj,ix,iy)
            [x0,y0] = ij2xy(obj,ix,iy);
            [x,y,z] = xy2xyz(obj,x0,y0);
        end

        function [x,y] = xyEdges(obj)
            % coordinates for pixels at the detector edge
            % this is useful for simulating scattering extents
            ix = [1:(obj.nx-1), obj.nx*ones(1,obj.ny-1), ...
                obj.nx:-1:2, ones(1,obj.ny-1)];
            iy = [ones(1,obj.nx-1), 1:(obj.ny-1), ...
                obj.ny*ones(1,obj.nx-1), obj.ny:-1:2];
            [x,y] = ij2xy(obj,ix,iy);
        end

        function [x,y,z] = xyzEdges(obj)
            % pixel coordinate in the laboratory frame
            [x0,y0] = obj.xyEdges;
            [x,y,z] = xy2xyz(obj,x0,y0);
        end

        function [x,y,z] = xyz(obj)
            % pixel coordinate in the laboratory frame
            [x0,y0] = obj.xy;
            [x,y,z] = xy2xyz(obj,x0,y0);
        end

        function [x,y,z] = xy2xyz(obj,x0,y0)
            x = obj.ed(1,1)*x0 + obj.ed(1,2)*y0 + obj.f*obj.ed(1,3);
            y = obj.ed(2,1)*x0 + obj.ed(2,2)*y0 + obj.f*obj.ed(2,3);
            z = obj.ed(3,1)*x0 + obj.ed(3,2)*y0 + obj.f*obj.ed(3,3);
        end

    end

end
