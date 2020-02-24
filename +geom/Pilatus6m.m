classdef Pilatus6m < geom.Detector
    methods
        function obj = Pilatus6m(varargin)
            obj@geom.Detector(varargin{:});
            a = defaultArgs();
            n = fieldnames(a);
            for j=1:length(n)
                nj = n{j};
                aj = a.(nj);
                if isempty(obj.(nj))
                    obj.(nj) = aj;
                else
                    if  ~isa(obj.(nj),class(aj)) || ...
                            (ischar(aj) && ~strcmpi(aj,obj.(nj))) ||...
                            (isnumeric(aj) && obj.(nj)~=aj)
                        warning(['argument ''%s'' incompatible with',...
                            'Pilatus6m detector and may cause problems'],nj);
                    end
                end
            end
        end
    end
    methods (Static=true)
        function vmask = isVirtualPixel()
            % Pilatus detectors have virtual pixels! see Kraft 2009
            vmask = false(2463,2527);
            [x,y] = virtualPixelCoords();
            vmask(x,:) = true;
            vmask(x+1,:) = true;
            vmask(x-1,:) = true;
            vmask(:,y) = true;
            vmask(:,y-1) = true;
            vmask(:,y+1) = true;
            vmask = vmask & ~geom.Pilatus6m.isGapPixel;
        end
        function vmask = isEdgePixel()
            % pixels on the boundary of detector panels
            [x,y] = panelBoundaryCoords();
            vmask = false(2463,2527);
            vmask(x,:) = true;
            vmask(:,y) = true;
            vmask = vmask & ~geom.Pilatus6m.isGapPixel;
        end
        function vmask = isGapPixel()
            [x,y] = panelGapCoords();
            vmask = false(2463,2527);
            vmask(x,:) = true;
            vmask(:,y) = true;
        end
        function ind = chipIndex(ix,iy,mode)
            if nargin==0
                mode = 'all';
            elseif nargin==2
                mode = 'pixels';
            end
            switch lower(mode)
                case 'all'
                    [chipx,chipy] = chipGridCoords;
                    ind = NaN*ones(2463,2527);
                    for yInd=1:size(chipx,1)
                        ind(chipx(yInd,1):chipx(yInd,2),chipy(yInd,1):chipy(yInd,2)) = yInd;
                    end
                case 'pixels'
                    [chipx,chipy] = chipCoords();            
                    xInd = NaN*zeros(size(ix));
                    yInd = NaN*zeros(size(ix));
                    for xInd=1:size(chipx,1)
                        xInd(ix(:) >= chipx(xInd,1) & ix(:) <= chipx(xInd,2)) = xInd;
                    end
                    for xInd=1:size(chipy,1)
                        yInd(iy(:) >= chipy(xInd,1) & iy(:) <= chipy(xInd,2)) = xInd;
                    end
                    ind = sub2ind([size(chipx,1),size(chipy,1)],xInd,yInd);
                case 'nearest'
                    [chipx,chipy] = chipCoords();
                    xedges = mean(reshape([1,reshape(chipx',1,[]),2463],2,[]),1);
                    yedges = mean(reshape([1,reshape(chipy',1,[]),2527],2,[]),1);
                    Nx = length(xedges) - 1;
                    Ny = length(yedges) - 1;
                    
                    [~,xInd] = histc(ix,[-Inf,xedges,Inf]);
                    [~,yInd] = histc(iy,[-Inf,yedges,Inf]);
                    xInd = xInd - 1;
                    yInd = yInd - 1;
                    
                    xInd(xInd<1) = 1;
                    yInd(yInd<1) = 1;
                    xInd(xInd>Nx) = Nx;
                    yInd(yInd>Ny) = Ny;
                    
                    ind = sub2ind([Nx,Ny],xInd,yInd);
                otherwise
                    error('did not recognize mode %s',mode);
            end
        end
 
        function Detector = binVirtualPixels(Detector)
            % in the 6m, virtual pixels lie in the following rows and
            % columns:
            [x,y] = virtualPixelCoords();
            
            % initialize xShift and yShift matrices if needed
            if isempty(Detector.xShift)
                Detector.xShift = zeros(2463,2527);
            elseif length(Detector.xShift) ==1
                Detector.xShift = Detector.xShift + zeros(2463,2527);
            end
            if isempty(Detector.yShift)
                Detector.yShift = zeros(2463,2527);
            elseif length(Detector.yShift) ==1
                Detector.yShift = Detector.yShift + zeros(2463,2527);
            end
            
            Detector.xShift(x+1,:) = Detector.xShift(x,:) -1;
            Detector.xShift(x-1,:) = Detector.xShift(x,:) +1;
            Detector.yShift(:,y+1) = Detector.yShift(:,y) -1;
            Detector.yShift(:,y-1) = Detector.yShift(:,y) +1;
            
            if ~isempty(Detector.mask)
                mx = true(2463,2527);
                my = true(2463,2527);
                mx(x,:) = Detector.mask(x,:) & ...
                    Detector.mask(x-1,:) & Detector.mask(x+1,:);
                my(y,:) = Detector.mask(y,:) & ...
                    Detector.mask(y-1,:) & Detector.mask(y+1,:);
                Detector.mask = Detector.mask & mx & my;
            end
        end
    end
end

function [x,y] = virtualPixelCoords()
y = (98:212:2430)';
[x1,x2] = ndgrid((0:4)*494,61:61:427);
x = x1(:) + x2(:);
end

function [x,y] = panelBoundaryCoords()
[y1,y2] = ndgrid((0:11)*212,[1,195]);
y = y1(:) + y2(:);
[x1,x2] = ndgrid((0:4)*494,[1,487]);
x = x1(:) + x2(:);
end

function [x,y] = panelGapCoords()
[y1,y2] = ndgrid((0:10)*212,196:212);
y = y1(:) + y2(:);
[x1,x2] = ndgrid((0:3)*494,488:494);
x = x1(:) + x2(:);
end

function [x,y] = chipCoords()
[y1,y2] = ndgrid([1,99],(0:11)*212);

[x1,x2] = ndgrid(1:61:428,(0:4)*494);

xstart = x1(:) + x2(:);
ystart = y1(:) + y2(:);

x = [xstart(:),xstart(:) + 59];
y = [ystart(:),ystart(:) + 96];
end

function [x,y] = chipGridCoords()
[y1,y2] = ndgrid([1,99],(0:11)*212);

[x1,x2] = ndgrid(1:61:428,(0:4)*494);

[xstart,ystart] = ndgrid(x1(:) + x2(:), y1(:) + y2(:));

x = [xstart(:),xstart(:) + 59];
y = [ystart(:),ystart(:) + 96];
end



function args = defaultArgs()
args = struct('name','Pilatus6m','nx',2463,'ny',2527,...
    'qx',0.172,'qy',0.172,'sensorMaterial','Silicon');
end