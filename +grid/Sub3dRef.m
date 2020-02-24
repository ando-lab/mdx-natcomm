classdef Sub3dRef < grid.Sub3d
    properties
        ref % reference indices
        voxelMask          % a list of indices = [ i(:), j(:) ]
                           % such that voxel(i,j) is flagged
        negateMask = true  % if true, a voxel(i,j) is excluded when
                           % voxelMask has a row [i,j]
    end
    properties(Dependent = true)
        numRows
        g0
        maskedIndices
    end
    methods
        function obj = Sub3dRef(varargin)
            obj@grid.Sub3d(varargin{:});
        end
        function val = get.g0(obj)
            val = grid.Sub3d(obj);
        end
        function val = get.maskedIndices(obj)
            [i,j] = find(obj.voxelMask);
            val = obj.sub2ind(i,j);
        end

        % overload some methods
        function [rowIndex1,varargout] = getRegion(obj,d,rowIndex,varargin)
            % [rowIndex1,fractionIndex1] = getRegion(obj,d,rowIndex,fractionIndex)
            roi2 = cell(1,nargout-1);
            varargout = cell(1,nargout-1);
            [roi1,roi2{:}] = obj.g0.getRegion(d,obj.ref(rowIndex),varargin{:});

            [isincluded1,rowIndex1] = ismember(roi1,obj.ref);

            rowIndex1(~isincluded1) = NaN;

            if ~isempty(varargout)
                varargout{1} = roi2{1};
                varargout{1}(~isincluded1) = NaN;
                [rowIndex1,varargout{1}] = obj.applyMask(rowIndex1,varargout{1});
            end

        end

        function [NewGrid,ia,ib] = union(obj,Grid2)
            [NewGrid,ia,ib] = mergeGrids('union',obj,Grid2);
        end

        function [NewGrid,ia,ib] = intersect(obj,Grid2)
            [NewGrid,ia,ib] = mergeGrids('intersect',obj,Grid2);
        end

        function NewGrid = subdivide(obj,ndiv2)
            if nargin==1 || isempty(ndiv2)
                ndiv2 = [2,2,2];
            elseif length(ndiv2)==1
                ndiv2 = [1,1,1]*ndiv2;
            elseif length(ndiv2)==3
                % default... do nothing
            else
                error('ndiv2 must be of length 1 or 3');
            end
            NewGrid = obj;
            NewGrid.ndiv = NewGrid.ndiv.*ndiv2;
            NewGrid.voxelMask = [];

            [i,j] = find(obj.voxelMask);

            [h0,k0,l0] = index2hkl(obj,i,j);
            hStep = 1/NewGrid.ndiv(1); %take fractional steps
            kStep = 1/NewGrid.ndiv(2);
            lStep = 1/NewGrid.ndiv(3);

            d = (ndiv2 - 1)/2;

            % find neighbors in a cartesian grid
            [dh,dk,dl] = ndgrid(-d(1):d(1),-d(2):d(2),-d(3):d(3));
            dh = dh(:)'*hStep;
            dk = dk(:)'*kStep;
            dl = dl(:)'*lStep;

            [wholeIndex,fractionIndex] = NewGrid.hkl2index(...
                repmat(h0,1,size(dh,2))+repmat(dh,size(h0,1),1),...
                repmat(k0,1,size(dk,2))+repmat(dk,size(k0,1),1),...
                repmat(l0,1,size(dl,2))+repmat(dl,size(l0,1),1));
            NewGrid.voxelMask = NewGrid.sub2mask(wholeIndex,fractionIndex);

        end

        function wholeIndex = hkl2ref(obj,h,k,l)
            wholeIndex = obj.g0.hkl2index(h,k,l);
        end

        function ind = sub2ind(obj,rowIndex,fractionIndex)
            arraySize = [obj.numRows,obj.arraysize(2)];
            ind = sub2ind(arraySize,rowIndex,fractionIndex);
        end

        function [rowIndex,fractionIndex] = ind2sub(obj,ind)
            arraySize = [obj.numRows,obj.arraysize(2)];
            [rowIndex,fractionIndex] = ind2sub(arraySize,ind);
        end

        function mask = sub2mask(obj,i,j)
            mask = sparse(i,j,true,obj.numRows,obj.arraysize(2));
        end

        function mask = ind2mask(obj,ind)
            [i,j] = obj.ind2sub(ind);
            mask = obj.sub2mask(i,j);
        end

        % wrapper for getRegion with a simplified linear-indexed interface
        % returns the indices in the input array (ind) which are adjacent
        % (empty values are assigned NaN) -> this is kind of a weird way of
        % doing it, but it works for now. An ajacency matrix is the
        % ultimate goal anyway.

        function neighborList = movingWindow(obj,window,ind,allowedInd)
            if nargin==3 || isempty(allowedInd)
                allowedInd = ind;
            end

            [rowIndex,fractionIndex] = obj.ind2sub(ind(:));

            [rowIndex1,fractionIndex1] = obj.getRegion(...
                window,rowIndex,fractionIndex);

            regionIndex = obj.sub2ind(rowIndex1,fractionIndex1);

            [isincl,rj] = ismember(regionIndex,allowedInd(:));
            rj(~isincl) = NaN;

            ri = repmat((1:size(rj,1))',1,size(rj,2));
            ri = ri(~isnan(rj(:)));
            rj = rj(~isnan(rj(:)));

            neighborList = accumarray(ri,rj,[numel(ind),1], @(x) {x});
        end

        function indexArray = mapArray(obj,ind1ref,ind2ref,varargin)
        % mapArray - given a reference list of indices and target h,k,l
        % arrays (see hkl2index), return an array a where 
        % h(i),k(i),l(i) maps to index a(i) in the reference list.
        % a is NaN for any h,k,l which does not map to the reference
        % list.
            [ind1array,ind2array] = obj.hkl2index(varargin{:});
            indarray = obj.sub2ind(ind1array,ind2array);
            indref = obj.sub2ind(ind1ref,ind2ref);
            [isincl,indexArray] = ismember(indarray,indref);
            indexArray(~isincl) = NaN;
        end

        function varargout = index2hkl(obj,rowIndex,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = obj.g0.index2hkl(obj.ref(rowIndex(:)),varargin{:});
        end

        function [rowIndex,fractionIndex,isExcl] = hkl2index(obj,varargin)
            % examples:
            %
            % [rowIndex,fractionIndex] = hkl2index(h,k,l);
            %
            %    rowIndex and fractionIndex will be equal to NaN if the
            %    h,k,l are not part of the refernece grid, or if subvoxels
            %    are specifically excluded by the grid object.
            %
            % [rowIndex,fractionIndex,isExcl] = hkl2index(h,k,l);
            %
            %    rowIndex and fractionIndex are equal to NaN only if h,k,l
            %    are not part of the reference grid. isExcl is false if
            %    the subvoxel is excluded by the grid object
            [wholeIndex,fractionIndex] = obj.g0.hkl2index(varargin{:});
            [isincluded,rowIndex] = ismember(wholeIndex,obj.ref);
            rowIndex(~isincluded) = NaN;
            fractionIndex(~isincluded) = NaN;
            if nargout < 3
                [rowIndex,fractionIndex] = obj.applyMask(rowIndex,fractionIndex);
            else %nargout == 3
                isExcl = obj.isExcluded(rowIndex,fractionIndex);
            end
        end

        function [rowIndex,fractionIndex] = applyMask(obj,rowIndex,fractionIndex)
            isExcl = obj.isExcluded(rowIndex,fractionIndex);
            rowIndex(isExcl) = NaN;
            fractionIndex(isExcl) = NaN;
        end

        function isExcl = isExcluded(obj,rowIndex,fractionIndex)
            isExcl = false(size(rowIndex));
            if ~isempty(obj.voxelMask) && nnz(obj.voxelMask) ~= 0
                thisInd = obj.sub2ind(rowIndex,fractionIndex);
                isExcl = ismember(thisInd,obj.maskedIndices);
            end
            if obj.negateMask
                % do nothing
            else
                isExcl = ~isExcl;
            end
        end

        function val = get.numRows(obj)
            val = length(obj.ref); % or size(obj.ref,1)? or numel(obj.ref)?
        end

        function [countHist,countOverflows,mask0] = countHistogram(obj,maxCount,im,varargin)

            countHistogramRows = obj.numRows*obj.arraysize(2);
            countHistogramCols = maxCount + 1;

            % return empty histograms if called without arguments
            if nargin==2
                countHist = sparse([],[],[],countHistogramRows,countHistogramCols);
                countOverflows = sparse([],[],[],countHistogramRows,1);
                mask0 = [];
                return;
            end

            [rowIndex,fractionIndex] = obj.hkl2index(varargin{:});
            mask0 = ~isnan(rowIndex) & im >= 0;
            mask1 = mask0 & im > maxCount;  % overflows
            mask2 = mask0 & im <= maxCount; % pixels to accumulate
            ind = obj.sub2ind(rowIndex(mask2),fractionIndex(mask2));
            ind1 = obj.sub2ind(rowIndex(mask1),fractionIndex(mask1));
            countHist = sparse(ind,double(im(mask2))+1,1,...
                countHistogramRows,countHistogramCols);
            countOverflows = sparse(ind1,1,double(im(mask1)),...
                countHistogramRows,1);
        end

    end
end

function [NewGrid,ia,ib] = mergeGrids(mode,Grid1,Grid2)
assert(isa(Grid2,class(Grid1)) && ...
    Grid1.nbits == Grid2.nbits && ...
    all(Grid1.ndiv == Grid2.ndiv) && ...
    Grid1.negateMask == Grid2.negateMask,...
    'grid specifications must be identical for union');

NewGrid = Grid1;
switch lower(mode)
    case 'intersect'
        [NewGrid.ref,ia,ib] = intersect(Grid1.ref,Grid2.ref,'stable');
    case 'union'
        [NewGrid.ref,ia,ib] = union(Grid1.ref,Grid2.ref,'stable');
    otherwise
        error('did not recognize merging mode');
end

NewGrid.voxelMask = NewGrid.sub2mask([],[]);  
if isempty(Grid1.voxelMask)
    Grid1.voxelMask = Grid1.sub2mask([],[]);
end
if isempty(Grid2.voxelMask)
    Grid2.voxelMask = Grid2.sub2mask([],[]);
end

if ~isempty(Grid1.ref) && nnz(Grid1.ref)~=0
    [isMember,newLocation] = ismember(Grid1.ref,NewGrid.ref);
    NewGrid.voxelMask(newLocation(isMember),:) = ...
        NewGrid.voxelMask(newLocation(isMember),:) | ...
        Grid1.voxelMask(isMember,:);
end
if ~isempty(Grid2.ref) && nnz(Grid2.ref)~=0
    [isMember,newLocation] = ismember(Grid2.ref,NewGrid.ref);
    NewGrid.voxelMask(newLocation(isMember),:) = ...
        NewGrid.voxelMask(newLocation(isMember),:) | ...
        Grid2.voxelMask(isMember,:);
end

end
