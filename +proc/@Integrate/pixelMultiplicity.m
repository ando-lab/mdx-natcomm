function [m,mExcl] = pixelMultiplicity(obj,Grid,binMode,binExcludedVoxels)
% returns a sparse matrix, m, which records the number of
% times a given pixel on the detector appears in a voxel.
%
% m(i,j) is the number of times pixel(j) is mapped to voxel(i)
% during the rotation. The voxels are defined according to the
% input Grid.
%
% If mode is 'coarse, the function uses the coarser
% parent of Grid (integer h,k,l).

if nargin<=2 || isempty(binMode)
    binMode = 'fine';
end
if nargin<=3 || isempty(binExcludedVoxels)
    binExcludedVoxels = true;
end

switch lower(binMode)
    case 'fine'
        numRows = Grid.numRows*Grid.arraysize(2);
    case 'coarse'
        numRows = Grid.numRows;
    otherwise
        error('did not recognize mode argument');
end

nx = obj.WedgeGeometry.Detector.nx;
ny = obj.WedgeGeometry.Detector.ny;
numCols = nx*ny;

m = sparse(numRows,numCols);
mExcl = sparse(numRows,numCols);

[X,Y] = ndgrid(1:nx,1:ny);
pixelInd = sub2ind([nx,ny],X,Y);

for frameNumber = obj.ImageSeries.frames
    if obj.verbose
        fprintf(1,'  mapping frame number %d\n',frameNumber);
    end
    [hh,kk,ll] = obj.WedgeGeometry.frame2hkl(frameNumber);
    
    [rowIndex,fractionIndex,isExcl] = Grid.hkl2index(hh,kk,ll);
    isMasked = isnan(rowIndex(:));   % mask for pixelMask
    mask0 = ~(isMasked | isExcl(:)); % mask for count and pixel
    mask1 = ~isMasked & isExcl(:);   % mask for countExcl and pixelExcl
    
    switch lower(binMode)
        case 'fine'
            voxelInd = Grid.sub2ind(rowIndex(mask0),fractionIndex(mask0));
            m = m + sparse(voxelInd,pixelInd(mask0),1,numRows,numCols);
            if binExcludedVoxels
                voxelInd = Grid.sub2ind(rowIndex(mask1),fractionIndex(mask1));
                mExcl = mExcl + sparse(voxelInd,pixelInd(mask1),1,numRows,numCols);
            end
        case 'coarse'
            m = m + sparse(rowIndex(mask0),pixelInd(mask0),1,numRows,numCols);
            if binExcludedVoxels
                mExcl = mExcl + sparse(rowIndex(mask1),pixelInd(mask1),1,numRows,numCols);
            end
    end
end

end