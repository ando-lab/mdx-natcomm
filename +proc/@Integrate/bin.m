function [count,pixel,n,countExcl,pixelExcl,nExcl] = bin(obj,Grid,...
    binMode,arrayType,binExcludedVoxels)

if nargin<=2 || isempty(binMode)
    binMode = 'fine';
end
if nargin <=3 || isempty(arrayType)
    arrayType = 'full';
end
if nargin <=4 || isempty(binExcludedVoxels)
    binExcludedVoxels = true;
end

switch lower(binMode)
    case 'fine'
        numRows = Grid.numRows;
        numCols = Grid.arraysize(2);
    case 'coarse'
        numRows = Grid.numRows;
        numCols = 1;
    otherwise
        error('did not recognize bin mode');
end

switch lower(arrayType)
    case 'full'
        count = zeros(numRows,numCols);
        pixel = zeros(numRows,numCols);
        n = zeros(numRows,numCols);
        if binExcludedVoxels
            countExcl = zeros(numRows,numCols);
            pixelExcl = zeros(numRows,numCols);
            nExcl = zeros(numRows,numCols);
        end
    case 'sparse'
        count = sparse(numRows,numCols);
        pixel = sparse(numRows,numCols);
        n = sparse(numRows,numCols);
        if binExcludedVoxels
            countExcl = sparse(numRows,numCols);
            pixelExcl = sparse(numRows,numCols);
            nExcl = sparse(numRows,numCols);
        end
    otherwise
        error('did not recognize array type');
end

ds = obj.ImageSeries.datastore;
while hasdata(ds)
    [im,iminfo] = ds.read;
    if obj.verbose, disp(iminfo); end

    frameNumber = iminfo.Label;
    [hh,kk,ll] = obj.WedgeGeometry.frame2hkl(frameNumber);

    [rowIndex,fractionIndex,isExcl] = Grid.hkl2index(hh,kk,ll);
    isMasked = isnan(rowIndex(:));   % mask for pixelMask
    mask0 = ~(isMasked | isExcl(:)); % mask for count and pixel
    mask1 = ~isMasked & isExcl(:);   % mask for countExcl and pixelExcl

    ind0 = [rowIndex(mask0),fractionIndex(mask0)];
    ind1 = [rowIndex(mask1),fractionIndex(mask1)];

    switch lower(binMode)
        case 'fine'
            % do nothing
        case 'coarse'
            ind0(:,2) = 1;
            ind1(:,2) = 1;
    end

    switch lower(arrayType)
        case 'full'
            count = count + accumarray(ind0,...
                double(im(mask0)),[numRows,numCols]);
            thisnpx = accumarray(ind0,1,[numRows,numCols]);
            pixel = pixel + thisnpx;
            n = n + thisnpx*frameNumber;

            if binExcludedVoxels
                countExcl = countExcl + accumarray(ind1,...
                    double(im(mask1)),[numRows,numCols]);
                thisnpx = accumarray(ind1,1,[numRows,numCols]);
                pixelExcl = pixelExcl + thisnpx;
                nExcl = nExcl + thisnpx*frameNumber;
            end

        case 'sparse'
            count = count + sparse(ind0(:,1),ind0(:,2),...
                double(im(mask0)),numRows,numCols);
            thisnpx = sparse(ind0(:,1),ind0(:,2),...
                1,numRows,numCols);
            pixel = pixel + thisnpx;
            n = n + thisnpx*frameNumber;

            if binExcludedVoxels
                countExcl = countExcl + sparse(ind1(:,1),ind1(:,2),...
                    double(im(mask1)),numRows,numCols);
                thisnpx = sparse(ind1(:,1),ind1(:,2),...
                    1,numRows,numCols);
                pixelExcl = pixelExcl + thisnpx;
                nExcl = nExcl + thisnpx*frameNumber;
            end
    end

end
end
