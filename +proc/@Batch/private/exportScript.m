function [diffuseTable,braggTable] = exportScript(...
    fid,options,AverageGeometry,Grid,bin,binExcl,corr,corrExcl)

    numWedges = length(bin);

    Vc = AverageGeometry.Crystal.UnitCell.vCell; 
    
    numRows = arrayfun(@(g) g.numRows,Grid);
    numCounts = cellfun(@(v) numel(v.pixels),bin);
    isCoarse = all(numRows(:) == numCounts(:));

    diffuseTable = table(); % initialize empty tables
    braggTable = table();
    
    % calculate the diffuse and bragg table
    for j=1:numWedges
        
        if isCoarse
            [h,k,l] = Grid(j).index2hkl(1:Grid(j).numRows);
        else
            [r,c] = ndgrid(1:Grid(j).numRows,1:Grid(j).arraysize(2));
            [h,k,l,dh,dk,dl] = Grid(j).index2hkl(r(:),c(:));
        end
        wedge = j*ones(size(bin{j}.counts));
        
        totalCorr = corr{j}.solidAngle.*...
                    corr{j}.polarization.*...
                    corr{j}.efficiency.*...
                    corr{j}.attenuation;
                
        sweptVolume = corr{j}.d3s.*bin{j}.pixels;
        
        % below is a hack to get things working with fine grids and no mask
        % need to work out the proper way to calculate in all cases.
        if isempty(Grid(j).ref)
            gridVolume = 0; % [] gives an error?
        elseif isempty(Grid(j).voxelMask)
            if isCoarse
                gridVolume = (1/Vc);
            else
                gridVolume = (1/Vc)/prod(Grid(j).ndiv);
            end
        else 
            if isCoarse
                gridVolume = (1/Vc)*(1 - sum(Grid(j).voxelMask,2)/prod(Grid(j).ndiv));
            else
                gridVolume = (1/Vc)/prod(Grid(j).ndiv);
            end
        end
        
        intensity  = bin{j}.counts./(bin{j}.pixels.*totalCorr.*corr{j}.dt);
        sigma =sqrt(bin{j}.counts)./(bin{j}.pixels.*totalCorr.*corr{j}.dt);
        sx = corr{j}.sx;
        sy = corr{j}.sy;
        sz = corr{j}.sz;
        ix = corr{j}.ix;
        iy = corr{j}.iy;
        iz = bin{j}.iz;
        bkgIntensity = corr{j}.bkg./(totalCorr.*corr{j}.bkgDt);
        bkgSigma = corr{j}.bkgErr./(totalCorr.*corr{j}.bkgDt);
        fraction = sweptVolume./gridVolume;
        
        if isCoarse
            thisDiffuseTable = table(h,k,l,...
                wedge,intensity,sigma,sx,sy,sz,ix,iy,iz,...
                bkgIntensity,bkgSigma,fraction);
        else
            thisDiffuseTable = table(h,k,l,dh,dk,dl,...
                wedge,intensity,sigma,sx,sy,sz,ix,iy,iz,...
                bkgIntensity,bkgSigma,fraction);
        end
        
        isInclDiffuse = bin{j}.pixels>=options.minimumPixels & ...
            bin{j}.counts>=options.minimumCounts;
        diffuseTable = [diffuseTable; thisDiffuseTable(isInclDiffuse,:)];
        
        if ~options.binExcludedVoxels || ~isCoarse
            continue;
        end
        
        % calculate bragg table
        totalCorr = corrExcl{j}.solidAngle.*...
                    corrExcl{j}.polarization.*...
                    corrExcl{j}.efficiency.*...
                    corrExcl{j}.attenuation;
        avgIntensity = binExcl{j}.counts./(totalCorr.*binExcl{j}.pixels.*corrExcl{j}.dt);
        avgSigma = sqrt(binExcl{j}.counts)./(totalCorr.*binExcl{j}.pixels.*corrExcl{j}.dt);
        
        sweptVolume = corrExcl{j}.d3s.*binExcl{j}.pixels;
        gridVolume = (1/Vc)*(sum(Grid(j).voxelMask,2)/prod(Grid(j).ndiv));
        
        intensity  = avgIntensity.*(Vc*corrExcl{j}.d3s).*binExcl{j}.pixels;
        sigma = avgSigma.*(Vc*corrExcl{j}.d3s).*binExcl{j}.pixels;
        sx = corrExcl{j}.sx;
        sy = corrExcl{j}.sy;
        sz = corrExcl{j}.sz;
        ix = corrExcl{j}.ix;
        iy = corrExcl{j}.iy;
        iz = binExcl{j}.iz;
        bkgIntensity = thisDiffuseTable.intensity.*(Vc*corrExcl{j}.d3s).*binExcl{j}.pixels;
        bkgSigma = thisDiffuseTable.sigma.*(Vc*corrExcl{j}.d3s).*binExcl{j}.pixels;
        fraction = sweptVolume./gridVolume;
        
        isInclBragg = isInclDiffuse & ...
            binExcl{j}.pixels >= options.minimumPixels & ...
            binExcl{j}.counts >= options.minimumCounts;
        
        thisBraggTable = table(h,k,l,...
                wedge,intensity,sigma,sx,sy,sz,ix,iy,iz,...
                bkgIntensity,bkgSigma,fraction);
        
        braggTable = [braggTable; thisBraggTable(isInclBragg,:)];
        
    end
end
