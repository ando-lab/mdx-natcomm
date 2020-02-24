function [tf, errorMessage] = reintegrate(varargin)

options = struct(...
    'workingDirectory','./',...
    'geometryIn','geom.mat',...
    'readBackground',true,... % added June 13, 2019
    'bkgGeometryIn','geomBkg.mat',...
    'gridIn','grid.mat',...
    'scaleIn','scale.mat',...
    'matOut','reintegrate.mat',...
    'logOut','reintegrate.log',...
    'exposureTimeField','Exposure_time',...
    'minimumCounts',0,...
    'minimumPixels',1,...
    'parallel',false,...
    'binExcludedVoxels',false);

BatchProcess = proc.Batch('proc.Batch.reintegrate',options,varargin{:});
options = BatchProcess.options; % post-processed options

% set some defaults:
options.binMode = 'fine';
options.readImageHeaders = true;

try % START MAIN SCRIPT

    BatchProcess.start();

    [InGeom,RefGeom,InIS,WedgeIS] = BatchProcess.readFromMatFile(...
        options.geometryIn,'InputGeometry','RefinedGeometry',...
        'InputImageSeries','WedgeImageSeries');

    [Grid,pixelMask] = BatchProcess.readFromMatFile(options.gridIn,...
        'Grid','pixelMask');

    ScalingModel = BatchProcess.readFromMatFile(options.scaleIn,'ScalingModel');
    
    if options.readBackground
        [BkgInIS,BkgWedgeIS] = BatchProcess.readFromMatFile(...
            options.bkgGeometryIn,'InputImageSeries','WedgeImageSeries');
    else
        BkgInIS = [];
        BkgWedgeIS = [];
    end

    % uses options: 'binMode','binExcludedVoxels','parallel'
    [bin,binExcl] = integrateScript(1,options,Grid,pixelMask,...
        InGeom,RefGeom,InIS,WedgeIS);

    % uses options: 'binMode', 'binExcludedVoxels', 'parallel', 'readBackground'
    % 'readImageHeaders', exposureTimeField'
    [corr,corrExcl] = correctScript(1,options,Grid,...
        InGeom,RefGeom,InIS,WedgeIS,BkgInIS,BkgWedgeIS);

    AverageGeometry = averageGeometry(InGeom,RefGeom);

    % set dzmax to be 1/2 the oscillation width per wedge
    dzmax = (1+diff(RefGeom(1).Spindle.seriesFrameRange))/2;
    [bin,corr] = proc.Batch.repartitionScript(1,struct('dzmax',dzmax),Grid,bin,corr);
    
    if options.binExcludedVoxels
        [binExcl,corrExcl] = proc.Batch.repartitionScript(1,struct('dzmax',dzmax),Grid,binExcl,corrExcl);
    end
    
    diffuseTable = proc.Batch.exportScaledScript(1,options,...
        Grid,bin,corr,ScalingModel,AverageGeometry.Detector);
    
    if options.binExcludedVoxels
        
        braggTable = proc.Batch.exportScaledScript(1,options,...
            Grid,binExcl,corrExcl,ScalingModel,AverageGeometry.Detector);
        
        % find nearest neighbors to Bragg (central) voxel
        Grid2 = Grid(1).g0;
        [wI,fINear] = Grid(1).g0.hkl2index(0,0,0,0,0,0);
        [~,fINear] = Grid(1).g0.getRegion(1,wI,fINear);
        
        % get neighbor values from diffuseTable
        [~,fI] = Grid(1).g0.hkl2index(diffuseTable.h,diffuseTable.k,diffuseTable.l,...
            diffuseTable.dh,diffuseTable.dk,diffuseTable.dl);
        diffuseTableNear = diffuseTable(ismember(fI,fINear),:);
        % calculate the weighted average of neighbors (one point per h,k,l)
        diffuseTableNear = proc.scale.merge(diffuseTableNear,{'h','k','l'},'I','sigma');
        
        % assign background measurement to each Bragg reflection
        [lia,locb] = ismember(braggTable(:,{'h','k','l'}),diffuseTableNear(:,{'h','k','l'}),'rows');
        braggTable.bkgI = NaN*ones([size(braggTable,1),1]);
        braggTable.bkgSigma = NaN*ones([size(braggTable,1),1]);
        braggTable.bkgI(lia) = diffuseTableNear.I(locb(lia));
        braggTable.bkgSigma(lia) = diffuseTableNear.sigma(locb(lia));
        
        % apply Lorentz correction
        Vc = AverageGeometry.Crystal.UnitCell.vCell;
        lfac = (Vc*braggTable.d3s)*prod(Grid(1).ndiv);
        braggTable.intI = (braggTable.I - braggTable.bkgI).*lfac;
        braggTable.intSigma = sqrt(braggTable.sigma.^2 + braggTable.bkgSigma.^2).*lfac;

    end

    % save data
    if ~options.binExcludedVoxels
        BatchProcess.saveToMatFile(options.matOut,...
            'options',options,...
            'diffuseTable',diffuseTable,...
            'AverageGeometry',AverageGeometry);
    else
        BatchProcess.saveToMatFile(options.matOut,...
            'options',options,...
            'diffuseTable',diffuseTable,...
            'braggTable',braggTable,...
            'AverageGeometry',AverageGeometry);
    end

    BatchProcess.finish; % done!

catch errorMessage
    BatchProcess.stop(errorMessage);
end
tf = BatchProcess.hasCompleted;
errorMessage = BatchProcess.errorMessage;
end
