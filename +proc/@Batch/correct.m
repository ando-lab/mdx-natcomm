function [tf, errorMessage] = correct(varargin)

options = struct(...
    'workingDirectory','./',...
    'geometryIn','geom.mat',...
    'gridIn','filter.mat',...
    'bkgGeometryIn','geomBkg.mat',...
    'matOut','correct.mat',...
    'logOut','correct.log',...
    'binMode','coarse',...  %'coarse','fine'
    'readBackground',true,...
    'readImageHeaders',true,...
    'exposureTimeField','Exposure_time',...
    'parallel',false,...
    'binExcludedVoxels',true);

BatchProcess = proc.Batch('proc.Batch.correct',options,varargin{:});
options = BatchProcess.options; % post-processed options

if isempty(options.bkgGeometryIn)
    options.readBackground = false;
end

try % START MAIN SCRIPT

    BatchProcess.start();

    [InGeom,RefGeom,InIS,WedgeIS] = BatchProcess.readFromMatFile(...
        options.geometryIn,'InputGeometry','RefinedGeometry',...
        'InputImageSeries','WedgeImageSeries');

    [Grid,pixelMask] = BatchProcess.readFromMatFile(options.gridIn,...
        'Grid','pixelMask');

    InGeom.Detector.mask = InGeom.Detector.mask & pixelMask;

    if options.readBackground

        [BkgInputIS,BkgWedgeIS] = BatchProcess.readFromMatFile(...
            options.bkgGeometryIn,'InputImageSeries','WedgeImageSeries');
    else
        BkgInputIS = [];
        BkgWedgeIS = [];
    end

    [corr,corrExcl] = correctScript(1,options,...
        Grid,InGeom,RefGeom,InIS,WedgeIS,BkgInputIS,BkgWedgeIS);

    % save results
    if options.binExcludedVoxels
        BatchProcess.saveToMatFile(options.matOut,...
            'options',options,...
            'corr',corr,...
            'corrExcl',corrExcl);
    else
        BatchProcess.saveToMatFile(options.matOut,...
            'options',options,...
            'corr',corr);
    end

    BatchProcess.finish; % done!

catch errorMessage
    BatchProcess.stop(errorMessage);
end

tf = BatchProcess.hasCompleted;
errorMessage = BatchProcess.errorMessage;

end
