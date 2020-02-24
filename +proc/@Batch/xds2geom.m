function [tf,errorMessage] = xds2geom(varargin)

options = struct(...
    'workingDirectory','./',...
    'xdsDir','',...
    'matOut','geom.mat',...
    'logOut','xds2geom.log');

BatchProcess = proc.Batch('proc.Batch.xds2geom',options,varargin{:});
options = BatchProcess.options; % post-processed options

try 

    BatchProcess.start();

    xdsDir = BatchProcess.checkInputDirectory(options.xdsDir);

    % START MAIN SCRIPT

    fprintf(1,'  reading XDS.INP\n');
    xdsInp = io.xds.Inp.read(xdsDir);

    fprintf(1,'  reading INTEGRATE.LP\n');
    xdsIntegrate = io.xds.Integrate.read(xdsDir);

    fprintf(1,'  creating geometry objects\n');
    InGeom = io.xds.Convert.inp2DiffractionExperiment(xdsInp);
    InIS = io.xds.Convert.inp2ImageSeries(xdsInp);
    RefGeom = io.xds.Convert.integrate2DiffractionExperiment(xdsIntegrate);
    WedgeIS = io.xds.Convert.integrate2ImageSeries(xdsIntegrate);

    % END MAIN SCRIPT

    BatchProcess.saveToMatFile(options.matOut,...
        'options',options,...
        'InputGeometry',InGeom,...
        'InputImageSeries',InIS,...
        'RefinedGeometry',RefGeom,...
        'WedgeImageSeries',WedgeIS);

    BatchProcess.finish; % done!

catch errorMessage
    BatchProcess.stop(errorMessage);
end

tf = BatchProcess.hasCompleted;
errorMessage = BatchProcess.errorMessage;

end
