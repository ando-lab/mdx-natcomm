function [Grid,pixelMask] = filterScript(fid,options,...
    InputGeometry,RefinedGeometry,InputImageSeries,WedgeImageSeries)

% VIRTUAL PIXEL CORRECTIONS (added October 14, 2017)
if isa(InputGeometry.Detector,'geom.Pilatus6m')
    fprintf('  masking virtual pixels (Pilatus 6m detector)\n');
    InputGeometry.Detector.mask = InputGeometry.Detector.mask & ...
        ~geom.Pilatus6m.isVirtualPixel;% & ~geom.Pilatus6m.isEdgePixel;
end

% MASK OUT BAD PIXELS (according to negative counts in first image)
fprintf(fid,'  reading first image (to update bad pixel mask)\n');
im1 = InputImageSeries.update(WedgeImageSeries(1)).datastore.preview;
InputGeometry.Detector.mask = InputGeometry.Detector.mask & im1>=0;

% for export
pixelMask = InputGeometry.Detector.mask;

% generate empty grid object
Grid = grid.Sub3dRef.empty();

numWedges = length(RefinedGeometry);

if options.parallel
    % start the parallel pool
    gcp();
end

try % shut down parallel tasks in case of error
    % loop over wedges
    for j=1:numWedges
        fprintf(fid,'  processing wedge %d of %d\n',j,numWedges);
        
        fprintf(fid,'    initializing the Integrater\n');
        Integrater = proc.Integrate('verbose',false);
        Integrater.WedgeGeometry = InputGeometry.update(RefinedGeometry(j));
        Integrater.ImageSeries = InputImageSeries.update(WedgeImageSeries(j));
        
        fprintf(fid,'    initializing the reference grid\n');
        [h,k,l] = Integrater.hklPredict(options.smax);
        Grid(j) = grid.Sub3dRef('ndiv',options.ndiv);
        Grid(j).ref = Grid(j).hkl2ref(h,k,l);
        
        if options.parallel
            fprintf(fid,'    Starting asynchronous filter process for this wedge\n');
            F(j) = parfeval(@run_filter,4,Integrater,Grid(j),options);
        else
            [tf,EM,Grid(j),wm] = run_filter(Integrater,Grid(j),options);
            if ~tf, rethrow(EM); end
        end
    end
    
    if options.parallel
        fprintf(fid,'  Waiting for processes to finish\n');
        
        for j=1:numWedges
            [completedIdx,tf,EM,thisGrid,wm] = fetchNext(F);
            
            if ~tf
                cancel(F);
                rethrow(EM); % pass the error message along
            end
            fprintf(fid,'    wedge %d process completed (%d remaining)\n',...
                completedIdx,numWedges-j);
            Grid(completedIdx) = thisGrid;
        end
    end
    
    % refine grid
    if options.refineGrid && numWedges > 1
        
        fprintf(fid,'  Refining the integration masks\n');
        RefinedGrid = grid.Sub3dRef.empty();
        
        RefinedGrid(1) = union(Grid(1),intersect(Grid(1),Grid(2)));
        for j=2:(numWedges-1)
            RefinedGrid(j) = union(Grid(j),...
                union(intersect(Grid(j),Grid(j-1)),...
                intersect(Grid(j),Grid(j-1))));
        end
        j = numWedges;
        RefinedGrid(j) = union(Grid(j),intersect(Grid(j),Grid(j-1)));
        
        Grid = RefinedGrid;
    end
catch EM
    if options.parallel && exist('F','var')
        cancel(F); % just in case
    end
    rethrow(EM);
end
end

function [tf,errorMessage,Grid,wm] = run_filter(Integrater,Grid,options)
tf = true;
errorMessage = [];
fid = 1;
try
    fprintf(fid,'    reading images and accumulating count histogram\n');
    [countHistogram,countOverflow] = ...
        Integrater.countHistogram(Grid,options.maxCount);
    
    fprintf(fid,'    initializing the Filter\n');
    occupiedIndex = find(sum(countHistogram,2));
    F = proc.Filter();
    F.countHistogram = full(countHistogram(occupiedIndex,:));
    
    fprintf(fid,'    designing moving window\n');
    F.neighborList = Grid.movingWindow(options.window,occupiedIndex);
    
    fprintf(fid,'    calculating moving weighted median\n');
    wm = F.medianFilter();
    
    fprintf(fid,'    calculating KL divergence of count histograms and Poisson distribution\n');
    [dkl,deltaDkl] = F.divergenceFromPoisson(wm);
    
    fprintf(fid,'    identifying outliers\n');
    [isOutlier,dklOut,meanCounts] = F.maskPoissonOutliers(wm,dkl,deltaDkl);
    
    fprintf(fid,'    generating integration mask\n');
    excludedIndices = union(occupiedIndex(isOutlier),find(countOverflow));
    Grid.voxelMask = Grid.ind2mask(excludedIndices);
catch errorMessage
    tf = false;
    wm = [];
end

end