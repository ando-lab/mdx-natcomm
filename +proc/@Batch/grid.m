function [tf, errorMessage] = grid(varargin)

options = struct(...
    'workingDirectory','./',...
    'geometryIn','geom.mat',...
    'matOut','grid.mat',...
    'logOut','grid.log',...
    'ndiv',[3,3,3],...
    'smax',Inf,...
    'hLim',[-Inf,Inf],...
    'kLim',[-Inf,Inf],...
    'lLim',[-Inf,Inf],...
    'h',[],... % - list of specific integer h values
    'k',[],... % - list of specific integer k values (same lenght as h)
    'l',[],... % - (...)
    'getSymmetryEquivalents',true,... % if true, expand h,k,l to symmetry equivalents
    'excludeBraggPosition',false... % if true, excludes all voxels with integer miller indices
    );

BatchProcess = proc.Batch('proc.Batch.grid',options,varargin{:});
options = BatchProcess.options; % post-processed options

try % START MAIN SCRIPT

    BatchProcess.start();

    [InGeom,RefGeom,InIS,WedgeIS] = BatchProcess.readFromMatFile(...
        options.geometryIn,'InputGeometry','RefinedGeometry',...
        'InputImageSeries','WedgeImageSeries');

    % MASK OUT BAD PIXELS (according to negative counts in first image)
    fprintf(1,'  reading first image (to update bad pixel mask)\n');
    im1 = InIS.update(WedgeIS(1)).datastore.preview;
    InGeom.Detector.mask = InGeom.Detector.mask & im1>=0;

    % make a reference grid (GridRef) if required
    useReference = ~isempty(options.h);

    if useReference
      hklref = [options.h(:),options.k(:),options.l(:)];
      if options.getSymmetryEquivalents
        % expand list to all symmetry equivalents
        sg = InGeom.Crystal.SpaceGroup;
        gp = sg.PointGroup.LaueGroup.generalPositionsMat;
        a = cellfun(@(m) m*(hklref'),gp,'UniformOutput',false);
        hklref = unique([a{:}]','rows');
      end
      GridRef = grid.Sub3dRef('ndiv',options.ndiv);
      GridRef.ref = GridRef.hkl2ref(hklref(:,1),hklref(:,2),hklref(:,3));
    end

    % generate empty grid object
    Grid = grid.Sub3dRef.empty();

    numWedges = length(RefGeom);

    % loop over wedges
    for j=1:numWedges
        fprintf(1,'  processing wedge %d of %d\n',j,numWedges);

        fprintf(1,'    initializing the Integrater\n');
        Integrater = proc.Integrate('verbose',false);
        Integrater.WedgeGeometry = InGeom.update(RefGeom(j));
        Integrater.ImageSeries = InIS.update(WedgeIS(j));

        fprintf(1,'    initializing the reference grid\n');
        [h,k,l] = Integrater.hklPredict(options.smax);

        isIncl = ...
            h >= options.hLim(1) & h <= options.hLim(2) & ...
            k >= options.kLim(1) & k <= options.kLim(2) & ...
            l >= options.lLim(1) & l <= options.lLim(2);

        Grid(j) = grid.Sub3dRef('ndiv',options.ndiv);
        Grid(j).ref = Grid(j).hkl2ref(h(isIncl),k(isIncl),l(isIncl));

        % intersect with reference grid if required
        if useReference
          Grid(j) = intersect(Grid(j),GridRef);
        end

        if options.excludeBraggPosition
          % punch out integer h,k,l voxels
          [h,k,l] = Grid(j).index2hkl(1:Grid(j).numRows);
          [ind1,ind2] = Grid(j).hkl2index(h,k,l);
          Grid(j).voxelMask = Grid(j).sub2mask(ind1,ind2);
        end
    end

    % save data
    BatchProcess.saveToMatFile(options.matOut,...
            'options',options,...
            'Grid',Grid,...
            'pixelMask',InGeom.Detector.mask);

    BatchProcess.finish; % done!

catch errorMessage
    BatchProcess.stop(errorMessage);
end

tf = BatchProcess.hasCompleted;
errorMessage = BatchProcess.errorMessage;

end
