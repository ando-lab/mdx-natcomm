function [corr,corrExcl] = correctScript(fid,options,...
    Grid,InGeom,RefGeom,InIS,WedgeIS,BkgInIS,BkgWedgeIS)

numWedges = length(RefGeom);
assert(length(Grid)==numWedges,'number of wedges should equal the number of grids');

corr = cell(numWedges,1);
corrExcl = cell(numWedges,1);

if options.parallel
    % start the parallel pool
    % parpool();
    gcp();
end
try % shut down parallel processes if error
    % loop over wedges
    for j=1:numWedges
        fprintf(fid,'  processing wedge %d of %d\n',j,numWedges);
        
        fprintf(fid,'    initializing the Integrater\n');
        Integrater = proc.Integrate('verbose',false);
        Integrater.WedgeGeometry = InGeom.update(RefGeom(j));
        Integrater.ImageSeries = InIS.update(WedgeIS(j));
        
        if options.readBackground
            BkgIS = BkgInIS.update(BkgWedgeIS(j));
        else
            BkgIS = [];
        end
        
        if options.parallel
            fprintf(fid,'    Starting asynchronous integration process for this wedge\n');
            F(j) = parfeval(@calculate_corrections,4,Integrater,Grid(j),...
                options,BkgIS);
        else
            [tf,EM,c,cE] = calculate_corrections(Integrater,Grid(j),options,BkgIS);
            if ~tf, rethrow(EM); end
            corr{j} = c;
            corrExcl{j} = cE;
        end
    end
    
    if options.parallel
        fprintf(fid,'  Waiting for processes to finish\n');
        
        for i=1:numWedges
            [j,tf,EM,c,cE] = fetchNext(F);
            
            if ~tf
                cancel(F);
                rethrow(EM); % pass the error message along
            end
            fprintf(fid,'    wedge %d process completed (%d remaining)\n',...
                j,numWedges-i);
            corr{j} = c;
            corrExcl{j} = cE;
        end
    end
catch EM
    if options.parallel && exist('F','var')
        cancel(F); % just in case
    end
    rethrow(EM);
end
end

function [tf,EM,corr,corrExcl] = calculate_corrections(...
    Integrater,Grid,options,BkgIS)
tf = true;
EM = [];
fid = 1;

corr = table([],[],[],[],[],[],[],[],[],[],[],[],[],[],...
    'VariableNames',...
    {'solidAngle','polarization','efficiency','attenuation',...
    'd3s','ix','iy','sx','sy','sz','dt','bkg','bkgErr','bkgDt'});
corrExcl = corr;

if isempty(Grid.ref) % nothing to integrate - save time and return!
    fprintf(fid,'    no grid points to load... skipping this wedge\n');
    return;
end

try
    fprintf(fid,'    calculating image corrections\n');
    
    im = struct('solidAngle',[],'polarization',[],'efficiency',[],...
        'attenuation',[],'d3s',[],'ix',[],'iy',[],'sx',[],'sy',[],'sz',[],...
        'dt',[],'bkg',[],'bkgDt',[]);
    
    % geometric intensity corrections
    [~,im.solidAngle,im.polarization,im.efficiency,im.attenuation] = ...
        Integrater.WedgeGeometry.intensityCorrections();
    
    % swept volume for each pixel (per frame)
    im.d3s = Integrater.WedgeGeometry.d3s();
    
    % location on detector surface in pixels
    [im.ix,im.iy] = ndgrid(...
        1:Integrater.WedgeGeometry.Detector.nx,...
        1:Integrater.WedgeGeometry.Detector.ny);
    
    % wavevector transfer in the lab frame
    [im.sx,im.sy,im.sz] = Integrater.WedgeGeometry.s();
    
    if options.readImageHeaders
        fprintf(fid,'    reading exposure times from crystal image headers\n');
        imh = readHeader(Integrater.ImageSeries);
        im.dt = mean([imh.(options.exposureTimeField)]);
    else
        im.dt = NaN;
    end
    
    if options.readBackground
        fprintf(fid,'    reading background images\n');
        [bkgim,bkgh] = readSum(BkgIS);
        im.bkg = double(bkgim);
        im.bkgDt = sum([bkgh.(options.exposureTimeField)]);
    else
        im.bkg = sparse(size(im.solidAngle,1),size(im.solidAngle,2));
        im.bkgDt = Inf; % was 0 - changed to Inf to prevent NaN from divide by zero
    end
    
    fprintf(fid,'    calculating the pixel multiplicity\n');
    [m,mExcl] = Integrater.pixelMultiplicity(Grid,...
        options.binMode,options.binExcludedVoxels);
    % note: mExcl will be empty if binExcludedVoxels is false
    
    fprintf(fid,'    accumulating image corrections\n');
    corr = accumulateTable(m);
    
    if options.binExcludedVoxels
        corrExcl = accumulateTable(mExcl);
    end
    
catch EM
    tf = false;
end

    function cc = accumulateTable(mm)
        
        cc = table();
        npx = full(sum(mm,2));
        cc.solidAngle   = (mm*im.solidAngle(:))./npx;
        cc.polarization = (mm*im.polarization(:))./npx;
        cc.efficiency   = (mm*im.efficiency(:))./npx;
        cc.attenuation  = (mm*im.attenuation(:))./npx;
        cc.d3s          = (mm*im.d3s(:))./npx;
        cc.ix           = (mm*im.ix(:))./npx;
        cc.iy           = (mm*im.iy(:))./npx;
        cc.sx           = (mm*im.sx(:))./npx;
        cc.sy           = (mm*im.sy(:))./npx;
        cc.sz           = (mm*im.sz(:))./npx;
        cc.dt           = im.dt*ones(size(npx));
        cc.bkg          = (mm*im.bkg(:))./npx;
        cc.bkgErr       = sqrt(mm.^2*im.bkg(:))./npx;
        cc.bkgDt        = im.bkgDt*ones(size(npx));
        
    end

end