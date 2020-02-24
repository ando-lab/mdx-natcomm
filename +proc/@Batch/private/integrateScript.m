function [bin,binExcl] = integrateScript(fid,options,Grid,pixelMask,...
    InputGeometry,RefinedGeometry,InputImageSeries,WedgeImageSeries)

numWedges = length(RefinedGeometry);
assert(length(Grid)==numWedges,'number of wedges should equal the number of grids');

bin = cell(numWedges,1);
binExcl = cell(numWedges,1);

% apply pixel mask
InputGeometry.Detector.mask = InputGeometry.Detector.mask & pixelMask;

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
        Integrater.WedgeGeometry = InputGeometry.update(RefinedGeometry(j));
        Integrater.ImageSeries = InputImageSeries.update(WedgeImageSeries(j));
        
        if options.parallel
            fprintf(fid,'    Starting asynchronous integration process for this wedge\n');
            F(j) = parfeval(@accumulate_counts,4,Integrater,Grid(j),options);
        else
            [tf,EM,bin{j},binExcl{j}] = accumulate_counts(Integrater,Grid(j),options);
            if ~tf, rethrow(EM); end
        end
        
    end
    
    if options.parallel
        fprintf(fid,'  Waiting for processes to finish\n');
        
        for i=1:numWedges
            [j,tf,EM,b,bE] = fetchNext(F);
            if ~tf
                cancel(F);
                rethrow(EM); % pass the error message along
            end
            fprintf(fid,'    wedge %d process completed (%d remaining)\n',...
                j,numWedges-i);
            bin{j} = b;
            binExcl{j} = bE;
        end
    end
catch EM
    if options.parallel && exist('F','var')
        cancel(F); % just in case
    end
    rethrow(EM);
end

end


function [tf,EM,bin,binExcl] = accumulate_counts(Integrater,Grid,options)
% wrapper function for proc.Integrate.bin
tf = true;
EM = [];
fid = 1;

bin = table([],[],[],'VariableNames',{'counts','pixels','iz'});
binExcl = table([],[],[],'VariableNames',{'counts','pixels','iz'});

if isempty(Grid.ref) % nothing to integrate - save time and return!
    fprintf(fid,'    no grid points to load... skipping this wedge\n');
    return;
end

try
    binopts = {options.binMode,'full',options.binExcludedVoxels};
    % arrayType is only allowed to be full in this version of the code
    fprintf(fid,'    reading images and accumulating counts\n');
    if options.binExcludedVoxels
        [counts,pixels,n,countsExcl,pixelsExcl,nExcl] = Integrater.bin(...
            Grid,binopts{:});
        binExcl = table(countsExcl(:),pixelsExcl(:),nExcl(:)./pixelsExcl(:),...
            'VariableNames',{'counts','pixels','iz'});
    else
        [counts,pixels,n] = Integrater.bin(Grid,binopts{:});
    end
    bin = table(counts(:),pixels(:),n(:)./pixels(:),...
        'VariableNames',{'counts','pixels','iz'});
    
catch EM
    tf = false;
end
end