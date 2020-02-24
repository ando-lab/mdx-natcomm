function [diffuseTable,braggTable,ScalingModel] = combineScript(fid,options,batch)

nBatches = length(batch);
diffuseTable = table();
braggTable = table();

% some default parameters for scaling model
nws = options.nws; % number of s points to use
nwx = options.nwx;
nwy = options.nwy;


for j=1:nBatches
    fprintf(fid,'  pre-processing batch %d\n',j);
    
    batch(j).numWedges = length(unique(batch(j).diffuseTable.wedge));
    batch(j).frameRange = batch(j).AverageGeometry.Spindle.seriesFrameRange;
    
    if options.mergeNeighbors
        framesPerWedge = (1 + batch(j).frameRange(2) - batch(j).frameRange(1))/batch(j).numWedges;
        izMax = options.neighborMaxSep*framesPerWedge;
        
        fprintf(fid,'    merging equivalent observations within %d frames of each other\n',izMax);
        
        newDiffuseTable = mergeDown(batch(j).diffuseTable,'mean',izMax);
        
        fprintf(fid,'      %d observations merged in diffuse table, %d remaining\n',...
            size(batch(j).diffuseTable,1) - size(newDiffuseTable,1),size(newDiffuseTable,1));
        
        batch(j).diffuseTable = newDiffuseTable;
        
        if options.combineBragg
            newBraggTable = mergeDown(batch(j).braggTable,'sum',izMax);
            numObs = size(newBraggTable,1);
            fprintf(fid,'      %d observations merged in Bragg table, %d remaining\n',...
                size(batch(j).braggTable,1) - numObs,numObs);
            
            batch(j).braggTable = newBraggTable;
        end
    end
    
    if options.combineBragg && ~isempty(options.fractionRangeBragg)

        fprintf(fid,'    removing observations from Bragg table with fraction outside range %s\n',...
             mat2str(options.fractionRangeBragg));
         
        isInRange = ...
            batch(j).braggTable.fraction > options.fractionRangeBragg(1) &...
            batch(j).braggTable.fraction < options.fractionRangeBragg(2);
        
        batch(j).braggTable = batch(j).braggTable(isInRange,:);
        
        fprintf(fid,'      %d observations removed, %d remaining\n',...
                sum(~isInRange),sum(isInRange));
               
    end
end

for j=1:nBatches
    
    fprintf(fid,'  loading diffuse data from batch %d into table\n',j);
    
    diffuseTableOut = calculateColumns(...
        batch(j).diffuseTable,...
        batch(j).AverageGeometry);
    
    diffuseTableOut.n = diffuseTableOut.n*j;
    diffuseTable = [diffuseTable; diffuseTableOut];
end
clear diffuseTableOut % free up some memory?

if options.combineBragg
    for j=1:nBatches
        
        fprintf(fid,'  loading Bragg data from batch %d into table\n',j);
        
        braggTableOut = calculateColumns(...
            batch(j).braggTable,...
            batch(j).AverageGeometry);
        
        braggTableOut.n = braggTableOut.n*j;
        braggTable = [braggTable; braggTableOut];
    end
    clear braggTableOut % free up some memory?
end

% initializing scaling models
smax = max(diffuseTable.s);
ScalingModel = initScalingModel(batch,smax,nws,nwx,nwy);

end

function ScalingModel = initScalingModel(batch,smax,nws,nwx,nwy)

ScalingModel = proc.scale.ScalingModel.empty();
nBatches = length(batch);
for j=1:nBatches
    numWedges = batch(j).numWedges;
    frameRange = batch(j).frameRange;
    nx = batch(j).AverageGeometry.Detector.nx;
    ny = batch(j).AverageGeometry.Detector.ny;
    np = max(reshape(batch(j).AverageGeometry.Detector.chipIndex,1,[]));
    sza = [nwx,nwy,numWedges+1];
    szb = [numWedges+1,1];
    szc = [nws,numWedges+1];
    szd = [np,1];
    
    ScalingModel(j) = proc.scale.ScalingModel(...
        'ixLim',[-0.5,0.5] + [1,nx],...
        'iyLim',[-0.5,0.5] + [1,ny],...
        'izLim',[-0.5,0.5] + frameRange,...
        'ipLim',[1,np],...
        'sLim',[0,ceil(smax*10000.1)/10000],...
        'sza',[nwx,nwy,numWedges+1],...
        'szb',numWedges+1,...
        'szc',[nws,numWedges+1],...
        'szd',np);
end

end

function tableOut = calculateColumns(tableIn,Geom)
h = tableIn.h;
k = tableIn.k;
l = tableIn.l;
[hasu,kasu,lasu] = Geom.Crystal.hkl2asu(h,k,l);
s = sqrt(tableIn.sx.^2 + tableIn.sy.^2 + tableIn.sz.^2);
ix = tableIn.ix;
iy = tableIn.iy;
p = Geom.Detector.chipIndex(ix,iy,'nearest');
I = tableIn.intensity - tableIn.bkgIntensity;
sigma = sqrt(tableIn.sigma.^2 + tableIn.bkgSigma.^2);
iz = tableIn.iz;
n = ones([size(tableIn,1),1]);
fraction = tableIn.fraction;
tableOut =  table(h,k,l,hasu,kasu,lasu,I,sigma,s,p,n,ix,iy,iz,fraction);
end

function hklMerge = mergeDown(tableIn,mode,izMax)

[hkl,~,ib] = unique(tableIn(:,{'h','k','l'}),'rows');
isEquiv = accumarray(ib,tableIn.iz,[size(hkl,1),1],...
    @(iz) length(iz)>1 && range(iz)<izMax);

tableIn.uniqueIndex = (1:size(tableIn,1))';
tableIn.uniqueIndex(isEquiv(ib)) = 0;
[hklMerge,~,ib] = unique(tableIn(:,{'h','k','l','uniqueIndex'}),'rows');
hklMerge.uniqueIndex = [];

sz = [size(hklMerge,1),1];

w = 1./(tableIn.sigma.^2 + tableIn.bkgSigma.^2);
wav = accumarray(ib,w,sz);

hklMerge.fraction = accumarray(ib,tableIn.fraction,sz);
hklMerge.iz = accumarray(ib,tableIn.iz.*w,sz)./wav;
hklMerge.ix = accumarray(ib,tableIn.ix.*w,sz)./wav;
hklMerge.iy = accumarray(ib,tableIn.iy.*w,sz)./wav;
hklMerge.iz = accumarray(ib,tableIn.iz.*w,sz)./wav;
hklMerge.sx = accumarray(ib,tableIn.sx.*w,sz)./wav;
hklMerge.sy = accumarray(ib,tableIn.sy.*w,sz)./wav;
hklMerge.sz = accumarray(ib,tableIn.sz.*w,sz)./wav;

switch lower(mode)
    case {'mean'}
        hklMerge.intensity     = accumarray(ib,w.*tableIn.intensity,sz)./wav;
        hklMerge.bkgIntensity  = accumarray(ib,w.*tableIn.bkgIntensity,sz)./wav;
        hklMerge.bkgSigma = sqrt(accumarray(ib,w.*tableIn.bkgSigma.^2,sz)./wav);
        hklMerge.sigma    = sqrt(accumarray(ib,w.*tableIn.sigma.^2,sz)./wav);
    case {'sum'}
        hklMerge.intensity     = accumarray(ib,tableIn.intensity,sz);
        hklMerge.bkgIntensity  = accumarray(ib,tableIn.bkgIntensity,sz);
        hklMerge.bkgSigma = sqrt(accumarray(ib,tableIn.bkgSigma.^2,sz));
        hklMerge.sigma    = sqrt(accumarray(ib,tableIn.sigma.^2,sz));
    otherwise
        error('did not recognize mode argument');
end

end