function diffuseTable = exportScaledScript(...
    fid,options,Grid,bin,corr,ScalingModel,Detector)

% assumes a fine grid

minimumPixels = options.minimumPixels;
minimumCounts = options.minimumCounts;

diffuseTable = table([],[],[],[],[],[],[],[],[],[],...
    'VariableNames',{'h','k','l','dh','dk','dl','I','sigma','d3s','iz'});

for j=1:length(Grid)

thisBin = bin{j};
thisCorr = corr{j};
thisGrid = Grid(j);

isInclDiffuse = thisBin.pixels>=minimumPixels & ...
    thisBin.counts>=minimumCounts;

thisBin = thisBin(isInclDiffuse,:);
thisCorr = thisCorr(isInclDiffuse,:);

[r,c] = ndgrid(1:thisGrid.numRows,1:thisGrid.arraysize(2));
r = r(isInclDiffuse);
c = c(isInclDiffuse);

[h,k,l,dh,dk,dl] = thisGrid.index2hkl(r,c);

totalCorr = thisCorr.solidAngle.*thisCorr.polarization.*...
    thisCorr.efficiency.*thisCorr.attenuation;

d3s = thisCorr.d3s.*thisBin.pixels;

intensity  = thisBin.counts./(thisBin.pixels.*totalCorr.*thisCorr.dt);
sigma =sqrt(thisBin.counts)./(thisBin.pixels.*totalCorr.*thisCorr.dt);
bkgIntensity = thisCorr.bkg./(totalCorr.*thisCorr.bkgDt);
bkgSigma = thisCorr.bkgErr./(totalCorr.*thisCorr.bkgDt);

s = sqrt(thisCorr.sx.^2 + thisCorr.sy.^2 + thisCorr.sz.^2);
ix = thisCorr.ix;
iy = thisCorr.iy;
iz = thisBin.iz;
p = Detector.chipIndex(ix,iy,'nearest');

[a,b,c,d] = getScales(ScalingModel,ix,iy,iz,p,s);

scale = a.*b.*d;
offset = c./b;

intensity = intensity./scale;
sigma = sigma./scale;
bkgIntensity = bkgIntensity./scale + offset;
bkgSigma = bkgSigma./scale;

I = intensity - bkgIntensity;
sigma = sqrt(sigma.^2 + bkgSigma.^2);

diffuseTable = [diffuseTable; table(h,k,l,dh,dk,dl,I,sigma,d3s,iz)];

end

end