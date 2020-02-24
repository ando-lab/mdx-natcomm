function [hklRef,refData] = mtzFindStrong(mtzFileName,numRef,resMin,resMax)

[mtzArray,mtzCols,mtzCryst] = io.mtz.read(mtzFileName);
[~,locb] = ismember({'H','K','L','F','SIGF'},mtzCols);
mtzTruncate = array2table(mtzArray(:,locb),'VariableNames',{'h','k','l','Fobs','sigma'});
[mtzTruncate.sx,mtzTruncate.sy,mtzTruncate.sz] = mtzCryst.hkl2s(mtzTruncate.h,mtzTruncate.k,mtzTruncate.l);

% if any are NaN, set the sigma to Inf and value to 0, so val/sigma = 0
mtzTruncate.sigma(isnan(mtzTruncate.Fobs)) = Inf;
mtzTruncate.Fobs(isinf(mtzTruncate.sigma)) = 0;

s = sqrt(mtzTruncate.sx.^2 + mtzTruncate.sy.^2 + mtzTruncate.sz.^2);

% choose reflections between 2.5 and 2.0A resolution
isIncl = 1./s >= resMin & 1./s <= resMax & ~isinf(mtzTruncate.sigma);

% get the most intense hkl
refData = mtzTruncate(isIncl,:);
[~,ixorder] = sort(refData.Fobs,'descend');
refData = refData(ixorder(1:numRef),:);

% make sure h,k,l are in the ASU according to scatterbrain
[h,k,l] = mtzCryst.hkl2asu(refData.h,refData.k,refData.l);

hklRef = [h,k,l];

end