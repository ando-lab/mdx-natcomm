%% calculate variational component of the diffuse intensity

%% load hkl tables and convert columns to arrays
load proc/mdx/mergeFine.mat

ndiv = [13,11,11];

hklMerge = hklMerge(~isinf(hklMerge.sigma),:);

hklMerge.h = hklMerge.h + hklMerge.dh/ndiv(1);
hklMerge.k = hklMerge.k + hklMerge.dk/ndiv(2);
hklMerge.l = hklMerge.l + hklMerge.dl/ndiv(3);

E2R = proc.script.ExpandTableToArray(...
        'hklcols',{'h','k','l','I','sigma'},... % first 3 must be miller indices (h,k,l)
        'SpaceGroup',symm.SpaceGroup(1),... % P1
        'symexpand',true,...
        'ndiv',ndiv);
    
[P,I,sigma] = E2R.run(hklMerge);

%% fill in missing data points

isOK = ~isnan(I);

ker = zeros(3,3,3);
ker([1,3],2,2) = 1;
ker(2,[1,3],2) = 1;
ker(2,2,[1,3]) = 1;

A0 = I;
A0(~isOK) = 0;
n0 = convn(isOK,ker,'same');
A0 = convn(A0,ker,'same')./n0;
isnew = ~isOK & n0>0;
I(isnew) = A0(isnew);
sigma(isnew) = Inf; % signals that the value was interpolated

%% subtract background

load('proc/mdx/unitCellInventory.mat','Crystal');
Basis = latt.Basis(Crystal.a,Crystal.b,Crystal.c,Crystal.alpha,Crystal.beta,Crystal.gamma);

[sx,sy,sz] = latt.LatticeGrid(P,Basis.invert).grid();
s = sqrt(sx.^2 + sy.^2 + sz.^2);

load calc/intensity_statistics.mat IL pu0

% total
isIncl = s(:) <= IL.xmax & s(:) >= IL.xmin;

IL.x = s(isIncl);
Ibkg = NaN*ones(P.N);
Ibkg(isIncl) = IL.interp(pu0);
%%
I = I - Ibkg;

ismissing = isnan(I(:)) & isIncl;
I(ismissing) = 0;
sigma(ismissing) = Inf;

%% get table

[h,k,l] = P.grid();
isASU = symm.SpaceGroup(1).LaueGroup.testASU(round(h,3),round(k,3),round(l,3));
isASU = isASU & ~isnan(I); % get rid of data outside the resolution limits

[h,k,l,dh,dk,dl] = grid.Sub3d('ndiv',ndiv).hkl2hkl(h(isASU),k(isASU),l(isASU));
I = single(I(isASU));
sigma = single(sigma(isASU));
hklTable = table(h,k,l,dh,dk,dl,I,sigma,'VariableNames',{'h','k','l','dh','dk','dl','I','sigma'});

%% write table

save calc/variational_intensity.mat hklTable