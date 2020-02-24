%% Regrid map from 13x11x11 to 7x7x7

%% set up calculation

% unit cell info
load('proc/mdx/unitCellInventory.mat','Crystal');
Basis = latt.Basis(Crystal.a,Crystal.b,Crystal.c,Crystal.alpha,Crystal.beta,Crystal.gamma);
SpaceGroup = symm.SpaceGroup('P1');
clear Crystal

% load diffuse data and expand to array
ndiv = [13,11,11];
hklmax = [22,27,30];
P = latt.PeriodicGrid((2*hklmax+1).*ndiv,-0.5*(2*hklmax + 1 - 1./ndiv),(2*hklmax+1)); % data grid

ngrid = [7,7,7];

% define a grid for interpolation (half-integer hkl)
Pgrid = latt.PeriodicGrid((2*hklmax + 1).*ngrid,[0,0,0],2*hklmax + 1).invert.invert;

% create reciprocal space grids
R = latt.LatticeGrid(P,Basis.invert);
Rout = latt.LatticeGrid(Pgrid,Basis.invert);

%% load data and interpolate

load('proc/mdx/mergeFine.mat','hklMerge');

% convert to fractional Miller indices
hklMerge.h = hklMerge.h + hklMerge.dh/ndiv(1);
hklMerge.k = hklMerge.k + hklMerge.dk/ndiv(2);
hklMerge.l = hklMerge.l + hklMerge.dl/ndiv(3);

% create 3D array from hkl table
E2A = proc.script.ExpandTableToArray(...
    'hklcols',{'h','k','l','I','sigma'},...
    'SpaceGroup',SpaceGroup,...
    'symexpand',true,...
    'ndiv',[]);

[~,I,sigma] = E2A.run(hklMerge,P);

[I,sigma] = regrid(I,sigma,SpaceGroup,R,Rout);

%% do the same for the random half datasets

load('proc/mdx/mergeFineSplit.mat','hklMerge');

% convert to fractional Miller indices
hklMerge.h = hklMerge.h + hklMerge.dh/ndiv(1);
hklMerge.k = hklMerge.k + hklMerge.dk/ndiv(2);
hklMerge.l = hklMerge.l + hklMerge.dl/ndiv(3);

E2A.hklcols = {'h','k','l','I1','sigma1'};

[~,I1,sigma1] = E2A.run(hklMerge,P);

[I1,sigma1] = regrid(I1,sigma1,SpaceGroup,R,Rout);

E2A.hklcols = {'h','k','l','I2','sigma2'};

[~,I2,sigma2] = E2A.run(hklMerge,P);

[I2,sigma2] = regrid(I2,sigma2,SpaceGroup,R,Rout);

%% save the results

R = Rout;

save calc/interpolated_7x7x7.mat I sigma I1 sigma1 I2 sigma2 R