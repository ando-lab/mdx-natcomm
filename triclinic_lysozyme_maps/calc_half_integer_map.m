%% Regrid to half-integer map

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

% define a grid for interpolation (half-integer hkl)
Pgrid = latt.PeriodicGrid((2*hklmax + 1),-0.5*(2*hklmax + 1),2*hklmax + 1);

% define the kernel for interpolation
[k1,k2,k3] = latt.PeriodicGrid([5,5,5],[-2,-2,-2],[5,5,5]).grid();
kernel = [k1(:),k2(:),k3(:)];
clear k1 k2 k3 % dont need these any more
nmin = 75; % need 75/(5*5*5)~60 percent coverage to attempt interpolation

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

[I,sigma] = regrid(I,sigma,SpaceGroup,R,Rout,nmin,kernel);

%% do the same for the random half datasets

load('proc/mdx/mergeFineSplit.mat','hklMerge');

% convert to fractional Miller indices
hklMerge.h = hklMerge.h + hklMerge.dh/ndiv(1);
hklMerge.k = hklMerge.k + hklMerge.dk/ndiv(2);
hklMerge.l = hklMerge.l + hklMerge.dl/ndiv(3);

E2A.hklcols = {'h','k','l','I1','sigma1'};

[~,I1,sigma1] = E2A.run(hklMerge,P);

[I1,sigma1] = regrid(I1,sigma1,SpaceGroup,R,Rout,nmin,kernel);

E2A.hklcols = {'h','k','l','I2','sigma2'};

[~,I2,sigma2] = E2A.run(hklMerge,P);

[I2,sigma2] = regrid(I2,sigma2,SpaceGroup,R,Rout,nmin,kernel);

%% save the results

R = Rout;

save calc/half_integer_map.mat I sigma I1 sigma1 I2 sigma2 R