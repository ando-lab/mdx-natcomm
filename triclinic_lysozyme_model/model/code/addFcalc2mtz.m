function mtz = addFcalc2mtz(mtzFileName,Atoms,Basis,SpaceGroup,varargin)

%% process optional inputs

opts = struct('dgrid',0.5);

for j=1:2:(length(varargin)-1)
    fn = varargin{j};
    val = varargin{j+1};
    assert(ismember(fn,fieldnames(opts)),'invalid name/value pair');
    opts.(fn) = val;
end

%% load the MTZ file

% load the experimental data from ctruncate mtz file
[mtzArray,mtzCols] = io.mtz.read(mtzFileName);

[~,locb] = ismember({'H','K','L','F','SIGF'},mtzCols);
mtz = array2table(mtzArray(:,locb),'VariableNames',{'h','k','l','Fobs','sigma'});

%clear locb mtzArray mtzCols 

% if any are NaN, set the sigma to Inf and value to 0, so val/sigma = 0
mtz.sigma(isnan(mtz.Fobs)) = Inf;
mtz.Fobs(isinf(mtz.sigma)) = 0;

%% do the Fourier transforms

nx = round(Basis.a/opts.dgrid);
ny = round(Basis.b/opts.dgrid);
nz = round(Basis.c/opts.dgrid);

% make the lattice grid object
LG = latt.LatticeGrid(latt.PeriodicGrid([nx,ny,nz],[0,0,0],[1,1,1]),Basis);

[S,R,LG] = latt.molecularTransform(...
    Atoms.xsol,Atoms.ysol,Atoms.zsol,Atoms.fsol,[],'LatticeGrid',LG);
F = latt.molecularTransform(...
    Atoms.x,Atoms.y,Atoms.z,Atoms.fatom,Atoms.occ,'LatticeGrid',LG);

% expand to unit cell
[h,k,l] = R.PeriodicGrid.grid();
Ops = SpaceGroup.generalPositions;
hklp0 = [h(:),k(:),l(:),zeros(numel(h),1)];
Fcell = zeros(R.N);
Scell = zeros(R.N);
for n=1:numel(Ops)
    hklp = hklp0*Ops(n);
    phaseFactor = exp(2i*pi*hklp(:,4));
    [n1,n2,n3] = R.PeriodicGrid.frac2ind(hklp(:,1),hklp(:,2),hklp(:,3));
    ind = sub2ind(R.N,n1,n2,n3);
    Fcell(:) = Fcell(:) + phaseFactor.*F(ind);
    Scell(:) = Scell(:) + phaseFactor.*S(ind);
end

%% assign Fcell and Scell columns to the mtzTruncate table

[mtz.sx,mtz.sy,mtz.sz] = R.Basis.frac2lab(...
    mtz.h,mtz.k,mtz.l);

[h,k,l] = R.PeriodicGrid.grid();
[ia,locb] = ismember(round([h(:),k(:),l(:)]),table2array(mtz(:,{'h','k','l'})),'rows');
mtz.Fsolv = NaN*ones(size(mtz,1),1);
mtz.Fcalc = NaN*ones(size(mtz,1),1);
mtz.Fsolv(locb(ia)) = Scell(ia);
mtz.Fcalc(locb(ia)) = Fcell(ia);

clear h k l ia locb

end