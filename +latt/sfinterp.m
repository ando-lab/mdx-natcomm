function Finterp = sfinterp(R,F,Rnew,SymOp)

minpad = 8;
interpMethod = 'spline';

%warning('off','MATLAB:griddedInterpolant:MeshgridEval3DWarnId');

if nargin < 4
    SymOp = symm.SymmetryOperator(eye(3),[0;0;0]);
end

numGrids = numel(Rnew);
numOps = numel(SymOp);

Finterp = cell(numGrids,numOps);

for j = 1:numGrids

thisR = Rnew(j);   
[h1,h2,h3] = thisR.PeriodicGrid.grid();

for n=1:numOps
Op = SymOp(n);

% apply unit cell operator
hkln = [h1(:),h2(:),h3(:),zeros(numel(h1),1)]*Op;
phaseFactor = exp(2i*pi*hkln(:,4));

% find points to interpolate
[k1,k2,k3] = thisR.Basis.frac2lab(hkln(:,1),hkln(:,2),hkln(:,3));
[h1i,h2i,h3i] = R.Basis.lab2frac(k1,k2,k3);

% make a grid for interpolation purposes
nx = ceil(range(h1i))+minpad;
ny = ceil(range(h2i))+minpad;
nz = ceil(range(h3i))+minpad;

ho1 = round(mean(h1i));
ho2 = round(mean(h2i));
ho3 = round(mean(h3i));

bzext = latt.PeriodicGrid([nx,ny,nz],[0,0,0],[1,1,1]).invert;
bzext.ori = bzext.ori + [ho1,ho2,ho3];

% make grid points for reference data and interpolation
[n1,n2,n3] = ndgrid(1:bzext.N(1),1:bzext.N(2),1:bzext.N(3));
[h1g,h2g,h3g] = bzext.ind2frac(n1,n2,n3); % hkl grid

% find indices of reference points
[n1,n2,n3] = R.PeriodicGrid.frac2ind(h1g,h2g,h3g);
ind = sub2ind(R.N,n1,n2,n3);

% perform the interpolation
Finterp{j,n} = reshape(phaseFactor.*interpn(h1g,h2g,h3g,F(ind),h1i,h2i,h3i,interpMethod),thisR.N);
end
end

%warning('on','MATLAB:griddedInterpolant:MeshgridEval3DWarnId');

end