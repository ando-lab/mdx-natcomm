function [F,R,LG,GA] = molecularTransform(x,y,z,GA,occupancy,varargin)

% can optionally pass a lattice grid, which will override the other options
% dgrid and osr
opts = struct('rmax',3,'dgrid',0.5,'osr',4,'LatticeGrid',[],'Uadd',[]);

if isempty(occupancy)
    occupancy = ones(size(x));
end

for j=1:2:(length(varargin)-1)
    fn = varargin{j};
    val = varargin{j+1};
    assert(ismember(fn,fieldnames(opts)),'invalid name/value pair');
    opts.(fn) = val;
end

if ~isempty(opts.LatticeGrid)
    LG = opts.LatticeGrid;
else
    
    a = opts.osr*range(x);
    b = opts.osr*range(y);
    c = opts.osr*range(z);
    
    % grid design (target grid spacing of dgrid)
    nx = round(a/opts.dgrid);
    ny = round(b/opts.dgrid);
    nz = round(c/opts.dgrid);
    
    % make the lattice grid object
    LG = latt.LatticeGrid(latt.PeriodicGrid([nx,ny,nz],[0,0,0],[1,1,1]),...
        latt.OrientedBasis(a,b,c,90,90,90));
end

if ~isempty(opts.Uadd)
    % add a B-factor
    %disp('adding B-factor');
    for j=1:numel(GA)
        GA(j) = GA(j).addU(eye(3)*opts.Uadd);
    end
    Tsharp = latt.GaussianDensitySum(1,-opts.Uadd);
else
    Tsharp = [];
end

% make a spherical kernel
[k1,k2,k3] = LG.sphkernel(opts.rmax);

% calculate electron density (splat algorithm)
rhoFun = @(X,Y,Z,G,occ) occ*electronDensity(G,X,Y,Z);

rho = LG.splat(k1,k2,k3,rhoFun,@plus,0,x,y,z,GA,occupancy);

[F,R] = LG.ifft(rho);

if ~isempty(Tsharp)
    [sx,sy,sz] = R.grid();
    F = F.*Tsharp.scatteringAmplitude(sx,sy,sz);
end

end