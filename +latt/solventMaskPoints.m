function [r1,r2,r3,LG] = solventMaskPoints(Basis,Operators,x,y,z,r,isIon,varargin)
%SOLVENTMASKPOINTS

opts = struct('dgrid',0.5,'vdwProb',1.2,'ionProb',0.8,'rShrink',0.8);

for j=1:2:(length(varargin)-1)
    fn = varargin{j};
    val = varargin{j+1};
    assert(ismember(fn,fieldnames(opts)),'invalid name/value pair');
    opts.(fn) = val;
end

% calculate coordinates after being acted upon by symmetry operators
G = struct('x',{},'y',{},'z',{});
B = symm.AffineTransformation(Basis.orthogonalizationMatrix,[0;0;0]);

if isempty(Operators)
    Operators = symm.SymmetryOperator(eye(3),[0;0;0]);
end

if isempty(isIon)
    isIon = false(size(r));
end

for n=1:numel(Operators)
    Op = Operators(n);
    xyzn = B*Op*inv(B)*[x(:)';y(:)';z(:)'];
    G(n).x = xyzn(1,:)';
    G(n).y = xyzn(2,:)';
    G(n).z = xyzn(3,:)';
end

nx = round(Basis.a/opts.dgrid);
ny = round(Basis.b/opts.dgrid);
nz = round(Basis.c/opts.dgrid);

% make the lattice grid object
LG = latt.LatticeGrid(latt.PeriodicGrid([nx,ny,nz],[0,0,0],[1,1,1]),Basis);

% set rProb = vdwr + prob, for each atom
rProb = (~isIon)*opts.vdwProb + (isIon)*opts.ionProb;

% make a spherical kernel
rMax = max(r+rProb) + max(sqrt(sum(LG.delta.^2,1))); % maximum radius to calculate

% define the distance function
distfun = @(xx,yy,zz,rr) sqrt(xx.*xx + yy.*yy + zz.*zz) - rr;

LG3 = LG;
LG3.PeriodicGrid.N = LG3.PeriodicGrid.N*3;
LG3.PeriodicGrid.P = LG3.PeriodicGrid.P*3;
[k1,k2,k3] = LG3.sphkernel(rMax);

distMapASU = Inf*ones(LG3.N);
distMapOther = Inf*ones(LG3.N);
for n=1:numel(G)
    distMap = LG3.splat(k1,k2,k3,distfun,@min,Inf,G(n).x,G(n).y,G(n).z,r+rProb);
    if n==1
        distMapASU = distMap;
    else
        distMapOther = min(distMapOther,distMap);
    end
end
distMap = min(distMapOther,distMapASU); % map of the unit cell

% generate supercell by unit cell shifts

distMapNeighbor = Inf*ones(LG3.N);

for n1=0:2
    for n2=0:2
        for n3=0:2
            shiftVec = LG3.N.*[n1,n2,n3]/3;
            if all([n1,n2,n3]==[0,0,0])
                cellMap = distMapOther;
            else
                cellMap = distMap;
            end
            cellMap = circshift(cellMap,shiftVec);
            distMapNeighbor = min(distMapNeighbor,cellMap);
        end
    end
end


clear cellMap shiftVec n1 n2 n3

% a voxel belongs to the central chain if it is closer to an atom of that
% chain than to any of the neighboring chain's atoms
isCen = distMapASU < distMapNeighbor;
% note: this assumption may break at special positions

% make the total solvent map
isAccessible = min(distMapNeighbor,distMapASU) > 0;

% design a kernel for dilating the map
[k1,k2,k3] = LG3.sphkernel(opts.rShrink);
n1 = max(k1); n2 = max(k2); n3 = max(k3);
nhood = false(n1*2 + 1,n2*2 + 1,n3*2 + 1);
ind = sub2ind(size(nhood),k1 + n1 + 1,k2 + n2 + 1,k3 + n3 + 1);
nhood(ind) = true;

% dilate the map

solventMask = imdilate(isAccessible,nhood);


% assign the solvent mask for the central chain
solventMask = solventMask | ~isCen;

% Ideally, the new solvent mask would have 1/8 as many points as
% the mask of the unit cell, but this is not always the case because of
% special positions and the number of grid points used. I really should fix
% this problem, but I think for now it is good enough.

% next, get a list of points corresponding to the "atom" representation of
% the solvent mask.
shiftVec = LG3.N.*[1,1,1]/3;
solventMaskShift = circshift(solventMask,shiftVec);
ix = find(~solventMaskShift);
[n1,n2,n3] = ind2sub(LG3.N,ix);
[x1,x2,x3] = LG3.PeriodicGrid.ind2frac(n1,n2,n3);
[r1,r2,r3] = LG3.Basis.frac2lab(x1-1,x2-1,x3-1); % shift back to cell [0,0,0]


end

