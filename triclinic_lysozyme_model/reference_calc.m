%% one-phonon structure factor calculations for reference data points

%% calculate one-phonon structure factors for PDB model

load model/atomic_model.mat Atoms Bsol Ufit

% define an origin (minimizes phase error in interpolated structure
% factors)
xcom = mean(Atoms.x);
ycom = mean(Atoms.y);
zcom = mean(Atoms.z);

% calculate structure factors for protein atoms
dx = Atoms.x-xcom;
dy = Atoms.y-ycom;
dz = Atoms.z-zcom;

opts = {'rmax',3,'dgrid',.5,'osr',3};

[F,R,LG] = latt.molecularTransform(dx,dy,dz,Atoms.fatom,Atoms.occ,opts{:});

opts = {'rmax',3,'LatticeGrid',LG}; % calculate on the same underlying grid

Fx = latt.molecularTransform(dx,dy,dz,Atoms.fatom,Atoms.x.*Atoms.occ,opts{:});
Fy = latt.molecularTransform(dx,dy,dz,Atoms.fatom,Atoms.y.*Atoms.occ,opts{:});
Fz = latt.molecularTransform(dx,dy,dz,Atoms.fatom,Atoms.z.*Atoms.occ,opts{:});

% calculate structure factors for solvent model
Uadd = -eye(3)*Bsol/(8*pi^2) - Ufit;

opts = {'rmax',3,'LatticeGrid',LG,'Uadd',Uadd};

dx = Atoms.xsol-xcom;
dy = Atoms.ysol-ycom;
dz = Atoms.zsol-zcom;
occSol = ones(size(Atoms,1),1);

S = latt.molecularTransform(dx,dy,dz,Atoms.fsol,occSol,opts{:});
Sx = latt.molecularTransform(dx,dy,dz,Atoms.fsol,Atoms.xsol,opts{:});
Sy = latt.molecularTransform(dx,dy,dz,Atoms.fsol,Atoms.ysol,opts{:});
Sz = latt.molecularTransform(dx,dy,dz,Atoms.fsol,Atoms.zsol,opts{:});

clear dx dy dz Biso Ufit opts Uadd occSol
%% interpolate structure factors around reference peaks

% load reference grid
load fit/reference_halos.mat hklGrid

% load lattice information
% load pdb_model.mat Basis SpaceGroup
load model/lattice_dynamics_model.mat N
Basis = N.Cell.Basis;
UCOps = N.Cell.UnitCellOperators;
clear N

B = symm.AffineTransformation(Basis.orthogonalizationMatrix,[0;0;0]);

% shift back to com
Opcom = inv(B)*symm.AffineTransformation(eye(3),[xcom;ycom;zcom])*B; 
Ops = arrayfun(@(o) o*Opcom,UCOps);

sGrid = arrayfun(@(g) latt.LatticeGrid(g,Basis.invert),hklGrid);

Fgrid = latt.sfinterp(R,F + S,sGrid,Ops);
Fxgrid = latt.sfinterp(R,Fx + Sx,sGrid,Ops);
Fygrid = latt.sfinterp(R,Fy + Sy,sGrid,Ops);
Fzgrid = latt.sfinterp(R,Fz + Sz,sGrid,Ops);

clear B Opcom Ops xcom ycom zcom F Fx Fy Fz R S Sx Sy Sz

%% calculate one-phonon structure factors

Wk = onePhononStructureFactors(Basis,UCOps,hklGrid,Fgrid,Fxgrid,Fygrid,Fzgrid);
K1 = hklGrid(1).invert.invert;

%%
save fit/reference_calc.mat Wk K1 hklGrid