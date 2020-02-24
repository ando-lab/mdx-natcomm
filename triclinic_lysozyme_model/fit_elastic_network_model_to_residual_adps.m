%% fit residual ADPs using an elastic network model

%% get experimental ADPs to fit
load model/atomic_model.mat Atoms P0 P SpaceGroup

isIncl = ~Atoms.isHet & ~strcmp(Atoms.element,'H'); % exclude heteroatoms and hydrogens
Atoms = Atoms(isIncl,:); % exclude heteroatoms

P0 = mat2cell(P0,3*ones(numel(isIncl),1),size(P0,2));
P0 = cell2mat(P0(isIncl));

P = mat2cell(P,3*ones(numel(isIncl),1),size(P,2));
P = cell2mat(P(isIncl));

% observed ADPs
Uobs = cell2mat(arrayfun(@(f) f.U,Atoms.fatom,'Uni',0));

%% calculate ADPs for the lattice dynamics model
load fit/fit_lattice_dynamics_model_to_reference.mat L
covMat = L.covarianceMatrix(L.normalModes(),[0,0,0]);
Ulatt = adpCalc(P0,covMat(1:6,1:6));

clear covMat L

%%
% projection operators to remove rigid-body component of motion
numOps = numel(SpaceGroup.generalPositions);

p1 = eye(size(P,2)) - full(P\P0)*full(P0\P); % get only the nonrigid part
p1 = kron(speye(numOps),p1);

clear numOps
%%
% set up ADP calculation from lattice model
load model/elastic_network_model.mat N M

% calculate which residue pairs the edges belong to
asu = N.Cell.AsymmetricUnit;
ix = cell2mat(arrayfun(@(g,n) n*ones(1,numel(g.x)),asu,1:numel(asu),'UniformOutput',false));
clear asu
g1 = ix(N.Edges.a1);
g2 = ix(N.Edges.a2);
clear ix
%%
% set up functions for least squares fitting

gammaFun = @(y) sqrt(y(g1).*y(g2)); % y is the coupling strength

Kfun = @(v1,v2) cellfun(@(k) p1'*k*p1,N.Hessian(v1,v2),'UniformOutput',false);

Pcell = kron(sparse(1,1,1,1,numel(SpaceGroup.generalPositions)),P);

weight = kron(ones(3,3),Atoms.occ);

fitfun = @(x) weight.*(adpFitFunction(Uobs,Pcell,M,Kfun,'parallel',gammaFun(x))-Ulatt);
%% test
tic;resid = fitfun(14*ones(129,1));toc % 1.3 seconds
sum(resid(:).^2) %  43.8744
%%

opts = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective',...
    'MaxFunctionEvaluations',10000,'Display','iter','UseParallel',true);
tic
xfit = lsqnonlin(fitfun,14*ones(129,1),zeros(129,1),Inf*ones(129,1),opts);
toc % 
%% console output
% (Feb 24, 2020)
%
%                                         Norm of      First-order 
%  Iteration  Func-count     f(x)          step          optimality
%      0        130         43.8744                         0.382
%  [...]
%     76      10010         20.4048    2.94017e-07       0.000723     
%% save result

N.springType = 'parallel';
N.springConstants = gammaFun(xfit);

K = cellfun(@(k) p1'*k*p1,N.Hessian(),'Uni',0);
L = nm.LatticeModel(M,K,[1,1,1]);

save fit/fit_elastic_network_model_to_residual_adps.mat L N
