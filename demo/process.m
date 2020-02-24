%% Process

%%
% check that mdx is on path
try
    geom.Detector();
catch
    error('mdx is not on the path! edit setup_environment.m and run');
end

%% Create input geometry objects
% hote: The normal workflow is to run XDS and import the geometry using
% proc.Batch.xds2geom() and proc.Batch.cbf2geom(). However, for this demo
% the simulated diffraction images are small and lack Bragg peaks, so
% XDS cannot be run beforehand. Running the input_geometry script
% creates the objects directly.

if ~isfolder('proc'), mkdir('proc'); end

input_geometry(); % creates proc/geom.mat and proc/geomBkg.mat

%% Integrate on a coarse grid (for scaling)

%% 
% *proc.Batch.filter*

[tf, errorMessage] = proc.Batch.filter(...
    'workingDirectory','proc',...
    'ndiv',[3,5,5],... % coarse grid divisions
    'maxCount',50,...  % histogram up to maxCount
    'parallel',false... % true requires parallel processing toolbox
    );  

if ~tf, rethrow(errorMessage);end
%% 
% *proc.Batch.integrate*

[tf,errorMessage] = proc.Batch.integrate(...
    'workingDirectory','proc',...
    'parallel',false... % requires parallel processing toolbox
    ); 

if ~tf, rethrow(errorMessage);end
%% 
% *proc.Batch.correct*

[tf,errorMessage] = proc.Batch.correct(...
    'workingDirectory','proc',...
    'parallel',false... % requires parallel processing toolbox
    );

if ~tf, rethrow(errorMessage);end
%% 
% *proc.Batch.export*
[tf,errorMessage] = proc.Batch.export(...
    'workingDirectory','proc');

if ~tf, rethrow(errorMessage);end
%% 
% *proc.Batch.combine*
[tf,errorMessage] = proc.Batch.combine(...
    'workingDirectory','proc',...
    'exportIn',{'export.mat'},...
    'mergeNeighbors',true,...
    'combineBragg',false,...
    'nwx',4,... % some default parameters for ScalingModel
    'nwy',4);

if ~tf, rethrow(errorMessage);end
%% 
% *proc.Batch.scale*

[tf,errorMessage] = proc.Batch.scale(...
    'workingDirectory','proc',...
    'scaleIn','combine.mat',...
    'bizMult',5,... % regularization parameter for smoothness of b
    'program',{'b',50,1E-4,'o',[],10,'a',20,1E-3}); 

if ~tf, rethrow(errorMessage);end

%% Plot of overall scale factor vs. frame number (b)

load proc/scale.mat ScalingModel

figure;plot(ScalingModel.Ib.w,ScalingModel.b);
xlabel('Frame Number');ylabel('Scale Factor (b)');

clear ScalingModel
%% 
% *proc.Batch.merge*

[tf,errorMessage] = proc.Batch.merge(...
    'workingDirectory','proc');

if ~tf, rethrow(errorMessage);end

%% Plot of coarse map intensity vs. resolution (1/d)
load proc/merge.mat hklMerge

figure;plot(hklMerge.s,hklMerge.I,'.');
xlabel('1/d (Å^{-1})');ylabel('Intensity');title('Coarse Map intensity vs. resolution');
clear hklMerge

%% Plot of coarse map intensity slice (l = 0)

load proc/merge.mat hklMerge

% create 3D array from hkl table
E2A = proc.script.ExpandTableToArray(...
    'hklcols',{'hasu','kasu','lasu','I'},...
    'SpaceGroup',symm.SpaceGroup(1),...
    'symexpand',true,...
    'ndiv',[1,1,1]);

[P,I] = E2A.run(hklMerge);
[~,~,nl0] = P.frac2ind(0,0,0);
[h,~,~] = P.ind2frac(1:P.N(1),0,0);
[~,k,~] = P.ind2frac(0,1:P.N(2),0);

figure;imagesc(h,k,I(:,:,nl0)');
xlabel('h');ylabel('k');title('Coarse Map Intensity (l = 0)');colormap hot;colorbar

clear h k nl0 P I E2A hklMerge
%% Integrate on a fine grid

%% 
% *proc.Batch.grid*
%
% make a fine grid

ndiv = [5,15,15]; % number of divisions for the fine grid

[tf,errorMessage] = proc.Batch.grid(...
    'workingDirectory','proc',...
    'matOut','grid.mat',...
    'smax',Inf,...
    'getSymmetryEquivalents',true,...
    'ndiv',ndiv,...
    'excludeBraggPosition',true,...
    'parallel',false... % true requires parallel computing toolbox
    );

if ~tf, rethrow(errorMessage);end

%% 
% *proc.Batch.reintegrate*

[tf,errorMessage] = proc.Batch.reintegrate(...
    'workingDirectory','proc',...
    'geometryIn','geom.mat',...
    'bkgGeometryIn','geomBkg.mat',...
    'gridIn','grid.mat',...
    'scaleIn','scale.mat',...
    'matOut','reintegrate.mat',...
    'logOut','reintegrate.log',...
    'minimumCounts',0,...
    'minimumPixels',1,...
    'parallel',false... % true requires parallel computing toolbox
    );

if ~tf, rethrow(errorMessage);end

%% 
% *MergeScaledDiffuse*

% get crystal info
load proc/export.mat AverageGeometry

M = proc.script.MergeScaledDiffuse(...
    'Grid',grid.Sub3d('ndiv',ndiv),...
    'Crystal',AverageGeometry.Crystal,...
    'nMin',2,...
    'workingDirectory','proc');

AverageGeometry

% run the script on reintegrate.mat
fn = M.mapToColumns({'reintegrate.mat'});

[hklMerge,isincl] = M.mergeColumns(fn);

% save table to mat file
save('proc/mergeFine.mat','hklMerge');

% split and merge for calculating cc1/2, cc*
rng(0,'twister'); % set random number generator for reproducibility
[hklMerge] = M.mergeRandomHalfSets(fn,isincl);

save('proc/mergeFineSplit.mat','hklMerge');

M.clearTmp(fn); % clear files from temporary directory

clear hklMerge isincl M

%% Plot of fine grid intensity vs. resolution (1/d)
load proc/mergeFine.mat hklMerge
load proc/export.mat AverageGeometry

h = hklMerge.h + hklMerge.dh/ndiv(1);
k = hklMerge.k + hklMerge.dk/ndiv(2);
l = hklMerge.l + hklMerge.dl/ndiv(3);

[sx,sy,sz] = AverageGeometry.hkl2s(h,k,l);
hklMerge.s = sqrt(sx.^2 + sy.^2 + sz.^2);

figure;plot(hklMerge.s,hklMerge.I,'.');
xlabel('1/d (Å^{-1})');ylabel('Intensity');title('Fine Map intensity vs. resolution');
clear hklMerge sx sy sz h k l AverageGeometry

%% Plot of fine grid intensity slice (l = 0)
load proc/mergeFine.mat hklMerge

hklMerge.hfrac = hklMerge.h + hklMerge.dh/ndiv(1);
hklMerge.kfrac = hklMerge.k + hklMerge.dk/ndiv(2);
hklMerge.lfrac = hklMerge.l + hklMerge.dl/ndiv(3);

% create 3D array from hkl table
E2A = proc.script.ExpandTableToArray(...
    'hklcols',{'hfrac','kfrac','lfrac','I'},...
    'SpaceGroup',symm.SpaceGroup(1),...
    'symexpand',true,...
    'ndiv',ndiv);

[P,I] = E2A.run(hklMerge);
[~,~,nl0] = P.frac2ind(0,0,0);
[h,~,~] = P.ind2frac(1:P.N(1),0,0);
[~,k,~] = P.ind2frac(0,1:P.N(2),0);

figure;imagesc(h,k,I(:,:,nl0)');colormap hot
xlabel('h');ylabel('k');title('Fine Map Intensity (l = 0)');colorbar

clear h k nl0 P I E2A hklMerge

%% Plot of CC* vs. resolution (1/d)

load proc/mergeFineSplit.mat hklMerge
load proc/export.mat AverageGeometry
C = AverageGeometry.Crystal;

hklMerge.hfrac = hklMerge.h + hklMerge.dh/ndiv(1);
hklMerge.kfrac = hklMerge.k + hklMerge.dk/ndiv(2);
hklMerge.lfrac = hklMerge.l + hklMerge.dl/ndiv(3);

ET2A = proc.script.ExpandTableToArray(...
        'hklcols',{'hfrac','kfrac','lfrac','I1','I2'},...
        'SpaceGroup',symm.SpaceGroup(1),...
        'symexpand',false,...
        'ndiv',ndiv);

[P,I1,I2] = ET2A.run(hklMerge);
    
B = latt.Basis(C.a,C.b,C.c,C.alpha,C.beta,C.gamma);
B = latt.OrientedBasis.realign(B,C.a_axis,C.b_axis,C.c_axis);
SVR = proc.script.StatisticsVsRadius(...
        'edges',0:.01:.6,... 
        'Basis',B.invert,...
        'PeriodicGrid',P);

sr12 = SVR.run(I1,I2);

% ccstar:
cc12 = sr12.cc;
ccstar = sqrt(2*cc12./(1+cc12));

figure;plot(sr12.r,ccstar,'-','DisplayName','CC*');hold on;
plot(sr12.r,cc12,'-','DisplayName','CC_{1/2}');
set(gca,'Ylim',[0,1]);
xlabel('1/d (Å^{-1})');ylabel('Correlation Coefficient (CC)');legend show

title('Fine map, correlation of random half-sets vs. resolution');

clear hklMerge C AverageGeometry P I1 I2 B SVR sr12 % save ccstar for next cell
