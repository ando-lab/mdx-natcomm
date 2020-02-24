%% load half integer map

load calc/half_integer_map.mat I I1 I2 R sigma
%%
% calculate CC*, variance, and mean

SVR = proc.script.StatisticsVsRadius(...
        'edges',0.04:.01:.9,...
        'Basis',R.Basis,...
        'PeriodicGrid',R.PeriodicGrid);

sr = SVR.run(I); % calculate mean and variance
sr12 = SVR.run(I1,I2); % calculate mean, variance, covariance, and cc12

%% define the isotropic intensity component

cc12 = sr12.cc;
cc12(cc12<0) = 0; % make sure ccstar is real
ccstar = sqrt(2*cc12./(1+cc12));

Iiso = sr.av - sqrt(sr.var).*ccstar;

% figure;plot(sr.r,Iiso)

%% make a smooth version of this for interpolation purposes

IL = proc.scale.InterpLin1('Nw',101,'xmin',SVR.edges(1),'xmax',SVR.edges(end),'x',sr.r);

% fit the weighted mean
lambda = 1E1; % regularization parameter smoothness
pu = (IL.A'*IL.A + lambda*(IL.B'* IL.B))\(IL.A'*sr.av);
pu0 = (IL.A'*IL.A + lambda*(IL.B'* IL.B))\(IL.A'*Iiso);

%figure(1);clf;plot(sr.r,sr.av,'.');hold on;plot(sr.r,Iiso,'.');plot(IL.w,pu);plot(IL.w,pu0);
%% save the result

save calc/intensity_statistics.mat IL pu0 pu sr sr12

