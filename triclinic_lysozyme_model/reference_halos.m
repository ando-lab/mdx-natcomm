%% create a reference dataset for fitting

%% identify strong Bragg peaks in a given resolution bin

mtzFileName = 'model/6o2h_aimless_truncate.mtz';

% get the 400 most-intense Bragg reflections in this resolution bin
numRef = 400;

% choose reflections between 2.5 and 2.0A resolution
resMin = 2; 
resMax = 2.5;

hklRef = mtzFindStrong(mtzFileName,numRef,resMin,resMax);

%% load diffuse data around those peaks

[hklGrid,I,sigma] = loadBrillouinZones('triclinic_lysozyme_maps.h5',...
    '/maps/processed/variational', hklRef);

%% export

save('fit/reference_halos.mat','hklGrid','I','sigma');