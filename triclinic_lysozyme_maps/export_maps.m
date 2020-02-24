%% write maps (and other metadata) to h5 file

%% 1) lay out the file

h5filename = 'export/triclinic_lysozyme_maps.h5';
delete(h5filename); % delete the file to clear

% make data groups
h5addgroup(h5filename,'/maps');
h5addgroup(h5filename,'/maps/raw/');
h5addgroup(h5filename,'/maps/processed/');
h5addgroup(h5filename,'/crystal');

%% 2) write some metadata

% write some header info
h5writeatt(h5filename,'/','creation_date',datestr(now));
h5writeatt(h5filename,'/','sample','hen lysozyme in triclinic space group (P1)');
%h5writeatt(h5filename,'/','author','Steve P. Meisburger');
%h5writeatt(h5filename,'/','citation','Meisburger, S. P., Case, D. A. & Ando, N. Diffuse X-ray Scattering from Correlated Motions in a Protein Crystal. bioRxiv 805424 (2019) doi:10.1101/805424');

% write some crystal info
load('proc/mdx/unitCellInventory.mat','Crystal');

loc = '/crystal';

h5writeatt(h5filename,loc,'spaceGroupNumber',uint32(Crystal.spaceGroupNumber));
h5writeatt(h5filename,loc,'a',single(Crystal.a));
h5writeatt(h5filename,loc,'b',single(Crystal.b));
h5writeatt(h5filename,loc,'c',single(Crystal.c));
h5writeatt(h5filename,loc,'alpha',single(Crystal.alpha));
h5writeatt(h5filename,loc,'beta',single(Crystal.beta));
h5writeatt(h5filename,loc,'gamma',single(Crystal.gamma));

clear Crystal

%% 3) write total intensity map

export_total(h5filename,'/maps/raw/total',22,27,30);

%% 4) write variational intensity map

export_variational(h5filename,'/maps/processed/variational',22,27,30);

%% 5) write interpolated intensity map

export_interpolated_7x7x7(h5filename,'/maps/processed/interpolated_7x7x7',22,27,30);

%% 5) write half-integer map

export_half_integer(h5filename,'/maps/processed/half_integer');
