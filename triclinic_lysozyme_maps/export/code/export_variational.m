function export_variational(h5filename,loc,hmax,kmax,lmax)


% script to export 7x7x7 interpolated map data to h5 file

%h5filename = 'test.h5';
%loc = '/variational';

%% load variational intensity map

load calc/variational_intensity.mat hklTable

ndiv = [13,11,11];

nh = hmax*2 + 1;
nk = kmax*2 + 1;
nl = lmax + 1;

P = latt.PeriodicGrid([nh,nk,nl].*ndiv,[-hmax,-kmax,0]-floor(ndiv/2)./ndiv,[nh,nk,nl]);

hklTable.h = hklTable.h + hklTable.dh/ndiv(1);
hklTable.k = hklTable.k + hklTable.dk/ndiv(2);
hklTable.l = hklTable.l + hklTable.dl/ndiv(3);

E2R = proc.script.ExpandTableToArray(...
        'hklcols',{'h','k','l','I','sigma'},... % first 3 must be miller indices (h,k,l)
        'symexpand',false);

[~,I,sigma] = E2R.run(hklTable,P);

clear hklTable
%%
% Export to h5 file

h5addgroup(h5filename,loc);

h5writeatt(h5filename,loc,'Description','Estimated variational intensity component');
h5writeatt(h5filename,loc,'units','electron scattering per unit cell');
h5writeatt(h5filename,loc,'grid_coordinates','fractional Miller indices (h,k,l)');
h5writeatt(h5filename,loc,'grid_size',uint32(P.N));
h5writeatt(h5filename,loc,'grid_ori',P.ori);
h5writeatt(h5filename,loc,'grid_delta',P.delta);

arraysize = P.N;
opts = {'Datatype','single','FillValue',single(NaN),'ChunkSize',ndiv,'Deflate',1};

fprintf(1,'writing data to group %s in file %s\n',loc,h5filename);
tic
f_old = dir(h5filename);

h5create(h5filename,fullfile(loc,'I'),arraysize,opts{:});
h5create(h5filename,fullfile(loc,'sigma'),arraysize,opts{:});

h5writeatt(h5filename,fullfile(loc,'I'),'Description','variational intensity component');
h5writeatt(h5filename,fullfile(loc,'sigma'),'Description','standard error of I. Inf means missing data were either interpolated or filled with zeros.');


% write data
h5write(h5filename,fullfile(loc,'I'),single(I));
h5write(h5filename,fullfile(loc,'sigma'),single(sigma));

f_new = dir(h5filename);

fprintf(1,'\b... %.3f Mb written in %.1f seconds\n',(f_new.bytes - f_old.bytes)/2^20,toc)

end