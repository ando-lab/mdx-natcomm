function export_total(h5filename,loc,hmax,kmax,lmax)
% export total scattering map to h5 file

%h5filename = 'test.h5';
%loc = '/total';

%%

load proc/mdx/mergeFine.mat hklMerge

ndiv = [13,11,11];

hklMerge = hklMerge(~isinf(hklMerge.sigma),:);

nh = hmax*2 + 1;
nk = kmax*2 + 1;
nl = lmax + 1;

P = latt.PeriodicGrid([nh,nk,nl].*ndiv,[-hmax,-kmax,0]-floor(ndiv/2)./ndiv,[nh,nk,nl]);

hklMerge.h = hklMerge.h + hklMerge.dh/ndiv(1);
hklMerge.k = hklMerge.k + hklMerge.dk/ndiv(2);
hklMerge.l = hklMerge.l + hklMerge.dl/ndiv(3);

E2R = proc.script.ExpandTableToArray(...
        'hklcols',{'h','k','l','I','sigma'},... % first 3 must be miller indices (h,k,l)
        'symexpand',false);

[~,I,sigma] = E2R.run(hklMerge,P);

clear hklMerge nh nk nl hmin hmax kmin kmax lmin lmax

% get random half datasets also
load proc/mdx/mergeFineSplit.mat hklMerge

hklMerge = hklMerge(~isinf(hklMerge.sigma1),:);

hklMerge.h = hklMerge.h + hklMerge.dh/ndiv(1);
hklMerge.k = hklMerge.k + hklMerge.dk/ndiv(2);
hklMerge.l = hklMerge.l + hklMerge.dl/ndiv(3);

E2R = proc.script.ExpandTableToArray(...
        'hklcols',{'h','k','l','I1','sigma1','I2','sigma2'},... % first 3 must be miller indices (h,k,l)
        'symexpand',false);

[~,I1,sigma1,I2,sigma2] = E2R.run(hklMerge,P);

clear hklMerge

%%

% write attributes
h5addgroup(h5filename,loc);

h5writeatt(h5filename,loc,'Description','Merged elastic scattering intensity on a fine grid');
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
h5create(h5filename,fullfile(loc,'I1'),arraysize,opts{:});
h5create(h5filename,fullfile(loc,'sigma1'),arraysize,opts{:});
h5create(h5filename,fullfile(loc,'I2'),arraysize,opts{:});
h5create(h5filename,fullfile(loc,'sigma2'),arraysize,opts{:});

h5writeatt(h5filename,fullfile(loc,'I'),'Description','total intensity');
h5writeatt(h5filename,fullfile(loc,'sigma'),'Description','standard error of I');

h5writeatt(h5filename,fullfile(loc,'I1'),'Description','total intensity, merged from random half-dataset #1');
h5writeatt(h5filename,fullfile(loc,'sigma1'),'Description','standard error of I1');

h5writeatt(h5filename,fullfile(loc,'I2'),'Description','total intensity, merged from random half-dataset #2');
h5writeatt(h5filename,fullfile(loc,'sigma2'),'Description','standard error of I2');


% write data
h5write(h5filename,fullfile(loc,'I'),single(I));
h5write(h5filename,fullfile(loc,'sigma'),single(sigma));
h5write(h5filename,fullfile(loc,'I1'),single(I1));
h5write(h5filename,fullfile(loc,'sigma1'),single(sigma1));
h5write(h5filename,fullfile(loc,'I2'),single(I2));
h5write(h5filename,fullfile(loc,'sigma2'),single(sigma2));

f_new = dir(h5filename);

fprintf(1,'\b... %.3f Mb written in %.1f seconds\n',(f_new.bytes - f_old.bytes)/2^20,toc)


end