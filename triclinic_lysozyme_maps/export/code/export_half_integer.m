function export_half_integer(h5filename,loc)
% write half integer map data to h5 file

%h5filename = 'test.h5';
%loc = '/half_integer';

%% load half integer map and get ASU only
load calc/half_integer_map.mat R I sigma I1 sigma1 I2 sigma2

P = R.PeriodicGrid;

% get the ASU

[h,k,l] = P.grid();

isASU = symm.SpaceGroup(1).LaueGroup.testASU(single(h),single(k),single(l));
isASU = isASU & ~isnan(I);

hklTable = table(h(isASU),k(isASU),l(isASU),...
    I(isASU),sigma(isASU),I1(isASU),sigma1(isASU),I2(isASU),sigma2(isASU),...
    'VariableNames',{'h','k','l','I','sigma','I1','sigma1','I2','sigma2'});

clear h k l isASU P I sigma I1 sigma1 I2 sigma2

E2R = proc.script.ExpandTableToArray(...
        'hklcols',{'h','k','l','I','sigma','I1','sigma1','I2','sigma2'},... % first 3 must be miller indices (h,k,l)
        'symexpand',false,...
        'ndiv',[1,1,1]);

[P,I,sigma,I1,sigma1,I2,sigma2] = E2R.run(hklTable);

% Export to h5 file

h5addgroup(h5filename,loc);

h5writeatt(h5filename,loc,'Description','Total intensity interpolated at half-integer Miller indices');
h5writeatt(h5filename,loc,'units','electron scattering per unit cell');
h5writeatt(h5filename,loc,'grid_coordinates','fractional Miller indices (h,k,l)');
h5writeatt(h5filename,loc,'grid_size',uint32(P.N));
h5writeatt(h5filename,loc,'grid_ori',P.ori);
h5writeatt(h5filename,loc,'grid_delta',P.delta);

arraysize = P.N;
opts = {'Datatype','single','FillValue',single(NaN)};

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