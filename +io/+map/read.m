function [mapdata,mapheader] = read(mapfname)
% READ binary map file
fid = fopen(mapfname,'r');

mapheader.nc = fread(fid,[1,1],'int32');
mapheader.nr = fread(fid,[1,1],'int32');
mapheader.ns = fread(fid,[1,1],'int32');
mapheader.mode = fread(fid,[1,1],'int32');
mapheader.ncstart = fread(fid,[1,1],'int32');
mapheader.nrstart = fread(fid,[1,1],'int32');
mapheader.nsstart = fread(fid,[1,1],'int32');
mapheader.nx = fread(fid,[1,1],'int32');
mapheader.ny = fread(fid,[1,1],'int32');
mapheader.nz = fread(fid,[1,1],'int32');
mapheader.x_length = fread(fid,[1,1],'real*4');
mapheader.y_length = fread(fid,[1,1],'real*4');
mapheader.z_length = fread(fid,[1,1],'real*4');
mapheader.alpha = fread(fid,[1,1],'real*4');
mapheader.beta = fread(fid,[1,1],'real*4');
mapheader.gamma = fread(fid,[1,1],'real*4');
mapheader.mapc = fread(fid,[1,1],'int32');
mapheader.mapr = fread(fid,[1,1],'int32');
mapheader.maps = fread(fid,[1,1],'int32');
mapheader.amin = fread(fid,[1,1],'real*4');
mapheader.amax = fread(fid,[1,1],'real*4');
mapheader.amean = fread(fid,[1,1],'real*4');
mapheader.ispg = fread(fid,[1,1],'int32');
mapheader.nsymbt = fread(fid,[1,1],'int32');
mapheader.lskflg = fread(fid,[1,1],'int32');
mapheader.skwmat = fread(fid,[9,1],'real*4'); % S11, S12, S13, S21 etc
mapheader.skwtrn = fread(fid,[3,1],'real*4');
mapheader.future = fread(fid,[15,1],'int32');
%mapheader.map = fread(fid,[1,1],'int32');
mapheader.map = fread(fid,[1,4],'char*1=>char');
mapheader.machst = fread(fid,[1,4],'char*1=>char');
mapheader.arms = fread(fid,[1,1],'real*4');
mapheader.nlabl = fread(fid,[1,1],'int32');
mapheader.label = fread(fid,[80,10],'char*1=>char')';
mapheader.sym = fread(fid,[80,mapheader.nsymbt/80],'char*1=>char')';

assert(mapheader.mode==2,'only mode 2 maps have been implemented');
npts = mapheader.nc*mapheader.nr*mapheader.ns; % number of data points to read
mapdata = fread(fid,npts,'real*4');

fclose(fid);

% reshape map data
mapdata = reshape(mapdata,[mapheader.nc,mapheader.nr,mapheader.ns]);

mapdata = ipermute(mapdata,[mapheader.mapc,mapheader.mapr,mapheader.maps]);
end

