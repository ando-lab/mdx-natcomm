function mapheader = initHeader()
% INITHEADER

mapheader.nc = 0;
mapheader.nr = 0;
mapheader.ns = 0;
mapheader.mode = 2;
mapheader.ncstart = 0;
mapheader.nrstart = 0;
mapheader.nsstart = 0;
mapheader.nx = 0;
mapheader.ny = 0;
mapheader.nz = 0;
mapheader.x_length = 0;
mapheader.y_length = 0;
mapheader.z_length = 0;
mapheader.alpha = 90;
mapheader.beta = 90;
mapheader.gamma = 90;
mapheader.mapc = 3;
mapheader.mapr = 1;
mapheader.maps = 2;
mapheader.amin = 0;
mapheader.amax = 0;
mapheader.amean = 0;
mapheader.ispg = 1;
mapheader.nsymbt = 80;
mapheader.lskflg = 0; %?
mapheader.skwmat = zeros([9,1]); % S11, S12, S13, S21 etc
mapheader.skwtrn = zeros([3,1]);
mapheader.future = zeros([15,1]);
mapheader.map = 'MAP '; 
mapheader.machst = ['DD',char([0,0])]; % little endian?
mapheader.arms = 0;
mapheader.nlabl = 0;
mapheader.label = repmat(' ',[10,80]);
mapheader.sym = repmat(' ',[1,80]);

end