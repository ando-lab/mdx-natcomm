function tf = write(mapfname,mh,mapdata)
%WRITE binary map file
fid = fopen(mapfname,'w');

fwrite(fid,[mh.nc,mh.nr,mh.ns,mh.mode,mh.ncstart,mh.nrstart,mh.nsstart,mh.nx,mh.ny,mh.nz],'int32');
fwrite(fid,[mh.x_length,mh.y_length,mh.z_length,mh.alpha,mh.beta,mh.gamma],'real*4');
fwrite(fid,[mh.mapc,mh.mapr,mh.maps],'int32');
fwrite(fid,[mh.amin,mh.amax,mh.amean],'real*4');
fwrite(fid,[mh.ispg,mh.nsymbt,mh.lskflg],'int32');
fwrite(fid,mh.skwmat,'real*4');
fwrite(fid,mh.skwtrn,'real*4');
fwrite(fid,mh.future,'int32');
fwrite(fid,[mh.map,mh.machst],'char*1');
fwrite(fid,mh.arms,'real*4');
fwrite(fid,mh.nlabl,'int32');
fwrite(fid,mh.label','char*1');
fwrite(fid,mh.sym','char*1');

% reshape map data
fwrite(fid,permute(mapdata,[mh.mapc,mh.mapr,mh.maps]),'real*4');

fclose(fid);
tf = true;
end

