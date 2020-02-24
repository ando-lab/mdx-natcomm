function h5addgroup(h5filename,groupname)

% adapted from:
% https://support.hdfgroup.org/ftp/HDF5/examples/examples-by-api/matlab/HDF5_M_Examples/h5ex_t_cmpd.m

% open h5 file
if isfile(h5filename) % if filename exists, append to it
    fid = H5F.open(h5filename,'H5F_ACC_RDWR','H5P_DEFAULT');
else
    fid = H5F.create(h5filename,'H5F_ACC_EXCL','H5P_DEFAULT','H5P_DEFAULT');
end

gid = H5G.create(fid,groupname,'H5P_DEFAULT','H5P_DEFAULT','H5P_DEFAULT');

% clean up
H5G.close(gid);
H5F.close(fid);

end
