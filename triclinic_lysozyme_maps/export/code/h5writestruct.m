function h5writestruct(h5filename,dataname,S)

% adapted from:
% https://support.hdfgroup.org/ftp/HDF5/examples/examples-by-api/matlab/HDF5_M_Examples/h5ex_t_cmpd.m

% translate from latlab datatype to h5 data type
dtypes = {...
    'double','H5T_NATIVE_DOUBLE';
    'single','H5T_NATIVE_FLOAT';
    'int32','H5T_NATIVE_INT';
    'uint32','H5T_NATIVE_UINT';
    'int16','H5T_NATIVE_SHORT';
    'uint16','H5T_NATIVE_USHORT';
    'int8','H5T_NATIVE_SCHAR';
    'uint8','H5T_NATIVE_UCHAR';...
    'cell',{'H5T_C_S1','H5T_VARIABLE'}... % must be a cell string
    }; 

% check struct input to verify that it can be written
[mlclass,sz,fn] = struct_info(S);
assert(all(ismember(mlclass,dtypes(:,1))),'unrecognized data type');

% all dimensions must be the same
assert(all(sz==sz(1)),'all array dimensions must be the same');
dims = sz(1);

% identify h5 data type and calculate size in bytes
h5type = cell(size(dtypes,1),1);
h5size = zeros(size(dtypes,1),1);
for j=1:size(dtypes,1)
    thistype = dtypes{j,2};
    if iscell(thistype) % special case for variable length types (cell strings)
        h5type{j} = H5T.copy(thistype{1});
        H5T.set_size(h5type{j},thistype{2});
    else % normal case
        h5type{j} = H5T.copy(thistype);
    end
    h5size(j) = H5T.get_size(h5type{j});
end

% calculate size and type for the struct fields
[~,ind] = ismember(mlclass,dtypes(:,1));
fieldsize = h5size(ind); % size in bytes
fieldtype = h5type(ind); % h5 data type
fieldoffset = cumsum([0;fieldsize(:)]); % offsets when writing data

% create compound data type in memory
filetype = H5T.create('H5T_COMPOUND',sum(fieldsize));

for n = 1:numel(fn)
    H5T.insert(filetype,fn{n},fieldoffset(n),fieldtype{n});
end

% open h5 file
if isfile(h5filename) % if filename exists, append to it
    fid = H5F.open(h5filename,'H5F_ACC_RDWR','H5P_DEFAULT');
else
    fid = H5F.create(h5filename,'H5F_ACC_EXCL','H5P_DEFAULT','H5P_DEFAULT');
end

space = H5S.create_simple(1,fliplr(dims),[]);

dset = H5D.create (fid, dataname, filetype, space, 'H5P_DEFAULT');

% clean up
H5D.write (dset, filetype, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', S);

H5D.close (dset);
H5S.close (space);
H5T.close (filetype);
H5F.close (fid);

end

function [mlclass,sz,fn] = struct_info(S)

fn = fieldnames(S);
mlclass = struct2cell(structfun(@class,S,'Uni',0));
sz = structfun(@numel,S,'Uni',1);

end