function img = read_cbf_image(fname,fpos,nbytes,nrows,ncols)
    cbf_data = read_cbf_data(fname,fpos,nbytes);
    img = cbf_decompress_mex(cbf_data,nrows,ncols);
end

function [cbf_data] = read_cbf_data(fname,fpos,nbytes)
    a = fopen(fname,'r');
    fseek(a,fpos,-1);
    cbf_data = uint8(fread(a,nbytes,'uint8'));
    fclose(a);
end