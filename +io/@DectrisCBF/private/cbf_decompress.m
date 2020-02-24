function img = cbf_decompress(fdata,nrows,ncols)
%
% the cpf files from CHESS are encoded using a "byte-offset" algorithm. the
% algorithm is described here;
% http://www.bernstein-plus-sons.com/software/CBF/doc/CBFlib.html#3.3
%
% I implemented a minimal version of the decompression algorithm
img = int32(zeros(nrows,ncols));

% takes ~ 1 minute to de-compress image
d8 = typecast(uint8(hex2dec('80')),'int8');
d16 = typecast(uint16(hex2dec('8000')),'int16');
%d32 = typecast(uint32(hex2dec('80000000')),'int32');

base_pixel = int32(0); % does the type make sense?
i = 1; j = 1;
%
while i <= length(fdata) % was i < length(fdata). Changed to <= on Sept 26, 2019
    delta8 = typecast(fdata(i),'int8');
    delta32 = int32(delta8);
    i = i+1;
    if delta8 == d8
        delta16 = typecast(fdata(i:(i+1)),'int16');
        delta32 = int32(delta16);
        %delta64 = int64(delta16);
        i = i+2;
        if delta16 == d16
            delta32 = typecast(fdata(i:(i+3)),'int32');
            %delta64 = int64(delta32);
            i = i + 4;
%             if delta32==d32
%                 delta64 = typecast(fdata(i:(i+7)),'int64');
%                 i = i+8;
%             end
        end
    end
    base_pixel = base_pixel + delta32;
    img(j) = base_pixel;
    j = j + 1;
end

end