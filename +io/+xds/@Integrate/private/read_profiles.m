function wedge = read_profiles(fileName,nptsAlphaBeta,nptsGamma)

fid = fopen(fileName,'r');

%xparm_in = read_header_xparm(fid);

j = 1;
wedge = struct('dataRange',{},'profile',{});

while ~feof(fid)
    newline = fgetl(fid);
    if iswedgestart(newline)
        tmp_cell= regexp(newline,...
            '^\s+PROCESSING OF IMAGES\s+(\d+) \.{3}\s+(\d+)\s*$','tokens');
        n1 = str2double(tmp_cell{1}{1});
        n2 = str2double(tmp_cell{1}{2});
        wedge(j).dataRange = [n1,n2];
    end
    
    if iswedgeend(newline)
        wedge(j).profile = read_profile(fid,nptsAlphaBeta,nptsGamma);
        j = j + 1;
    end
end
end

function tf = iswedgestart(newline)
tf = ~isempty(regexp(newline,'^\s+PROCESSING OF IMAGES','once'));
end

function tf = iswedgeend(newline)
tf = ~isempty(regexp(newline,'^\s*\*{5} AVERAGE','once'));
end


function profile = read_profile(fid,npts_alpha_beta,npts_gamma)
profile = zeros(npts_alpha_beta,npts_alpha_beta,npts_gamma);
k_start = 1;
while ~feof(fid) && k_start<= npts_gamma;
    newline = fgetl(fid);
    if ~isempty(newline)
        for j=1:npts_alpha_beta
            regexp_split = regexp(newline,...
                '(?<=\s{3})([\s\d-]{2}\d)+(?=(\s{3}|$))','match');
            vals = cellfun(@(x) str2num(reshape(x,3,[])')',...
                regexp_split,'UniformOutput',0);
            for k = 1:length(vals);
                profile(j,:,k + k_start - 1) = vals{k};
            end
            next_k_start = k_start + length(vals);
            newline = fgetl(fid);
        end
        k_start = next_k_start;
    end
end
end