function idx_ok = struct_search(struct_in,varargin)
idx_ok = ones(size(struct_in));
nsearches = floor(length(varargin)/2);
for i=1:nsearches
    searchfield = varargin{2*i - 1};
    searchstring = varargin{2*i};
    rexpr = ['^' regexptranslate('wildcard',searchstring) '$'];
    search_result = regexp({struct_in.(searchfield)}, ...
        rexpr, 'match');
    idx_ok = idx_ok & ~cellfun(@isempty,search_result);
end

end