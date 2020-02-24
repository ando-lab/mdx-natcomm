function tf = has_pilatus_cbf_convention(obj)
% Check for the _array_data.header_convention data item
pilatus_header_pattern = ['_array_data.header_convention ',...
    '+["'']?(SLS|PILATUS)_\d+(\.?\d*)*["'']?'];
tf = logical(regexp(obj.rawheader,pilatus_header_pattern));
end