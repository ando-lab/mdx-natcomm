function setup_environment()

mdx_lib_path = '../'; % can be absolute or relative
code_directories = {'fit/code','model/code'};

% If you have already exported the map with mdx, then you can use this path.
% If not, download the map h5 file at CXIDB ID 128 and enter the path here:
map_data = '../triclinic_lysozyme_maps/export/triclinic_lysozyme_maps.h5';

% 1) add mdx-kit to (global) MATLAB path
addpath(mdx_lib_path);

% 1) add code directories to (global) MATLAB path
addpath(code_directories{:});

% 2) make a symlink to the map data for easy access
unix(sprintf('rm triclinic_lysozyme_maps.h5; ln -s %s triclinic_lysozyme_maps.h5',map_data));

end
