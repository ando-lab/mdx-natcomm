function setup_environment()

mdx_lib_path = '../'; % can be absolute or relative
data_directory = '~/data/lysozyme/nitrate/ubatch4/'; % must be absolute
code_directories = {'calc/code','export/code'};

% 1) add mdx-kit to (global) MATLAB path
addpath(mdx_lib_path);

% 1) add code directories to (global) MATLAB path
addpath(code_directories{:});

% 2) make a symlink to data directory for easy access
unix(sprintf('rm data; ln -s %s data',data_directory));

end
