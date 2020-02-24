function tf = setup_environment()

mdx_lib_path = '../'; % change if demo folder is moved elsewhere (recommended)

% 1) add mdx-kit to (global) MATLAB path
addpath(mdx_lib_path);

% 2) add code directory to (global) MATLAB path
addpath('code');

tf = true;

end
