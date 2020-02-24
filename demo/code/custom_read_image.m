function [im,h] = custom_read_image(fn)
im = imread(fn);
if nargout==2
    h = custom_read_header(fn);
end
end