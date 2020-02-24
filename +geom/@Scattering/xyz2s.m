function [sx,sy,sz] = xyz2s(x,y,z,Source)
% s = (s1-s0)
wavelength = Source.wavelength;
s0 = Source.wavevector;

r = sqrt(x.*x + y.*y + z.*z);

sx = ((1/wavelength)*x./r - s0(1));
sy = ((1/wavelength)*y./r - s0(2));
sz = ((1/wavelength)*z./r - s0(3));
end
