function coeffs = fitSumOfGaussians(x,y,ngauss,interptf,xmax)

if nargin<=3 || isempty(interptf)
    interptf = false;
end
if nargin<=4 || isempty(xmax)
    xmax = max(x);
end

ixmax = find(x<=xmax,1,'last');

if interptf
    dx = min(diff(x));
    xx = min(x):dx:x(ixmax);
    yy = interp1(x,y,xx,'pchip');
else
    xx = x(1:ixmax);
    yy = y(1:ixmax);
end
    
a0 = yy(1)*ones(1,ngauss)/ngauss;
b0 = logspace(2,0,ngauss);
c0 = 0;

coeffs0 = [a0;b0];
coeffs0 = [coeffs0(:);c0]';

coeffs = lsqnonlin(@(x) yy-sumOfGaussians(xx,x),coeffs0,0*coeffs0);

end