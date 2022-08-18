function y = sumOfGaussians(x,coeffs)
ncoeffs = length(coeffs);

if floor(ncoeffs/2)==ncoeffs/2 % number of coeffs is even: no offset
    aj = coeffs(1:2:(end-1));
    bj = coeffs(2:2:end);
    c = 0;
else % number of coeffs is odd: last one is the offset
    aj = coeffs(1:2:(end-2));
    bj = coeffs(2:2:(end-1));
    c = coeffs(end);
end

ngauss = length(aj);
y = c*ones(size(x));
x2 = x.*x;
for j=1:ngauss
    y = y + aj(j)*exp(-bj(j)*x2);
end

end