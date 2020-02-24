function [Iout,sigmaOut] = regrid(Iin,sigmaIn,SG,Rin,Rout,nmin,kernel)
% REGRID
%
% Before running:
% 
% 1) load the fine map (such as 7x7x5)
%
% 2) symmetry-expand the map onto a full grid
%
% REGRID does the following:
%
% 3) create a coarse grid (such as 4x4x4)
%
% 4) create a kernel
%
% 5) for all points in the coarse grid, find nearest point in the fine
% grid, calculate x,y,z to all points in the kernel, get intensities and
% sigmas from the fine grid, fit to a polynomial, return the value of the
% polynomial at the center
%
% 6) return the result

if nargin< 6 || isempty(nmin)
    % minimum number of points to attempt an interpolation
    nmin = 21;
end

if nargin >= 7 && ~isempty(kernel)
    k1 = kernel(:,1);
    k2 = kernel(:,2);
    k3 = kernel(:,3);
else
    % 3x3x3 kernel
    [k1,k2,k3] = latt.PeriodicGrid([3,3,3],[-1,-1,-1],[3,3,3]).grid();
end

% do the interpolation

% make the output grid
[h,k,l] = Rout.PeriodicGrid.grid();
[sx,sy,sz] = Rout.Basis.frac2lab(h,k,l);

% test if is ASU
isASU = SG.LaueGroup.testASU(round(h,10),round(k,10),round(l,10));

Iout = NaN*ones(Rout.N);
sigmaOut = Inf*ones(Rout.N);

% interpolate
fprintf(1,'INTERPOLATING ON A COARSE GRID: %5.1f%%\n',0);
for j=1:numel(sx)
    if ~isASU(j)
        continue;
    end
    if ~mod(j,round(numel(sx)/1000))
        fprintf(1,'\b\b\b\b\b\b\b%5.1f%%\n',100*j/numel(sx));
    end

% calculate the neighborhood of grid points
[n1,n2,n3] = Rin.lab2ind(sx(j),sy(j),sz(j),false);
n123 = [n1 + k1(:),n2 + k2(:),n3 + k3(:)];
isincl = n123(:,1) > 1 & n123(:,2) >= 1 & n123(:,3) >= 1 & ...
    n123(:,1) <= Rin.N(1) & n123(:,2) <= Rin.N(2) & n123(:,3) <= Rin.N(3);
n123 = n123(isincl,:);
ind = sub2ind(Rin.N,n123(:,1),n123(:,2),n123(:,3));

% remove any points where I is NaN
isincl = ~isnan(Iin(ind));
n123 = n123(isincl,:);
ind = ind(isincl);

if numel(ind) < nmin
    continue;
end

% calculate the distance
[x,y,z] = Rin.ind2lab(n123(:,1),n123(:,2),n123(:,3));
x = x - sx(j);
y = y - sy(j);
z = z - sz(j);

B = zeros(numel(x),10);
B(:,1) = 1;
B(:,2) = x;
B(:,3) = y;
B(:,4) = z;
B(:,5) = x.*x;
B(:,6) = x.*y;
B(:,7) = x.*z;
B(:,8) = y.*y;
B(:,9) = y.*z;
B(:,10) = z.*z;

Ik = Iin(ind);
sigmak = sigmaIn(ind);

B1 = B.*repmat(1./sigmak,1,10);

cfit = B1\(Ik./sigmak); % the solution
Iout(j) = cfit(1); % the value at x=y=z=0

% error propagation
V = inv(B1'*B1);
sigmaOut(j) = real(sqrt(V(1,1))); % the variance at x=y=z=0

end

% expand asu to symmetry related points

T = table(h(isASU),k(isASU),l(isASU),Iout(isASU),sigmaOut(isASU),...
    'VariableNames',{'h','k','l','I','sigma'});

E2A = proc.script.ExpandTableToArray(...
    'hklcols',{'h','k','l','I','sigma'},...
    'SpaceGroup',SG,...
    'symexpand',true,...
    'ndiv',[]... % unused
    );

[~,Iout,sigmaOut] = E2A.run(T,Rout.PeriodicGrid);

end