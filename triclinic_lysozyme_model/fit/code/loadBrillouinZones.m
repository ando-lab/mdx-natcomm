function [hklGrid,I,sigma] = loadBrillouinZones(h5filename,h5dataset,hklRef)

grid_delta = double(h5readatt(h5filename,h5dataset,'grid_delta'))';
grid_ori = double(h5readatt(h5filename,h5dataset,'grid_ori'))';
grid_size = double(h5readatt(h5filename,h5dataset,'grid_size'))';

R = latt.PeriodicGrid(grid_size,grid_ori,grid_delta.*grid_size);
ndiv = 1./grid_delta;

I = NaN*ones([size(hklRef,1),ndiv(1),ndiv(2),ndiv(3)]);
sigma = NaN*ones([size(hklRef,1),ndiv(1),ndiv(2),ndiv(3)]);

% load data from h5 file
for j=1:size(hklRef,1)
    [n1,n2,n3] = R.frac2ind(hklRef(j,1),hklRef(j,2),hklRef(j,3),false);
    start = [n1,n2,n3] - 0.5*(ndiv-1);
    count = ndiv;
    try
        I1 = h5read(h5filename,fullfile(h5dataset,'I'),start,count);
        sigma1 = h5read(h5filename,fullfile(h5dataset,'sigma'),start,count);
    catch EM
        warning('error thrown by h5read for h,k,l %d,%d,%d:\n message = ''%''',...
            hklRef(j,1),hklRef(j,2),hklRef(j,3),EM.message);
        continue;
    end
    I(j,:,:,:) = I1;
    sigma(j,:,:,:) = sigma1;
end

sigma(isnan(I) | isnan(sigma)) = Inf;
I(isinf(sigma)) = 0;


% make grids for output
R1 = latt.PeriodicGrid(ndiv,-.5*(ndiv-1)./ndiv,[1,1,1]);
hklGrid = latt.PeriodicGrid.empty();
for j=1:size(hklRef,1)
    hklGrid(j,1) = R1;
    hklGrid(j,1).ori = hklGrid(j).ori + hklRef(j,:);
end

end


