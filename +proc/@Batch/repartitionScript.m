function [bin, corr] = repartitionScript(fid,options,Grid,bin,corr)

dzmax = options.dzmax;

for j=1:(length(bin)-1)

bin1 = bin{j};
bin2 = bin{j+1};
Grid1 = Grid(j);
Grid2 = Grid(j+1);

[~,ia,ib] = intersect(Grid1,Grid2);

[x,y] = ndgrid(ia,1:Grid1.arraysize(2));
ind1 = Grid1.sub2ind(x,y);
[x,y] = ndgrid(ib,1:Grid2.arraysize(2));
ind2 = Grid2.sub2ind(x,y);

hasPixels1 = reshape(bin1.pixels > 0,Grid1.numRows,[]);
hasPixels2 = reshape(bin2.pixels > 0,Grid2.numRows,[]);
hasPixels12 = hasPixels1(ia,:) & hasPixels2(ib,:);
ind1 = ind1(hasPixels12);
ind2 = ind2(hasPixels12);

pixels1 = bin1.pixels(ind1);
pixels2 = bin2.pixels(ind2);

counts1 = bin1.counts(ind1);
counts2 = bin2.counts(ind2);

iz1 = bin1.iz(ind1);
iz2 = bin2.iz(ind2);

pixels = pixels1 + pixels2;
counts = counts1 + counts2;
iz = (iz1.*pixels1 + iz2.*pixels2)./(pixels1 + pixels2);

t12 = table(counts,pixels,iz);
tz = table(counts*0,pixels*0,iz*NaN);

corr12 = corr{j}(ind1,:);
corr12nan = corr12;
corr1 = corr{j}(ind1,:);
corr2 = corr{j+1}(ind2,:);
for col=corr12.Properties.VariableNames
    corr12.(col{1}) = (corr1.(col{1}).*pixels1 + ...
        corr2.(col{1}).*pixels2)./(pixels1+pixels2);
    corr12nan.(col{1}) = corr12nan.(col{1})*NaN;
end
corr12.bkgErr = sqrt((corr1.bkgErr.^2.*pixels1 + ...
        corr2.bkgErr.^2.*pixels2)./(pixels1+pixels2));

move1to2 = pixels1 < pixels2 & abs(iz2 - iz1) < dzmax;
move2to1 = pixels1 > pixels2 & abs(iz2 - iz1) < dzmax;

bin{j+1}(ind2(move1to2),:) = t12(move1to2,:);
bin{j}(ind1(move1to2),:) = tz(move1to2,:);
corr{j+1}(ind2(move1to2),:) = corr12(move1to2,:);
corr{j}(ind1(move1to2),:) = corr12nan(move1to2,:);

bin{j}(ind1(move2to1),:) = t12(move2to1,:);
bin{j+1}(ind2(move2to1),:) = tz(move2to1,:);

corr{j}(ind1(move2to1),:) = corr12(move2to1,:);
corr{j+1}(ind2(move2to1),:) = corr12nan(move2to1,:);

end

end