function [countHist,countOverflow] = countHistogram(obj,Grid,maxCount)
ds = obj.ImageSeries.datastore;
[countHist,countOverflow] = Grid.countHistogram(maxCount);
while hasdata(ds)
    [im,iminfo] = ds.read;
    if obj.verbose, disp(iminfo); end
    frameNumber = iminfo.Label;
    [hh,kk,ll] = obj.WedgeGeometry.frame2hkl(frameNumber);
    [ch,co] = Grid.countHistogram(maxCount,im,hh,kk,ll);
    countHist = countHist + ch;
    countOverflow = countOverflow + co;
end
end