function [hklMerge,ih,isOutlier] = mergeScript(fid,options,ScalingModel,hklTable)

% this might be unused?
hklMerge = unique(hklTable(:,options.mergeByColumns),'rows');

nBatches = length(ScalingModel);

% calculate scales
fprintf(fid,'  calculating scale factors\n');

if numel(ScalingModel) > 1
    [a,b,c,d] = ScalingModel.getScales(hklTable.ix,hklTable.iy,hklTable.iz,...
        hklTable.p,hklTable.s,hklTable.n);
else % numel(ScalingModel) = 1, omit n
    [a,b,c,d] = ScalingModel.getScales(hklTable.ix,hklTable.iy,hklTable.iz,...
        hklTable.p,hklTable.s);
end

fprintf(fid,'  scaling intensities\n');
% apply scaling model
scale = a.*b.*d;
offset = c./b;
hklTable.Isc = hklTable.I./scale - offset;
hklTable.sigmasc = hklTable.sigma./scale;

% merge
fprintf(fid,'  merging intensities\n');
[hklMerge,ih] = proc.scale.merge(hklTable,options.mergeByColumns,...
    {'Isc','I'},{'sigmasc','sigma'},'s');

% residual errors
if options.sigmaCutoff
    fprintf(fid,'  removing outliers\n');
    resid = (hklMerge.I(ih) - hklTable.Isc)./hklTable.sigmasc;
    stdevEquiv = sqrt(mean(resid.^2));
    fprintf(fid,'    mean reduced chi-squared of equivalent reflections (1 is ideal): %g\n',stdevEquiv);
    sigma_cutoff = options.sigmaCutoff*stdevEquiv;
    fprintf(fid,'    removing observations with abs(resid/sigma) > %g\n',sigma_cutoff);
    isOutlier = abs(resid) > sigma_cutoff;
    fprintf(fid,'    %d outliers removed, %d observations remaining\n',nnz(isOutlier),nnz(~isOutlier));
    hklTable.sigmasc(isOutlier) = Inf; % effectively removes these points

    % merge again
    fprintf(fid,'  merging intensities\n');
    [hklMerge,ih] = proc.scale.merge(hklTable,options.mergeByColumns,...
        {'Isc','I'},{'sigmasc','sigma'},'s');
else
    isOutlier = [];
end

end
