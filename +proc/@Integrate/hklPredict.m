function [h,k,l] = hklPredict(obj,smax)
frameRange = obj.ImageSeries.frameRange;
if nargin<2 || isempty(smax) || isinf(smax)
    [x,y,z] = obj.WedgeGeometry.Detector.xyzEdges;
    [sx,sy,sz] = geom.Scattering.xyz2s(x,y,z,obj.WedgeGeometry.Source);
    smax = sqrt(max(sx.*sx + sy.*sy + sz.*sz));
end
allhkl = zeros(0,3);
for j = frameRange(1):frameRange(2)
    phiMid = obj.WedgeGeometry.Spindle.phiMid(j);
    % find h,k,l of cells intersecting ewald sphere
    [h,k,l] = geom.Scattering.ewald2hkl(...
        obj.WedgeGeometry.Source,obj.WedgeGeometry.Crystal,...
        obj.WedgeGeometry.Spindle, phiMid,smax);
    allhkl = unique([allhkl;[h,k,l]],'rows');
end
h = allhkl(:,1); k = allhkl(:,2); l = allhkl(:,3);
end