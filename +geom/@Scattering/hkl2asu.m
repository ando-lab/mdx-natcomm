function [hnew,knew,lnew] = hkl2asu(h,k,l,Crystal)
% map h, k, l to reciprocal asymmetric unit
% using symmetry operators for the Laue group

lg = Crystal.SpaceGroup.PointGroup.LaueGroup;
Ops = lg.generalPositionsMat;

hnew = NaN*zeros(size(h));
knew = NaN*zeros(size(h));
lnew = NaN*zeros(size(h));

for j=1:length(Ops)
    newhkl = [h(:),k(:),l(:)]*Ops{j};
    isasu = lg.testASU(newhkl(:,1),newhkl(:,2),newhkl(:,3));
    hnew(isasu) = newhkl(isasu,1);
    knew(isasu) = newhkl(isasu,2);
    lnew(isasu) = newhkl(isasu,3);
end

end