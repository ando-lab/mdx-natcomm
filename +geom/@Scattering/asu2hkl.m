function [hnew,knew,lnew,ind] = asu2hkl(h,k,l,Crystal)

lg = Crystal.SpaceGroup.PointGroup.LaueGroup;
Ops = lg.generalPositionsMat;
nOps = length(Ops);

newhkl = zeros(nOps,numel(h),3);
ind = repmat((1:numel(h)),nOps,1);

for j=1:length(Ops)
    newhkl(j,:,:) = [h(:),k(:),l(:)]*Ops{j};
end

hnew = reshape(newhkl(:,:,1),[],1);
knew = reshape(newhkl(:,:,2),[],1);
lnew = reshape(newhkl(:,:,3),[],1);
ind = reshape(ind,[],1);

uniqueVals = unique([ind,hnew,knew,lnew],'rows','sorted');

ind = uniqueVals(:,1);
hnew = uniqueVals(:,2);
knew = uniqueVals(:,3);
lnew = uniqueVals(:,4);

end