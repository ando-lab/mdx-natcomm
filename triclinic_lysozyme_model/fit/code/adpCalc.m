% OLD VERSION
function Ucalc = adpCalc(P,covMat)
% calculate the ADPs predicted by an elastic model

numPts = size(P,1)/3;

Ucalc = zeros(3*numPts,3);

for j=1:numPts
    ind = (1:3) + (j-1)*3;
    Pj = P(ind,:);
    Ucalc(ind,:) = Pj*covMat*Pj';
end

end
% 
% 
% function Ucalc = adpCalc2(P,covMat)
% %ew version, modeled after getCovariances (faster for multi-group
% % systems)
% 
% numPoints = size(P,1)/3;
% numGroups = size(P,2)/6;
% 
% assert(size(P,2)==size(covMat,1) && size(P,2)==size(covMat,2) && size(covMat,3)==1);
% 
% Pn = mat2cell(P,3*ones(numPoints,1),6*ones(numGroups,1));
% 
% [tmpAtomIndex,groupIndex] = find(cellfun(@nnz,Pn)>0);
% assert(issorted(tmpAtomIndex) & numel(tmpAtomIndex) == numPoints);
% 
% % flatten Pn
% Pn = Pn(sub2ind(size(Pn),tmpAtomIndex,groupIndex));
% 
% Ucalc = zeros(3*numPoints,3);
% 
% for j=1:numPoints
%     
%     covjj = covMat((1:6) + 6*(groupIndex(j)-1),(1:6) + 6*(groupIndex(j)-1));
% 
%     Ucalc((1:3) + 3*(j-1),:) = Pn{j}*covjj*Pn{j}';
% end
% 
% 
% 
% end
% 
% 
% 
% 
% 
