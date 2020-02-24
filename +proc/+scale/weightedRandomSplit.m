function isgrp1 = weightedRandomSplit(h,w)
% weightedRandomSplit - randomly splits the observations into two groups
%
% returns a logical vector isgrp1 which is true if an observation is in
% group 1, and false if it is in group2
%
% which group an observation finds itself in is random, except that the
% observational weights are approximately equal. In other words, if we
% calculate
%
% w1 = accumarray(h(isgrp1),w(isgrp1),Nh,1);
% w2 = accumarray(h(~isgrp1),w(~isgrp1),Nh,1);
%
% then for each type of observation, w1 ~ w2

% pre-index by h
[sortedIndex,indexOrder] = sort(h);
[hvals,indFirst] = unique(sortedIndex,'stable');
Nh = size(hvals,1);
indLast = [indFirst(2:end)-1; length(h)];

% indicates which group a reflection belongs to (0,1)
isgrp1 = false([length(h),1]); 

for thish = 1:Nh
ix = indexOrder(indFirst(thish):indLast(thish));

thisw = w(ix);
thisorder = randperm(length(thisw));

% the following is a little inelegant, but it's reasonably fast so I won't
% worry about improving it.
nsplit = find(cumsum(thisw(thisorder))/sum(thisw) < 0.5, 1, 'last');
if isempty(nsplit)
    nsplit = 1; % covers the case where the first data point has more than half of the weight
end
isincl = false(length(thisorder),1);
isincl(thisorder(1:nsplit)) = true;
isincl = xor(randi(2)==1,isincl); % randomly flip which group is included

isgrp1(ix(isincl)) = true;

end

end