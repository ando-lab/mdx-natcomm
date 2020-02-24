classdef Filter
    properties
        countHistogram
        neighborList
    end
    properties(Dependent = true)
        numRows
        maxCount
    end
    methods
        function obj = Filter()
        end
        function val = get.numRows(obj)
            val = size(obj.countHistogram,1);
        end
        function val = get.maxCount(obj)
            val = size(obj.countHistogram,2)-1;
        end
        
        function wm = medianFilter(obj)

            wm = NaN*ones(obj.numRows,1);   % weighted median
            
            countsPerVoxel = obj.countHistogram*(0:obj.maxCount)';
            pixelsPerVoxel = sum(obj.countHistogram,2);
            
            for i=1:obj.numRows
                % calculate weighted median
                cph = countsPerVoxel(obj.neighborList{i});
                npx = pixelsPerVoxel(obj.neighborList{i});
                wm(i) = obj.weightedMedian(cph./npx,npx);
            end
            
        end
        
        function [dkl,deltaDkl] = divergenceFromPoisson(obj,wm)

            % initialize outputs
            dkl = zeros(obj.numRows,1);     % kl divergence of region count histogram vs poisson distribution
            deltaDkl = zeros(obj.numRows,1);% change in total KL divergence of dataset when masking this voxel
            
            % loop over all window regions
            for i=1:obj.numRows
                neighbors = obj.neighborList{i};
                
                % count histogram of all voxels in region
                ch = obj.countHistogram(neighbors,:);
                
                % calculate log(Q(lambda)) where Q is the poisson distribution and
                % lambda = weighted median of neighbors
                logq = -wm(i) + [0,cumsum(log(wm(i)) - log(1:obj.maxCount))];
                
                % calculate KL divergence of full histogram
                
                chtot = sum(ch,1);
                pj = chtot/sum(chtot,2);
                dklterms = pj.*log(pj) - pj.*logq;
                dkl(i) = sum(dklterms(~isnan(dklterms)));
                
                % calculate KL divergence with each member removed
                cprime = repmat(chtot,length(neighbors),1) - ch;
                pj = cprime.*repmat(1./sum(cprime,2),1,size(cprime,2));
                dklprime = pj.*log(pj) - pj.*repmat(logq,size(cprime,1),1);
                dklprime(isnan(dklprime)) = 0;
                dklprime = sum(dklprime,2);
                
                % add change in KL divergence to deltaDkl
                deltaDkl(neighbors) = deltaDkl(neighbors) + dklprime - dkl(i);
            end

        end

        function [isOutlier,dklOut,meanCounts] = maskPoissonOutliers(obj,...
                wm,dkl,deltaDkl)
            
            % order of voxels from worst to best
            [~,ixorder] = sort(deltaDkl,'ascend');
            
            chtot = zeros(obj.numRows,obj.maxCount+1);
            
            % calculate chtot
            for i=1:obj.numRows
                thisind = obj.neighborList{i};
                ch = obj.countHistogram(thisind,:);
                chtot(i,:) = sum(ch,1);
            end
            
            % sequential masking
            for n = 1:length(ixorder)
                i = ixorder(n);
                thisind = obj.neighborList{i};
                
                ch = obj.countHistogram(i,:);
                
                cprime = chtot(thisind,:) - repmat(ch,length(thisind),1);
                pj = cprime.*repmat(1./sum(cprime,2),1,size(cprime,2));
                
                % this is complicated-looking only because I vectorized it
                logq = repmat(-wm(thisind),1,obj.maxCount+1) + [zeros(length(thisind),1),...
                    cumsum(repmat(log(wm(thisind)),1,obj.maxCount) - ...
                    repmat(log(1:obj.maxCount),length(thisind),1),2)];
                
                dklprime = pj.*log(pj) - pj.*logq;
                dklprime(isnan(dklprime)) = 0;
                dklprime = sum(dklprime,2);
                
                % note: without updating, sum(dklprime - dkl(thisind)) = deltaDkl(i)
                if sum(dklprime - dkl(thisind)) < 0 % that is, Dkl decreases upon masking
                    % mask and update!
                    
                    dkl(thisind) = dklprime;
                    chtot(thisind,:) = cprime;
                else
                    nmax = n - 1;
                    break;
                    
                end
            end
            
            isOutlier = false(obj.numRows,1);
            isOutlier(ixorder(1:nmax)) = true;
            dklOut = dkl;
            meanCounts = (chtot*(0:obj.maxCount)')./sum(chtot,2);
        end
    end
    
    methods(Static = true)
        
        function wm = weightedMedian(countRate,numPixels)
            
            if length(countRate)==1
                wm = countRate;
            else
                
                [countRate,ix] = sort(countRate,'ascend');
                numPixels = numPixels(ix);
                cw = cumsum([0;numPixels(:)]);
                cw = cw/cw(end);
                cwmid = 0.5*cw(1:(end-1)) + 0.5*cw(2:end);
                wm = interp1q(cwmid,countRate,0.5);
            end
        end
       
    end
    
end

