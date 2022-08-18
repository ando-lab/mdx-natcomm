classdef ScaleToReference < util.propertyValueConstructor
    % SCALETOREFERENCE - find a single scale factor that brings intensities
    % from an mtz file into best agreement with reference data (in a least
    % squares sense, with outlier rejection).
    properties
        refcols = {'hasu','kasu','lasu','I','sigma'};
        
        mtzIn = 'XDS_ASCII_aimless.mtz';
        mtzcols = {'H','K','L','IMEAN','SIGIMEAN'};
        
        outcols = {'h','k','l','I','sigma'}; % column names for output table
            
        map2asu = true; % re-map hkl values of the mtz to the ASU
        
        nIter = 10 % number of iterations for robust linear regression
        outlierCutoff = 15 % cutoff in a maximum absolute deviation, weighted by sigma of the reference
        
        fid = 1; % where to print stats. 1 is command window. can also print to log file.
    end
    
    methods
        function obj = ScaleToReference(varargin)
            obj@util.propertyValueConstructor(varargin{:});
        end
        
        function [mtzOut,ref,scaleFactor] = run(obj,ref)
            
            
            ref = ref(:,obj.refcols);
            ref.Properties.VariableNames = obj.outcols;
            
            % get table from mtzIn 
            [mtzData, colNames, c] = io.mtz.read(obj.mtzIn);
            Crystal = geom.Crystal(c);
            [ism,ix] = ismember(obj.mtzcols,colNames);
            assert(all(ism),'mtz cols not found');
            mtzOut = array2table(mtzData(:,ix),'VariableNames',obj.outcols);
            
            if obj.map2asu
                [mtzOut.(1),mtzOut.(2),mtzOut.(3)] = Crystal.hkl2asu(...
                    mtzOut.(1),mtzOut.(2),mtzOut.(3));
            end
            
            % add mtz columns to ref table
            [~,ia,ic] = intersect(...
                table2array(ref(:,1:3)),...
                table2array(mtzOut(:,1:3)),'rows');
            
            nrows = size(ref,1);
            ref.I_mtz = NaN*ones(nrows,1);
            ref.sigma_mtz = NaN*ones(nrows,1);
            ref.I_mtz(ia) = mtzOut.I(ic);
            ref.sigma_mtz(ia) = mtzOut.sigma(ic);
            
            ref.isIncluded = false(nrows,1);
            ref.isIncluded(ia) = true;
            ref.isIncluded(isnan(ref.I)) = false; % just in case
            
            x = ref.I(ref.isIncluded);
            y = ref.I_mtz(ref.isIncluded);
            scaleFactor = 1./(y\x);
            
            % refine the scale factor if Niter > 0
            fprintf(obj.fid,'   %-8s %-8s %-8s %-8s %-8s\n','iter','cc','scale','nbad','ngood');
            nbad = update_stats(obj.fid,0,ref.isIncluded,x,y,scaleFactor);
           
            for j=1:obj.nIter
                
                resid = (ref.I - ref.I_mtz/scaleFactor)./ref.sigma;
                medianAbsDev = median(abs(resid(ref.isIncluded)));
        
                isOutlier = abs(resid) > obj.outlierCutoff*medianAbsDev;
                ref.isIncluded = ref.isIncluded & ~isOutlier;
                
                x = ref.I(ref.isIncluded);
                y = ref.I_mtz(ref.isIncluded);
                
                scaleFactor = 1./(y\x);
                
                nbad_new = update_stats(obj.fid,j,ref.isIncluded,x,y,scaleFactor);
                
                if nbad_new == nbad
                    break;
                end
                nbad = nbad_new;
            end
            
            function nbad = update_stats(fid,iter,isGood,x,y,scale) % a status update
                cc = sum((x-mean(x)).*(y-mean(y)))/sqrt(sum((x-mean(x)).^2).*sum((y-mean(y)).^2));
                nbad = sum(~isGood);
                ngood = sum(isGood);
                fprintf(fid,'   %-8g %-8g %-8g %-8g %-8g\n',iter,cc,scale,nbad,ngood);
            end

            % apply scale factor to original data table
           mtzOut.I = mtzOut.I/scaleFactor;
           mtzOut.sigma = mtzOut.sigma/scaleFactor;
        end
    end
end