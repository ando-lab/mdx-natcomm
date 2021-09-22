function [tf, errorMessage] = scale(varargin)

options = struct(...
    'workingDirectory','./',...
    'combineIn','combine.mat',...
    'scaleIn','combine.mat',...
    'matOut','scale.mat',...
    'logOut','scale.log',...
    'mergeByColumns',[],...
    'program',[],...
    'cMin',0,...
    'bizMult',1,...
    'csMult',100,...
    'cizMult',10,...
    'cMult',0.001,...
    'aixiyMult',10,...
    'aizMult',1,...
    'dipMult',0.01);

options.mergeByColumns = {'hasu','kasu','lasu'};
options.program = {'b',50,1E-4,'o',[],10,'cb',20,1E-3,'o',[],5,'c',20,1E-3,'a',10,1E-4,'d',5,1E-4};
% o,[],5 means remove outliers > 5 sigma

BatchProcess = proc.Batch('proc.Batch.scale',options,varargin{:});
options = BatchProcess.options; % post-processed options

try % START MAIN SCRIPT
        
    BatchProcess.start();

    hklTable = BatchProcess.readFromMatFile(options.combineIn,...
        'diffuseTable');
    
    ScalingModel = BatchProcess.readFromMatFile(options.scaleIn,'ScalingModel');

    fid =1; % direct output
    
    hklMerge = unique(hklTable(:,options.mergeByColumns),'rows');
    nBatches = length(ScalingModel);
    
    fprintf(fid,'  initializing proc.scale.MultiBatchScaler\n');
    B = proc.scale.MultiBatchScaler.empty();
    for j=1:nBatches
        t1 = hklTable(hklTable.n==j,:);
        [~,ih] = ismember(t1(:,options.mergeByColumns),hklMerge,'rows');
        B(j) = proc.scale.MultiBatchScaler('I',t1.I,'sigma',t1.sigma,...
            'ih',ih,'ihmax',size(hklMerge,1),...
            'Model',proc.scale.ScaleFactors('sigma',t1.sigma,...
            'ix',t1.ix,'iy',t1.iy,'iz',t1.iz,'ip',t1.p,'s',t1.s,...
            'ScalingModel',ScalingModel(j)));
    end
    % free up some memory
    clear hklTable
    
    fprintf(fid,'  calculating regularization parameters\n');
    
    % initial merge
    u = B.merge();
    
    % initialize regularization constants
    [diagAA,diagBB] = B.bCalcDiag(u);
    bizLambda = options.bizMult * sum(diagAA)/sum(diagBB);
    
    [diagAA,diagBBx,diagBBy,diagH] = B.cCalcDiag(u);
    csLambda = options.csMult * sum(diagAA)/sum(diagBBx);
    cizLambda = options.cizMult * sum(diagAA)/sum(diagBBy);
    cLambda = options.cMult*max(diagAA)/mean(diagH); % put pressure on hightest point.
    
    [diagAA,diagBBxy,diagBBz] = B.aCalcDiag(u);
    aixiyLambda = options.aixiyMult * sum(diagAA)/sum(diagBBxy);
    aizLambda = options.aizMult * sum(diagAA)/sum(diagBBz);
    
    [diagAA,diagH] = B.dCalcDiag(u);
    dipLambda = options.dipMult * max(diagAA)/mean(diagH); % put pressure on hightest point
    
    fprintf(fid,'  refining scales\n');
    
    for k = 1:3:(length(options.program)-2)
        refParam = options.program{k}; % a,b,c,d,cb
        nIter = options.program{k+1};
        x2tol = options.program{k+2};
        switch lower(refParam)
            case 'a' % refine a
                fprintf(fid,'    refining scale vs. frame number and detector position\n');
                x2init = Inf;
                for j=1:nIter
                    u = B.merge();
                    x2 = mean(B.aFit(u,aixiyLambda,aizLambda));
                    fprintf(fid,'      %d of %d: x2 = %g\n',j,nIter,x2);
                    if x2init - x2 < x2tol
                        fprintf(fid,'\b <-- reached tolerance\n');
                        break;
                    else
                        x2init = x2;
                    end
                end
            case 'b' % refine b
                fprintf(fid,'    refining scale vs. frame number\n');
                x2init = Inf;
                for j=1:nIter
                    u = B.merge();
                    x2 = mean(B.bFit(u,bizLambda));
                    fprintf(fid,'      %d of %d: x2 = %g\n',j,nIter,x2);
                    if x2init - x2 < x2tol
                        fprintf(fid,'\b <-- reached tolerance\n');
                        break;
                    else
                        x2init = x2;
                    end
                end
            case 'c' % refine c
                fprintf(fid,'    refining constant offset vs. frame number and resolution\n');
                x2init = Inf;
                for j=1:nIter
                    u = B.merge();
                    x2 = mean(B.cFit(u,csLambda,cizLambda,cLambda));
                    
                    if ~isempty(options.cMin)
                        for n=1:length(B)
                            B(n).Model.c(~(B(n).Model.c>=options.cMin)) = options.cMin;
                            B(n).Model.c = B(n).Model.c - min(B(n).Model.cVal) + options.cMin;
                        end
                    end
                    fprintf(fid,'      %d of %d: x2 = %g\n',j,nIter,x2);
                    if x2init - x2 < x2tol
                        fprintf(fid,'\b <-- reached tolerance\n');
                        break;
                    else
                        x2init = x2;
                    end
                end
            case 'd' % refine d
                fprintf(fid,'    refining scale vs. detector panel\n');
                x2init = Inf;
                for j=1:nIter
                    u = B.merge();
                    x2 = B.dFit(u,dipLambda);
                    fprintf(fid,'      %d of %d: x2 = %g\n',j,nIter,x2);
                    if x2init - x2 < x2tol
                        fprintf(fid,'\b <-- reached tolerance\n');
                        break;
                    else
                        x2init = x2;
                    end
                end
            case {'cb','bc'} % refine c and b
                fprintf(fid,'    refining constant offset vs. frame number and resolution\n');
                fprintf(fid,'    alternating with scale vs. frame number\n');
                x2init = Inf;
                for j=1:nIter
                    u = B.merge();
                    x2c = B.cFit(u,csLambda,cizLambda,cLambda);
                    if ~isempty(options.cMin)
                        for n=1:length(B)
                            B(n).Model.c(~(B(n).Model.c>=options.cMin)) = options.cMin;
                            B(n).Model.c = B(n).Model.c - min(B(n).Model.cVal) + options.cMin;
                        end
                    end
                    x2b = B.bFit(u,bizLambda);
                    %x2 = mean([x2c,x2b]);
                    x2 = mean(x2b); % exclude x2c
                    fprintf(fid,'      %d of %d: x2 = %g\n',j,nIter,x2);
                    if x2init - x2 < x2tol
                        fprintf(fid,'\b <-- reached tolerance\n');
                        break;
                    else
                        x2init = x2;
                    end
                end
            case {'o','outlier','outliers'} % remove outliers
                fprintf(fid,'    removing outliers\n');
                sigma_cutoff_0 = x2tol; % passed along this variable
                u = B.merge();
                for n=1:length(B)
                    resid = (B(n).I - B(n).predict(u))./B(n).sigma;
                    % scale sigma cutoff according to actual
                    % statistics
                    sigma_cutoff = sigma_cutoff_0*sqrt(mean(resid.^2));
                    outliers = find(abs(resid)>sigma_cutoff);
                    B(n).sigma(outliers) = Inf; % this effectively removes these data points
                    fprintf(fid,'      %d outliers removed from batch %d\n',length(outliers),n);
                end
            otherwise
                error('did not recognize refinement parameter %s\n',refParam);
        end
        
    end
    
    fprintf(fid,'  done refining scales\n');
    for j=1:nBatches
        ScalingModel(j) = B(j).Model.ScalingModel;
    end
    
    BatchProcess.saveToMatFile(options.matOut,...
            'options',options,...
            'ScalingModel',ScalingModel);
        
    BatchProcess.finish; % done!
    
catch errorMessage
    BatchProcess.stop(errorMessage);
end

tf = BatchProcess.hasCompleted;
errorMessage = BatchProcess.errorMessage;

end