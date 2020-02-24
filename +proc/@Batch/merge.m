function [tf, errorMessage] = merge(varargin)

options = struct(...
    'workingDirectory','./',...
    'combineIn','combine.mat',...
    'scaleIn','scale.mat',...
    'matOut','merge.mat',...
    'logOut','merge.log',...
    'mergeBragg',false,... % if true, use braggTable instead of diffuseTable. Turns off constant offset in ScalingModel
    'sigmaCutoff',5);
options.mergeByColumns = {'hasu','kasu','lasu'};

BatchProcess = proc.Batch('proc.Batch.merge',options,varargin{:});
options = BatchProcess.options; % post-processed options

try % START MAIN SCRIPT
        
    BatchProcess.start();

    ScalingModel = BatchProcess.readFromMatFile(options.scaleIn,'ScalingModel');
    
    if options.mergeBragg
        hklTable = BatchProcess.readFromMatFile(...
            options.combineIn,'braggTable');
        for j=1:length(ScalingModel)
            ScalingModel(j).c = []; % remove constant offset
        end
    else
        hklTable = BatchProcess.readFromMatFile(...
            options.combineIn,'diffuseTable');
    end
    
    [hklMerge,ih,isOutlier] = mergeScript(1,options,ScalingModel,hklTable);
    
    BatchProcess.saveToMatFile(options.matOut,...
            'options',options,...
            'hklMerge',hklMerge,...
            'ih',ih,...
            'isOutlier',isOutlier);
            
    BatchProcess.finish; % done!
    
catch errorMessage
    BatchProcess.stop(errorMessage);
end

tf = BatchProcess.hasCompleted;
errorMessage = BatchProcess.errorMessage;

end