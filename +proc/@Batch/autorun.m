function [tf, errorMessage] = autorun(varargin)
% autorun - run several batch processes in sequence
%
% this handles case where a struct array is passed as the only argument.
% autorun executes separately on each instance of the struct array, and
% returns a logical array (tf) and cell array (errorMessage) with one entry
% for each element of the struct array
if length(varargin)==1 && isstruct(varargin{1}) && numel(varargin{1})>1
    nBatches = length(varargin{1});
    tf = false(nBatches,1);
    errorMessage = cell(nBatches,1);
    for j=1:nBatches
        [tf(j),errorMessage{j}] = autorun1(varargin{1}(j));
    end
else % else, run normally
    [tf, errorMessage] = autorun1(varargin{:});
end
end

function [tf, errorMessage] = autorun1(varargin)

options = struct(...
    'run','',...   % set below to prevent struct array
    'workingDirectory','./',... % for all
    'xdsDir','../xds/',...      % for xds2geom
    'matrixOperator',[],...            % for reindex
    'fileNameTemplate','',...   % for cbf2geom
    'frameRange',[],...
    'startingAngleField','Start_angle',...
    'oscillationRangeField','Angle_increment',...
    'ndiv',[3,3,3],...          % for filter
    'window',2,...
    'maxCount',20,...
    'smax',Inf,...
    'refineGrid',true,...
    'parallel',false,...        % for integrate
    'binMode','coarse',...      %'arrayType','full',...
    'binExcludedVoxels',true,...
    'readBackground',true,...
    'readImageHeaders',true,...
    'exposureTimeField','Exposure_time',...
    'minimumCounts',0,...       % for export
    'minimumPixels',1);

BatchProcess = proc.Batch('proc.Batch.autorun',options,varargin{:});
options = BatchProcess.options; % post-processed options

try 
    BatchProcess.start();
    
    for j=1:length(options.run)
        switch lower(options.run{j})
            case 'xds2geom'
                [tf, errorMessage] = proc.Batch.xds2geom(...
                    getFields(options,'workingDirectory','xdsDir'));
            case 'reindex'
                if isempty(options.matrixOperator)
                    fprintf(1,'\nproc.Batch.reindex (skipping since matrixOperator is empty)\n');
                    continue; % skip
                end
                [tf, errorMessage] = proc.Batch.reindex(...
                    getFields(options,'workingDirectory','matrixOperator'));
            case 'cbf2geom'
                [tf, errorMessage] = proc.Batch.cbf2geom(...
                    getFields(options,'workingDirectory',...
                    'fileNameTemplate','frameRange',...
                    'startingAngleField','oscillationRangeField'));
            case 'filter'
                [tf, errorMessage] = proc.Batch.filter(...
                    getFields(options,'workingDirectory',...
                    'ndiv','window','maxCount','smax','refineGrid','parallel'));
            case 'integrate'
                [tf, errorMessage] = proc.Batch.integrate(...
                    getFields(options,'workingDirectory',...
                    'binMode','parallel','binExcludedVoxels'));
            case 'correct'
                [tf, errorMessage] = proc.Batch.correct(...
                    getFields(options,'workingDirectory',...
                    'binMode','readBackground','readImageHeaders',...
                    'exposureTimeField','parallel','binExcludedVoxels'));
            case 'export'
                [tf, errorMessage] = proc.Batch.export(...
                    getFields(options,'workingDirectory',...
                    'minimumCounts','minimumPixels','binExcludedVoxels'));
            otherwise
                tf = false;
                try % throw a silent error
                    error('did not recognize program %s',options.run{j});
                catch errorMessage
                end
                
        end
        if ~tf, rethrow(errorMessage); end
    end
    BatchProcess.finish;
catch errorMessage
    BatchProcess.stop(errorMessage);
end

tf = BatchProcess.hasCompleted;
errorMessage = BatchProcess.errorMessage;

end

function [ sOut ] = getFields( sIn, varargin )
for j=1:length(varargin)
    sOut.(varargin{j}) = sIn.(varargin{j});
end
end
