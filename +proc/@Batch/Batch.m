classdef Batch < handle
    % Scripts for batch-processing diffuse scattering data
    %
    % All scripts read and write files, and have a common interface:
    %
    %   [tf,errorMessage] = proc.Batch.scriptName(<name value pairs>)
    %
    % Scripts fail quietly (rather than throw an error, tf is set to false
    % and the caught error message appears as the second output argument.
    %
    % The error message can be retrieved by calling rethrow(errorMessage)
    %
    % proc.Batch Methods:
    %     autorun - chain batch scripts together
    %    xds2geom - harvest geometry information from xds log files
    %    cbf2geom - harvest geometry information from cbf file headers
    %       merge - apply scales and merge equivalent reflections
    %      filter - use Poisson statistics to classify voxels
    %        grid - design a grid from scratch (instead of using filter)
    %   integrate - sum intensities in 3D
    % reintegrate - integrate, correct, apply scales, merge, and export
    %     correct - calculate geometry and background corrections
    %      export - write corrected intensities to a table
    %     combine - put multiple batches into one table (for scaling)
    %       scale - refine a scale model against multiple batches of data
    %       merge - apply a scale model to a table of intensities
    %     rescale - absolute intensity scaling
    
    properties(SetAccess = immutable)
        name = 'proc.Batch'
        options
        workingDirectory = './'   % may be changed by options
        logOut = ''               % keeps a log if not empty
    end
    properties(SetAccess = private)
        startTime
        endTime
        hasCompleted = false
        errorMessage = []
        logFileWriter = util.Logger.empty();
    end
    methods(Access = private) % may be used only by static methods
        function obj = Batch(name,defaultOpts,varargin)
            obj.name = name;
            obj.options = parseVarargin(defaultOpts,varargin{:});

            % assign any options pertaining to proc.Batch
            optFields = intersect(fieldnames(obj.options),...
                {'workingDirectory','logOut'});
            for j=1:length(optFields)
                obj.(optFields{j}) = obj.options.(optFields{j});
            end

            if ~isempty(obj.logOut)
                logFileName = fullfile(obj.workingDirectory,obj.logOut);
                path = fileparts(logFileName);
                if ~isdir(path)
                    mkdir(path);
                end
                obj.logFileWriter = util.Logger(logFileName);
            end
        end

        function start(obj)
            obj.startTime = datestr(now);
            if ~isempty(obj.logFileWriter)
                obj.logFileWriter.start();
            end
            fprintf(1,'\n%s\ndate: %s\n',obj.name,obj.startTime);
            % print options
            fprintf(1,'%s\n',opt2str(obj.options));
        end

        function finish(obj)
            obj.endTime = datestr(now);
            fprintf(1,'\n%s\ndate: %s\n',obj.name,obj.endTime);
            fprintf(1,'Done!\n');
            obj.hasCompleted = true;
            
            if ~isempty(obj.logFileWriter)
                delete(obj.logFileWriter);
            end
        end

        function stop(obj,errorMessage)
            obj.endTime = datestr(now);
            obj.errorMessage = errorMessage;
            warning('%s ended with the error message:\n\t%s',...
                obj.name,errorMessage.message);
            if ~isempty(obj.logFileWriter)
                delete(obj.logFileWriter);
            end
        end

        function delete(obj)
            if ~isempty(obj.logFileWriter)
                delete(obj.logFileWriter);
            end
        end

        function fullFilePath = prepareOutputFile(obj,fileName,varargin)
            fileName = fullfile(obj.workingDirectory,fileName);
            fullFilePath = prepareOutputFile(fileName,varargin{:});
        end

        function fileName = saveToMatFile(obj,fileName,varargin)
            if length(varargin)==1 && isstruct(varargin{1})
                dataToSave = varargin{1};
            else
                dataToSave = struct();
                isFlag = strncmp(varargin,'-',1); % pass through '-append',etc
                nvpairs = varargin(~isFlag);
                for j=1:2:length(nvpairs)
                    dataToSave.(nvpairs{j}) = nvpairs{j+1};
                end
            end
            fprintf(1,'  checking output file %s\n',fileName);
            fileName = fullfile(obj.workingDirectory,fileName);
            fileName = prepareOutputFile(fileName,'mat');

            fn = fieldnames(dataToSave);
            fprintf(1,'  saving workspace variables to file: %s\n',fileName);
            if any(isFlag)
                fprintf(1,'  with options: %s\n',strjoin(varargin(isFlag)));
            end
            fprintf(1,'    %s\n',fn{:});

            save(fileName,'-struct','dataToSave',varargin{isFlag});
        end

        function fullPath = checkInputDirectory(obj,dirPath)
            fprintf(1,'  checking input directory %s\n',dirPath);
            dirPath = fullfile(obj.workingDirectory,dirPath);
            fullPath = checkInputFile(dirPath,'dir');
        end

        function [varargout] = readFromMatFile(obj,fileName,varargin)
            fprintf(1,'  checking input file %s\n',fileName);
            fileName = fullfile(obj.workingDirectory,fileName);
            [fileName,fileVars] = checkInputFile(fileName,'mat');
           
            if ~isempty(varargin) % read specific variables
                varargout = cell(length(varargin),1);

                fprintf(1,'  reading workspace variables from file %s\n',fileName);
                varsToRead = intersect(varargin,fileVars);
                fileContents = load(fileName,varsToRead{:});
                for j=1:length(varargin)
                    if ismember(varargin{j},fileVars)
                        fprintf(1,'    %s\n',varargin{j});
                        varargout{j} = fileContents.(varargin{j});
                    else
                        fprintf(1,'    %s (not found. assigned as [])\n',varargin{j});
                        varargout{j} = [];
                    end
                end
            else % load all variables as a struct (lazy way)
                varargout{1} = load(fileName);
            end
        end
    end

    methods(Static=true)
        
        % regular scripts: 
        [bin,corr] = repartitionScript(fid,options,Grid,bin,corr);
        
        diffuseTable = exportScaledScript(fid,options,Grid,bin,corr,ScalingModel,Detector)
        
        % batch scripts:

        [tf, errorMessage] = autorun(varargin)

        [tf, errorMessage] = cbf2geom(varargin)

        [tf, errorMessage] = scale(varargin)

        [tf, errorMessage] = merge(varargin)

        [tf, errorMessage] = combine(varargin)

        [tf, errorMessage] = grid(varargin)

        [tf, errorMessage] = filter(varargin)

        [tf, errorMessage] = xds2geom(varargin)

        [tf, errorMessage] = integrate(varargin)

        [tf, errorMessage] = reintegrate(varargin)
        
        [tf, errorMessage] = correct(varargin)

        [tf, errorMessage] = export(varargin)
        
        [tf, errorMessage] = rescale(varargin)
        
    end
end

% some private methods

function [fileOut,varargout] = checkInputFile(fileIn,type)
% if mat file, varargout{1} will contain a list of variables
varargout = {};
[path,name,ext] = fileparts(fileIn);
switch lower(type)
    case {'dir','directory'}
        if ~isempty(name) && ~isempty(ext)
            error('input directory %s is not a directory',fileIn);
        elseif ~isdir(fileIn)
            error('input directory %s does not exist or is not in the path',fileIn);
        else
            fileOut = fileIn;
        end
    case {'mat','matfile'}
        if isempty(ext)
            ext = '.mat';
        elseif strcmp(ext,'.mat')
            % fine, do nothing
        else
            error('expected .mat file extension on input file %s',fileIn);
        end
        fileOut = fullfile(path,[name ext]);
        if isempty(dir(fileOut))
            error('input file %s does not exist or is not in the path',fileOut);
        end
        vars = whos('-file',fileOut);
        varargout{1} = {vars.name};
    otherwise
        error('did not recognize input type');
end

end

function opt = parseVarargin(opt,varargin)
if length(varargin)==1 && isstruct(varargin{1})
    fn = fieldnames(varargin{1});
    for j=1:length(fn)
        if isfield(opt,fn{j})
            opt.(fn{j}) = varargin{1}.(fn{j});
        else
            warning('%s is not a valid option',fn{j});
        end
    end
else % name,value pairs
    for j=1:2:length(varargin)
        if isfield(opt,varargin{j})
            opt.(varargin{j}) = varargin{j+1};
        else
            warning('%s is not a valid option',varargin{j});
        end
    end
end
end

function matOut = prepareOutputFile(matIn,fileExt)
if nargin < 2 || isempty(fileExt)
    fileExt = 'mat';
end
[path,name,ext] = fileparts(matIn);
if isempty(ext)
    ext = ['.' fileExt];%
elseif strcmp(ext,['.' fileExt])
    % fine - do nothing
else
    error('expected %s file extension',fileExt);
end
if ~isdir(path)
    mkdir(path);
end
matOut = fullfile(path,[name ext]);
if ~isempty(dir(matOut))
    fprintf(1,'\nWARNING: Output file already exists and will be overwritten.\n\n');
end
end
