classdef ImageSeries < util.propertyValueConstructor
    % ImageSeries is a convenient wrapper function for reading in image
    % data, and can be extended to deal with changing formats
    properties
        fileNameTemplate
        fileReadFunction = @io.DectrisCBF.read
        % fileReadFunction must accept a file name as argument, and return
        % two values: the image data and the image header
        headerReadFunction = @io.PilatusHeader.read
        frameRange = [] % optional. if not empty, will check that frames are within range
        excludedFrames = [] % optional
    end
    properties(Dependent = true)
        frames
        files
        datastore % export to datastore object, which has similar functionality but supports mapreduce
    end
        
    methods
        function obj = ImageSeries(varargin)
            obj@util.propertyValueConstructor(varargin{:});
        end
        
        function val = get.frames(obj)
            if ~isempty(obj.frameRange) && length(obj.frameRange)==2
            val = setdiff(obj.frameRange(1):obj.frameRange(end),obj.excludedFrames);
            else
                val = [];
            end
        end
        
        function fn = get.files(obj)
            imageNumbers = obj.frames;
            if ~isempty(imageNumbers) && ~isempty(obj.fileNameTemplate)
                    fn = arrayfun(@(n) sprintf(obj.fileNameTemplate,n),...
                    imageNumbers,'uniformOutput',false);
            else
                fn = {};
            end
        end
        
        function ds = get.datastore(obj)
            fileNames = obj.files;
            if ~isempty(fileNames)
                [~,~,ext] = fileparts(obj.fileNameTemplate);
                ds = imageDatastore(fileNames,...
                    'ReadFcn',obj.fileReadFunction,...
                    'FileExtensions',ext,...
                    'Labels',obj.frames);
            else
                error('file names required to create datastore');
            end
        end
        
        function header = readHeader(obj,frameList)
            if nargin==1
                frameList = obj.frames;
            end
            nFrames = length(frameList);
            assert(nFrames>=1,'one or more frames expected');
            
            for j=1:nFrames
                currentFrame = frameList(j);
                thisFileName = sprintf(obj.fileNameTemplate,currentFrame);
                header(j) = obj.headerReadFunction(thisFileName);
            end
        end
        
        function [data,header] = readStack(obj,frameList)
            if nargin==1
                frameList = obj.frames;
            end
            
            nFrames = length(frameList);
            assert(nFrames>=1,'one or more frames expected');
            
            currentFrame = frameList(1);
            thisFileName = sprintf(obj.fileNameTemplate,currentFrame);
            [thisImage,header(1)] = obj.fileReadFunction(thisFileName);
            
            imageSize = size(thisImage);
            imageType = class(thisImage);
            
            data = zeros([imageSize,nFrames],imageType);
            
            data(:,:,1) = thisImage;
            
            for j=2:length(frameList)
                currentFrame = frameList(j);
                thisFileName = sprintf(obj.fileNameTemplate,currentFrame);
                [data(:,:,j),header(j)] = obj.fileReadFunction(thisFileName);
            end

        end
        
        function [data,header] = readSum(obj,frameList)
            if nargin==1
                frameList = obj.frames;
            end
            
            nFrames = length(frameList);
            assert(nFrames>=1,'one or more frames expected');
            
            currentFrame = frameList(1);
            thisFileName = sprintf(obj.fileNameTemplate,currentFrame);
            
            [data,header(1)] = obj.fileReadFunction(thisFileName);
            
            for j=2:length(frameList)
                currentFrame = frameList(j);
                thisFileName = sprintf(obj.fileNameTemplate,currentFrame);
                [thisImage,header(j)] = obj.fileReadFunction(thisFileName);
                data = data + thisImage;
            end
            
        end
        
    end
end
