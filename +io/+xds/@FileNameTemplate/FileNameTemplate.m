classdef FileNameTemplate
    properties(Dependent)
        xds
        format
        wildcard
    end
    
    properties(Access = protected)
        fileTemplate
        nDigits
        isZeroPadded
    end
    
    methods
        
        function obj = FileNameTemplate(nameTemplate)
    
            if isxds(nameTemplate)
                % it's XDS-formatted string (question marks)
                [fileTemplate,nDigits,isZeroPadded] = xds2template(nameTemplate);
            elseif isformat(nameTemplate)
                % it's a format string (e.g. for sprintf)
                [fileTemplate,nDigits,isZeroPadded] = format2template(nameTemplate);
            elseif isexample(nameTemplate)
                % it's an example file name
                [fileTemplate,nDigits,isZeroPadded] = example2template(nameTemplate);
            else
                error('did not recognize file name template');
            end
            
            obj.fileTemplate = fileTemplate;
            obj.nDigits = nDigits;
            obj.isZeroPadded = isZeroPadded;

        end
        
        function value = get.xds(obj)
            value = template2xds(obj.fileTemplate,obj.nDigits);
        end
        
        function value = get.format(obj)
            value = template2format(obj.fileTemplate,obj.nDigits,obj.isZeroPadded);
        end
        
        function value = get.wildcard(obj)
            value = template2wildcard(obj.fileTemplate);
        end
        
        function fileName = getFileName(obj,varargin)
            assert(length(varargin)==length(obj.nDigits),...
                'number of arguments was different from expected');
            
            fileName = sprintf(obj.format,varargin{:});

        end
    end
end