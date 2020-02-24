classdef DectrisCBF
    properties (SetAccess = protected)
        file
        header
        image
    end
    methods
        function obj = DectrisCBF(fileName)
            %obj.file = io.DataFile(fileName);
            obj.file = fileName;
            obj.header = io.PilatusHeader(fileName);
            obj.image = readImageData(obj,fileName);           
        end
    end
    methods(Static = true)
        function [img,header] = read(fileName)
            imobj = io.DectrisCBF(fileName);
            img = imobj.image;
            header = imobj.header.HeaderContents;
        end
    end
    methods (Access = protected)
        function imageData = readImageData(obj,fileName)
            %disp(obj.header.CIFHeader); % for debugging
            fpos = obj.header.binary_start;
            nbytes = obj.header.CIFHeader.X_Binary_Size;
            nrows = obj.header.CIFHeader.X_Binary_Size_Fastest_Dimension;
            ncols = obj.header.CIFHeader.X_Binary_Size_Second_Dimension;
            
            imageData = read_cbf_image(fileName,fpos,nbytes,nrows,ncols);
            
        end
    end
end