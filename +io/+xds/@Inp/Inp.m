classdef Inp
    methods (Static)
        function xds_inp = read(fileName)
            if nargin==0 || isempty(fileName)
                fileName = 'XDS.INP';
            elseif isdir(fileName)
                fileName = fullfile(fileName,'XDS.INP');
            end
            xds_inp = read_xds_inp(fileName);
        end
        
        function xds_inp = new()
            xds_inp = read_xds_inp();
        end
        
        function write(xds_inp,fileName)
            if nargin==0 || isempty(fileName)
                fileName = 'XDS.INP';
            elseif isdir(fileName)
                fileName = fullfile(fileName,'XDS.INP');
            end
            write_xds_inp(xds_inp,fileName);
        end
        
        function m = getMask(xds_inp)
            nx = xds_inp.nx;
            ny = xds_inp.ny;
            m = true(nx,ny);
            if ~isempty(xds_inp.untrusted_rectangle)
                m = m & rectangle2mask(xds_inp.untrusted_rectangle,nx,ny);
            end
            if ~isempty(xds_inp.untrusted_ellipse)
                m = m & ellipse2mask(xds_inp.untrusted_ellipse,nx,ny);
            end
            if ~isempty(xds_inp.untrusted_quadrilateral)
                warning(['untrusted_quadrilateral has not been ',...
                    'implemented in this version, and will not be ',...
                    'included in the mask']);
                
            end
        end
        
        function print(xds_inp)
            if nargin==0 || isempty(xds_inp)
                error('empty xds input structure');
            end
            write_xds_inp(xds_inp);
        end
        
        function draw(xds_inp)
            if nargin==0 || isempty(xds_inp)
                error('empty xds input structure');
            end
            draw_detector(xds_inp);            
        end
        
    end
end