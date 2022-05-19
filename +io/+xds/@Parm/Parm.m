classdef Parm
    methods(Static)
        function xparm = read(fileName)
            if nargin==0 || isempty(fileName)
                fileName = 'XPARM.XDS';
            elseif isdir(fileName)
                fileName = fullfile(fileName,'XPARM.XDS');
            end
            xparm = load_xparm(fileName);
        end
        function tf = write(xparm,fileName)
            if nargin==1 || isempty(fileName)
                fileName = 'XPARM.XDS';
            elseif isdir(fileName)
                fileName = fullfile(fileName,'XPARM.XDS');
            end
            tf = write_xparm(xparm,fileName);
        end
    end
end