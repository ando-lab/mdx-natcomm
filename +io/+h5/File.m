classdef File
    %io.h5.File - convenient interface for reading and writing h5 files
    %
    %  Should work with any h5 file, with caveats:
    %
    %  The root level attributes mdx_version and mdx_creation_date are
    %  written by default when creating new files.
    %
    %  When opening an existing file, mdx_version is checked against the
    %  version property of the FileInterface class. If they don't match (or
    %  if mdx_version does not exist) then warnings are issued.
    
    properties(Constant)
        version = 1; % written as root Attribute: mdx_version 
    end
    
    properties(SetAccess = immutable)
        h5filename
    end
    
    methods
        function obj = File(h5filename)
            obj.h5filename = h5filename;
            if ~isfile(obj.h5filename)
                obj.init();
            else % file exists, check version
                obj.versionCheck(); % issue warnings if version differs
            end
        end
        
        function versionCheck(obj)
            if ~isfile(obj.h5filename)
                error('file does not exist');
            end
            try
                attr = obj.readatt('/',{'mdx_version'});
                if attr.mdx_version ~= obj.version
                    warning('file version (v%g) does not match software version (v%g)... proceed with caution',attr.mdx_version,obj.version);
                end
            catch
                warning('mdx_version not found. h5 file might not work with FileInterface');
            end
        end
        
        function init(obj,varargin)
            if isfile(obj.h5filename)
                warning('%s exists and will be overwritten\n',obj.h5filename);
                delete(obj.h5filename);
            end
            
            h5loc = '/';
            
            fid = H5F.create(obj.h5filename,'H5F_ACC_EXCL','H5P_DEFAULT','H5P_DEFAULT');
                        
            % clean up
            H5F.close(fid);
            
            % write creation date attribute and version
            h5writeatt(obj.h5filename,h5loc,'mdx_creation_date',datestr(now));
            h5writeatt(obj.h5filename,h5loc,'mdx_version',obj.version);
            
            % write optional attributes
            obj.writeatt(h5loc,varargin{:});
            
        end
        
        function [attr,otherattr] = readatt(obj,h5loc,mandatory_fields)
            
            info = h5info(obj.h5filename,h5loc);
            attr = cell2struct({info.Attributes.Value}',{info.Attributes.Name});
            
            if nargin > 2 && ~isempty(mandatory_fields)
                assert(all(ismember(mandatory_fields,fieldnames(attr))),...
                    'mandatory fields not found');
                if nargout == 2
                    otherattr = rmfield(attr,mandatory_fields);
                end
                other_fields = setdiff(fieldnames(attr),mandatory_fields,'stable');
                attr = rmfield(attr,other_fields);
            end
            
        end
        
        function writeatt(obj,h5loc,varargin)
            % varargin is either:
            % cell array (key, value pairs)
            % or a scalar struct
            
            if isempty(varargin)
                % do nothing
                return;
            elseif isstruct(varargin{1})
                assert(numel(varargin) == 1 && numel(varargin{1}) == 1,...
                    'expected scalar struct');
                fn = fieldnames(varargin{1});
                for j=1:numel(fn)
                    h5writeatt(obj.h5filename,h5loc,fn{j},varargin{1}.(fn{j}));
                end
            else % key, value pairs
                for j=2:2:numel(varargin)
                    h5writeatt(obj.h5filename,h5loc,varargin{j-1},varargin{j});
                end
            end
            
        end
        
        function tf = groupexists(obj,groupname,fid)
            
            closeFileWhenFinished = false;
            
            if nargin < 3 || isempty(fid)
                fid = H5F.open(obj.h5filename,'H5F_ACC_RDONLY','H5P_DEFAULT');
                closeFileWhenFinished = true;
            end
            
            try
                gid = H5G.open(fid,groupname);
                tf = true;
                H5G.close(gid);
            catch
                tf = false;
            end
            
            if closeFileWhenFinished, H5F.close(fid); end
            
        end
        
        function tf = addgroup(obj,groupname,varargin)
            % adapted from:
            % https://support.hdfgroup.org/ftp/HDF5/examples/examples-by-api/matlab/HDF5_M_Examples/h5ex_t_cmpd.m
            
            tf = false;
            
            groupPath = strsplit(groupname,'/');
            assert(isempty(groupPath{1}),'group path must begin with root (/).');
            
            % open h5 file
            fid = H5F.open(obj.h5filename,'H5F_ACC_RDWR','H5P_DEFAULT');
            
            for j=2:numel(groupPath)
                thisGroup = strjoin(groupPath(1:j),'/');
                if ~obj.groupexists(thisGroup,fid)
                    % if group requested, add it and close
                    gid = H5G.create(fid,thisGroup,'H5P_DEFAULT','H5P_DEFAULT','H5P_DEFAULT');
                    % clean up
                    H5G.close(gid);
                    tf = true; % at least one group was added
                end
            end

            % write attributes if requested
            obj.writeatt(groupname,varargin{:});
            
            H5F.close(fid);
            
        end
        
        function create(obj,dsname,grid_size,varargin)
            % wrapper for h5fcreate
            %
            % varargin is in the form of keyword, value pairs.

            h5create(obj.h5filename,dsname,grid_size,varargin{:});
            
        end
        
        function write(obj,dsname,data,varargin)
            % wrapper for h5write
            %
            % varargin is in the form of keyword, value pairs.
            
            h5write(obj.h5filename,dsname,data,varargin{:});
        end
        
        function val = read(obj,dsname,varargin)
            % wrapper for h5read
            %
            % varargin is in the form of keyword, value pairs.
            
            val = h5read(obj.h5filename,dsname,varargin{:});
        end
        
        function disp(obj)
            % replaces h5disp
            
            %h5disp(obj.h5filename,'/','min');
            fprintf(1,'%s:\n',obj.h5filename);
            fileMeta = h5info(obj.h5filename);
            printTree(fileMeta,1);
            
            function printTree(f,l)
                indentation = repmat('  ',[1,l-1]);
                
                if isfield(f,'Groups') % folder
                    [~,endName] = fileparts(f.Name);
                    fprintf(1,'%s%s/\n',indentation,endName);
                else % a dataset
                    sz = strjoin(arrayfun(@num2str,f.Dataspace.Size,'Uni',0),'x');
                    type = f.Datatype.Class;
                    fprintf(1,'%s%s: [%s %s]\n',indentation,f.Name,sz,type);
                end
                
                indentation = repmat('  ',[1,l]);
                
                for j=1:numel(f.Attributes)
                    fprintf(1,'%s%s%s: %s\n',indentation,'@',f.Attributes(j).Name,formatAttributeValue(f.Attributes(j)));
                end
                if isfield(f,'Datasets')
                    for j=1:numel(f.Datasets)
                        printTree(f.Datasets(j),l+1);
                    end
                end
                if isfield(f,'Groups')
                    for j=1:numel(f.Groups)
                        printTree(f.Groups(j),l+1);
                    end
                end
                
                function str = formatAttributeValue(att)
                    maxlen = 60;
                    switch att.Datatype.Class
                        case 'H5T_STRING'
                            str = att.Value;
                        case {'H5T_INTEGER','H5T_FLOAT'}
                            str = strjoin(arrayfun(@num2str,att.Value,'Uni',0),', ');
                        otherwise
                            str = '[other]';
                    end
                    if numel(str) > maxlen
                        str = [str(1:(maxlen-3)), '...'];
                    end
                end
                
            end
            
        end
    end
end
