classdef FileInterface
    %FILEINTERFACE convenient read/write access to h5 files
    
    properties(SetAccess = private)
        fileName
    end
    
    properties
        verbose = false
    end
    
    methods
        function obj = FileInterface(fileName)
            %FILEINTERFACE 
            obj.fileName = fileName;
            
            % TODO: init function?
            %
            % check if file exists
            % check permissions (read / write?)
        end
        
        function tf = init(obj)
            % create a new h5 file if it does not exist
            
            if obj.file_exists() % if filename exists, do nothing
                warning('File %s already exists, so nothing was done. If you want to reinitialize the file, delete it and run FileInterface.init() again.',obj.fileName);
                tf = false;
            else
                fid = H5F.create(obj.fileName,'H5F_ACC_EXCL','H5P_DEFAULT','H5P_DEFAULT');
                H5F.close(fid);
                tf = true;
            end
        end
        
        function tf = file_exists(obj)
            tf = isfile(obj.fileName);
        end
        
        function tf = group_exists(obj,loc)
            assert(obj.file_exists(),'h5 file does not exist');
            fid = H5F.open(obj.fileName,'H5F_ACC_RDONLY','H5P_DEFAULT');
            
            try    
                gid = H5G.open(fid,loc);
                H5G.close(gid);
                tf = true;
            catch
                tf = false;
            end
            
            H5F.close(fid);
            
        end
        
        function create(obj,loc,varargin) % create a new dataset
            if obj.verbose
                fprintf(1,'Creating dataset in %s\n\tLocation: %s\n\tOptions: %s\n',...
                    obj.fileName,...
                    loc,...
                    varargin2string(varargin{:}));
            end
            h5create(obj.fileName,loc,varargin{:});
        end
        
        function write(obj,loc,data,varargin)
            if obj.verbose
                fprintf(1,'Writing data to %s\n\tLocation: %s\n\tOptions: %s\n',...
                    obj.fileName,...
                    loc,...
                    varargin2string(varargin{:}));
            end
            
            h5write(obj.fileName,loc,data,varargin{:});
            
            if obj.verbose, fprintf(1,'Done.\n'); end
        end
        
        function data = read(obj,loc,varargin)
            if obj.verbose
                fprintf(1,'Reading data from %s\n\tLocation: %s\n\tOptions: %s\n',...
                    obj.fileName,...
                    loc,...
                    varargin2string(varargin{:}));
            end
            
            data = h5read(obj.fileName,loc,varargin{:});
            
            if obj.verbose, fprintf(1,'Done.\n'); end
        end
        
        function a = readatt(obj,loc,att)
            %READATT read all attributes as a struct
            
            if nargin == 3 % read one specififfic attribute
                a = h5readatt(obj.fileName,loc,att);
                
            else % read all attributes as a struct
                att_names = obj.listatt(loc);
                a = struct();
                for j=1:numel(att_names)
                    a.(att_names{j}) = h5readatt(obj.fileName,loc,att_names{j});
                end
                
%                 hinfo = h5info(obj.fileName,loc);
%                 if ~isempty(hinfo.Attributes)
%                     a = cell2struct({hinfo.Attributes.Value}',{hinfo.Attributes.Name});
%                 else
%                     a = struct();
%                 end
            end
        end
        
        function a = listatt(obj,loc)
            a = list_attributes(obj.fileName,loc);
        end
        
        
        function writeatt(obj,loc,varargin)
            % varargin: struct or name/value pairs
            
            if nargin==3 && isstruct(varargin{1})
                attNames = fieldnames(varargin{1});
                attVals = struct2cell(varargin{1});
            else
                attNames = varargin(1:2:(end-1));
                attVals = varargin(2:2:end);
            end

            for j = 1:numel(attNames)
                h5writeatt(obj.fileName,loc,attNames{j},attVals{j});
            end
        end
        
        function addgroup(obj,loc)

            if obj.group_exists(loc)
                warning('group already exists, so none was added');
                return;
            end
            
            % see:
            % https://support.hdfgroup.org/ftp/HDF5/examples/examples-by-api/matlab/HDF5_M_Examples/h5ex_t_cmpd.m
            
            fid = H5F.open(obj.fileName,'H5F_ACC_RDWR','H5P_DEFAULT');
            gid = H5G.create(fid,loc,'H5P_DEFAULT','H5P_DEFAULT','H5P_DEFAULT');

            % clean up
            H5G.close(gid);
            H5F.close(fid);
            
        end
        
        function disp(obj)
            if ~isfile(obj.fileName)
                builtin('disp',obj);
                return;
            end
            
            hinfo = h5info(obj.fileName,'/');
            
            fprintf(1,'%s:\n',obj.fileName);
            print_h5_tree(hinfo);
            
            function print_h5_tree(s,depth)
                if nargin == 1, depth = 0; end
                printGroup(depth,s);
                for j=1:numel(s.Attributes)
                    printAttribute(depth + 1,s.Attributes(j));
                end
                for j=1:numel(s.Datasets)
                    printDataset(depth + 1,s.Datasets(j));
                end
                for j=1:numel(s.Groups)
                    print_h5_tree(s.Groups(j),depth + 1);
                end
            end
            
            function printAttribute(depth,attInfo)
                name = attInfo.Name;
                val = attInfo.Value;
                if isnumeric(val)
                    if numel(val) == 1
                        val = num2str(val);
                    else
                        sz = strjoin(arrayfun(@num2str,size(val),'Uni',0),'x');
                        type = attInfo.Datatype.Class;
                        val = sprintf('[%s %s]',sz,type);
                    end
                end
                if ~ischar(val)
                    val = '?';
                end
                fprintf(1,'%s@%s: %s\n',repmat('   ',[1,depth]),name,val);
            end
            
            function printDataset(depth,dsetInfo)
                name = dsetInfo.Name;
                sz = strjoin(arrayfun(@num2str,dsetInfo.Dataspace.Size,'Uni',0),'x');
                type = dsetInfo.Datatype.Class;
                
                fprintf(1,'%s%s: [%s %s]\n',repmat('   ',[1,depth]),name,sz,type);
                for j=1:numel(dsetInfo.Attributes)
                    printAttribute(depth + 1,dsetInfo.Attributes(j));
                end
            end
            
            function printGroup(depth,groupInfo)
                fprintf(1,'%s%s\n',repmat('   ',[1,depth]),groupInfo.Name);
            end
            
        end
        

    end
end


function s = varargin2string(varargin)

if isempty(varargin)
    s = '{}';
else
    s = strtrim(regexprep(formattedDisplayText(varargin),{'}\s+{'},{', '}));
end

end



function attr_names = list_attributes(h5filename,loc)

attr_names = {};

try

    fid = H5F.open(h5filename);
    
    try
        gid = H5G.open(fid,loc);
        H5A.iterate(gid,[],@iterator_func);
    catch EM1
        H5G.close(gid);
        rethrow(EM)
    end
catch EM2
    H5F.close(fid);
    rethrow(EM)
end

H5G.close(gid);
H5F.close(fid);

function status = iterator_func(~,attr_name)
    status = 0;
    attr_names = [attr_names,{attr_name}];
end

end