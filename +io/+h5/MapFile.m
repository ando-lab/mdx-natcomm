classdef MapFile < io.h5.FileInterface
    %MAP functions to read and write gridded data in H5 files
    %
    
    methods
        function obj = MapFile(varargin)
            %MAP Construct an instance of this class
            obj@io.h5.FileInterface(varargin{:});
        end
        
        function init(obj,Basis,SpaceGroup)
            init@io.h5.FileInterface(obj);
            obj.writeatt('/','creation_date',datestr(now),'mdx_class',class(obj));
            obj.write_crystal(Basis,SpaceGroup);
        end
        
        function new_grid(obj,loc,Grid,varargin)
            obj.write_grid(loc,Grid)
            if nargin > 3 && ~isempty(varargin)
                obj.writeatt(loc,varargin{:});
            end
        end
        
        function new_dataset(obj,loc,data,varargin)
            Grid = obj.read_grid(fileparts(loc));
            sz = Grid.N;
            
            opts = struct('Datatype',class(data),'ChunkSize',sz,'Deflate',1);
            for j=2:2:numel(varargin)
                opts.(varargin{j-1}) = varargin{j};
            end
            kvopts = [fieldnames(opts),struct2cell(opts)]';
            
            obj.create(loc,sz,kvopts{:});
            
            if ~isempty(data)
                obj.write(loc,data);
            end
        end
        
        function [data,Grid,Basis,SpaceGroup] = import(obj,loc,mode)
            
            if nargin < 3
                mode = 'f'; % return LG as latt.PeriodicGrid (fractional coordinates)
            end
            
            switch lower(mode)
                case {'r','reciprocal'} % reciprocal space
                    Grid = obj.read_grid(fileparts(loc));
                    [Basis,SpaceGroup] = obj.read_crystal();
                    Basis = Basis.invert;
                case {'d','direct'} % direct space
                    Grid = obj.read_grid(fileparts(loc));
                    [Basis,SpaceGroup] = obj.read_crystal();
                case {'f','fractional'} % fractional coordinates
                    Grid = obj.read_grid(fileparts(loc));
                case '' % just read the data
                otherwise
                    error('mode argument not recognized');
            end
            
            data = obj.read(loc);
        end
    end
    
    methods(Hidden)
        
        function write_crystal(obj,Basis,SpaceGroup)
            
            loc = '/crystal';
            
            xtal = struct(...
                'spaceGroupNumber',uint32(SpaceGroup.number),...
                'a',single(Basis.a),...
                'b',single(Basis.b),...
                'c',single(Basis.c),...
                'alpha',single(Basis.alpha),...
                'beta',single(Basis.beta),...
                'gamma',single(Basis.gamma));
                        
            if ~obj.group_exists(loc)
                obj.addgroup(loc);
            end
            
            obj.writeatt(loc,xtal);
            
        end
        
        function [Basis,SpaceGroup] = read_crystal(obj)
            
            loc = '/crystal';
            
            xtal = num2double(obj.readatt(loc));
            
            Basis = latt.Basis(xtal.a,xtal.b,xtal.c,xtal.alpha,xtal.beta,xtal.gamma);
            SpaceGroup = symm.SpaceGroup(xtal.spaceGroupNumber);
            
        end
        
        function Grid = read_grid(obj,loc)
            att = num2double(obj.readatt(loc));
            
            sz = att.grid_size;
            ori = att.grid_ori;
            p = att.grid_delta.*att.grid_size;
            
            Grid = latt.PeriodicGrid(sz(:)', ori(:)', p(:)');
        end
        
        function write_grid(obj,loc,Grid)
            
            g = struct(...
                'grid_size',uint32(Grid.N)',...
                'grid_ori',Grid.ori',...
                'grid_delta',Grid.delta');

            if ~obj.group_exists(loc)
                obj.addgroup(loc);
            end
            
            obj.writeatt(loc,g);
            
        end
    end
end

        
%         function export(obj,loc,mode,data,Grid,Basis,SpaceGroup,varargin)
%             
%             assert(obj.file_exists(),'h5 file does not exist. run init() method to create');
%             
%             if isempty(mode)
%                 mode = 'a'; % append
%             end
% 
%             parentloc = fileparts(loc);
%             
%             switch lower(mode)
%                 case {'a','append'}
%                     assert(obj.group_exists(parentloc))
%                     assert(nargin == 4)
%                     
%                 case {'r','reciprocal'}
%                     assert(nargin == 7)
%                     
%                     obj.write_crystal(Basis.invert,SpaceGroup);
%                     obj.write_grid(parentloc,Grid);
%                     
%                 case {'d','direct'}
%                     assert(nargin == 7)
%                     
%                     obj.write_crystal(Basis,SpaceGroup);
%                     obj.write_grid(parentloc,Grid);
%                     
%                 case {'f','fractional'}
%                     assert(nargin == 5)
%                     
%                     obj.write_grid(parentloc,Grid);
%                     
%                 otherwise
%                     error('mode type not recognized');
%             end
%             
%             obj.create(loc,size(data),'Datatype',class(data),varargin{:});
%             obj.write(loc,data);
% 
% 
%         end



function S = num2double(S)
% convert numeric fields to double precision

fn = fieldnames(S);
for j = 1:numel(fn)
    if isnumeric(S.(fn{j}))
        S.(fn{j}) = double(S.(fn{j}));
    end
end

end
