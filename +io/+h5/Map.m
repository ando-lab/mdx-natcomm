classdef Map
    %MAP read and write reciprocal space maps to hdf5 format
    
    properties(SetAccess = private)
        h5file
        mapLocation = '/maps'
        crystalLocation = '/crystal'; % default for mdx-lib
        %mapSpace = 'recip' % choices: direct, recip / reciprocal
    end
    
    methods
        
        function obj = Map(h5file,varargin)
            for j=2:2:nargin
                obj.(varargin{j-1}) = varargin{j};
            end
            
            obj.h5file = h5file;
            
            if ischar(obj.h5file)
                obj.h5file = io.h5.File(h5file);
            else
                obj.h5file = h5file;
            end
        end
        
        function write_crystal(obj,Basis,SpaceGroup)
            
            crystal = struct(...
                'spaceGroupNumber',uint32(SpaceGroup.number),...
                'a',single(Basis.a),...
                'b',single(Basis.b),...
                'c',single(Basis.c),...
                'alpha',single(Basis.alpha),...
                'beta',single(Basis.beta),...
                'gamma',single(basis.gamma));
            
            obj.h5file.addgroup(obj.crystalLocation,crystal);
            
        end
        
        function [Basis,SpaceGroup] = read_crystal(obj)

            assert(obj.h5file.groupexists(obj.crystalLocation));
            
            mandatory_fields = {'spaceGroupNumber','a','b','c','alpha','beta','gamma'};
            
            att = obj.h5file.readatt(obj.crystalLocation,mandatory_fields);

            spaceGroupNumber = double(att.spaceGroupNumber);
            a = double(att.a);
            b = double(att.b);
            c = double(att.c);
            alpha = double(att.alpha);
            beta = double(att.beta);
            gamma = double(att.gamma);  
            
            Basis = latt.Basis(a,b,c,alpha,beta,gamma);
            SpaceGroup = symm.SpaceGroup(spaceGroupNumber);
        end
        
        function write_grid(obj,P)
            
            grid = struct('grid_size',uint32(P.N),...
                'grid_ori',P.ori,...
                'grid_delta',P.delta);
            
            obj.h5file.addgroup(obj.mapLocation,grid);
            
        end
        
        function P = read_grid(obj)
            
            assert(obj.h5file.groupexists(obj.mapLocation));
            
            mandatory_fields = {'grid_size','grid_ori','grid_delta'};
            
            att = obj.h5file.readatt(obj.mapLocation,mandatory_fields);
            
            N = double(att.grid_size(:))';
            ori = double(att.grid_ori(:))';
            delta = double(att.grid_delta(:))';
            
            P = latt.PeriodicGrid(N,ori,N.*delta);
            
        end
        
        function mapdata = read_map(obj,name)
            assert(obj.h5file.groupexists(obj.mapLocation));
            mapdata = obj.h5file.read(fullfile(obj.mapLocation,name));
        end
        
        function info = write_map(obj,name,mapdata)
            
            P = obj.read_grid();
            ndiv = round(1./P.delta);
            arraysize = P.N;
            
            opts = {'Datatype','single',...
                'FillValue',single(NaN),...
                'ChunkSize',ndiv,...
                'Deflate',1};
            
            obj.h5file.create(fullfile(obj.mapLocation,name),arraysize,opts{:});
            
            sz0 = file_size_mb(obj.h5file); % initial file size
            t0 = tic; % start timer
            
            obj.h5file.write(fullfile(obj.mapLocation,name),single(mapdata));

            sz1 = file_size_mb(obj.h5file);
            dt = toc(t0);
             
            info = struct('size_mb',sz1-sz0,'time_s',dt);
            
            function sz = file_size_mb(h5f)
                f = dir(h5f.h5filename);
                sz = f.bytes/2^20;
            end

        end
        
        function [mapdata,P,Basis,SpaceGroup] = read(obj,name)
            mapdata = obj.read_map(name);
            if nargout >= 2
                P = obj.read_grid();
            end
            if nargout >= 3
                [Basis,SpaceGroup] = obj.read_crystal();
            end
        end
        
        function info = write(obj,name,mapdata,P,Basis,SpaceGroup)
            assert(nargin >= 3 && ~isempty(mapdata) && ~isempty(name))
            
            if nargin >= 4 && ~isempty(P)
                obj.write_grid(P);
            end
            
            if nargin == 6 && ~isempty(Basis) && ~isempty(SpaceGroup)
                obj.write_crystal(Basis,SpaceGroup)
            end
            
            info = obj.write_map(name,mapdata);
            
        end
        
    end

%     methods(Static)
%         function [P,varargout] = dhkl2map(hklTable,cols,ndiv,hklRange)
%             % hklTable has columns with integer indices h, k, l, dh, dk, dl, ...
%             % or with fractional indices h, k, l, .... 
%             % the program makes a determination based on whether dh is
%             % present or not.
%             %
%            
%             if nargin < 3 || isempty(hklRange)
%                 hklRange = [...
%                     min(hklTable.h),max(hklTable.h),...
%                     min(hklTable.k),max(hklTable.k),...
%                     min(hklTable.l),max(hklTable.l)];
%             end
%             
%             % prep grid
%             N = hklRange([2,4,6]) - hklRange([1,3,5]) + 1;
%             P = latt.PeriodicGrid(N.*ndiv,hklRange([1,3,5])-floor(ndiv/2)./ndiv,N);
%             
%         end
%     end
end