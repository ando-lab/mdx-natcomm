classdef MapTools < util.propertyValueConstructor
    %MAPTOOLS
    
    properties
        Grid (1,1) latt.PeriodicGrid = latt.PeriodicGrid([1,1,1],[0,0,0],[1,1,1])
        Basis (1,1) latt.Basis = latt.Basis(1,1,1,90,90,90)
        SpaceGroup (1,1) symm.SpaceGroup = symm.SpaceGroup(1)
        type {mustBeMember(type,{'density','structurefactor','patterson','intensity'})} = 'density'
        isPeriodic (1,1) logical = false
    end
    properties(Dependent)
        Symmetry
    end
    
    methods
        function obj = MapTools(varargin)
            %MAPTOOLS
            obj@util.propertyValueConstructor(varargin{:});
        end
        
        function S = get.Symmetry(obj)
            switch obj.type
                case 'density'
                    S = obj.SpaceGroup;
                case 'structurefactor'
                    S = obj.SpaceGroup;
                case 'patterson'
                    S = obj.SpaceGroup.LaueGroup;
                case 'intensity'
                    S = obj.SpaceGroup.LaueGroup;
            end
        end
        
        function [newobj,resizefun] = resize(obj,mode,varargin)
            
            switch lower(mode)
                case 'symexpand'
                    assert(ismember(obj.type,{'intensity'})); % other types have not been implemented / tested yet
                    Ops = obj.Symmetry.generalPositions;
                    [NewGrid,resizefun] = resize_for_symexpand(obj.Grid,Ops);
                case 'asu'
                    assert(ismember(obj.type,{'intensity'})); % other types have not been implemented / tested yet
                    [nmin,nmax] = mask_extents(obj.isASU);
                    [NewGrid,resizefun] = crop_grid(obj.Grid,nmin,nmax);
                case 'radius'
                    r = varargin{1};
                    if numel(varargin)==1
                        ori = [0,0,0];
                    else
                        ori = varargin{2};
                    end
                    B = obj.Basis.orthogonalizationMatrix();
                    [fmin,fmax] = bounding_box(B,ori,r);
                    [NewGrid,resizefun] = resize_grid_to_bounds(obj.Grid,fmin,fmax);
                case 'factor'
                    if numel(varargin)==1
                        if numel(varargin{1})==1
                            f = [1,1,1]*varargin{1};
                        else
                            assert(numel(varargin{1})==3)
                            f = varargin{1}(:)';
                        end
                    else
                        assert(numel(varargin)==3)
                        f = [varargin{1},varargin{2},varargin{3}];
                    end
                    [fmin,fmax] = grid_extents(obj.Grid);
                    [NewGrid,resizefun] = resize_grid_to_bounds(obj.Grid,fmin.*f,fmax.*f);
                case 'roi'
                    if numel(varargin)==1
                        farg = varargin{1};
                        assert(numel(farg)==6);
                        fmin = farg([1,3,5]);
                        fmax = farg([2,4,6]);
                    elseif numel(varargin)==2
                        fmin = varargin{1};
                        fmax = varargin{2};
                    else
                        error('incorrect number of arguments');
                    end
                    [NewGrid,resizefun] = resize_grid_to_bounds(obj.Grid,fmin,fmax);
                    % roi of the form [f1min,f1max,f2min,f2max,f3min,f3max]
                case 'grid'
                    assert(numel(varargin)==1 && isa(varargin{1},'latt.PeriodicGrid'));
                    [fmin,fmax] = grid_extents(varargin{1});
                    [NewGrid,resizefun] = resize_grid_to_bounds(obj.Grid,fmin,fmax);
                otherwise
                    error('resize mode not recognized');
            end
            
            newobj = obj;
            newobj.Grid = NewGrid;
            
            % resize to smax
            %
            % 1. symmetry expand the grid corners
            
        end
        
        function newobj = invert(obj)
            newobj = obj;
            newobj.Grid = obj.Grid.invert();
            newobj.Basis = obj.Basis.invert();
            switch obj.type
                case 'density'
                    newobj.type = 'structurefactor';
                case 'structurefactor'
                    newobj.type = 'density';
                case 'patterson'
                    newobj.type = 'intensity';
                case 'intensity'
                    newobj.type = 'patterson';
            end
        end
        
        function [newdata,newobj] = fourier_transform(obj,data,densityIsReal)
            
            if nargin < 3 || isempty(densityIsReal)
                densityIsReal = true;
            end
            
            LG = latt.LatticeGrid(obj.Grid,obj.Basis);
            
            newdata = data;
            newdata(isnan(newdata)) = 0;
            
            newobj = obj.invert();
            
            switch obj.type
                case 'density'
                    newdata = LG.ifft(newdata);
                case 'structurefactor'
                    newdata = LG.fft(newdata);
                    if densityIsReal, newdata = real(newdata); end
                case 'patterson'
                    newdata = LG.ifft(newdata);
                    if densityIsReal, newdata = real(newdata); end
                case 'intensity'
                    newdata = real(LG.fft(newdata));
            end
        end
        
        function [newdata,newobj] = fourier_interpolate(obj,data,ffactor)
            [d1,M1] = obj.fourier_transform(data);
            [M2,resizefun] = M1.resize('factor',ffactor);
            [newdata,newobj] = M2.fourier_transform(resizefun(d1));
        end
        
        function msk = spherical_mask(obj,radius)
            % spherical_mask
            %
            % if radius is of the form [rmin,rmax] then the mask is true
            % when rmin <= r <= rmax
            %
            [x,y,z] = obj.Grid.grid();
            [x,y,z] = obj.Basis.frac2lab(x,y,z);
            r2 = x.*x + y.*y + z.*z;
            if numel(radius)==1
                msk = r2  < radius^2;
            else
                assert(numel(radius)==2);
                msk = r2 < radius(2)^2 & r2 >= radius(1)^2;
            end
        end
        
        function asumask = isASU(obj,h,k,l)
            S = obj.Symmetry;
            assert(isa(S,'symm.LaueGroup')); % only implemented asus on this so far
            if nargin == 1
                [h,k,l] = obj.Grid.grid();
            end
            ndiv = obj.Grid.invert.P;
            h = round(ndiv(1)*h);
            k = round(ndiv(2)*k);
            l = round(ndiv(3)*l);
            asumask = S.testASU(h,k,l);
        end
        
        function h = export_mrc(obj,mrcfilename,M)
            
            % symmetry info
            S = obj.Symmetry;
            spaceGroupNumber = S.Info.number;
            symops = {S.generalPositions.xyzForm};

            % compute map extents
            n1 = [1,obj.Grid.N(1)+1];
            n2 = [1,obj.Grid.N(2)+1];
            n3 = [1,obj.Grid.N(3)+1];
            [f1,f2,f3] = obj.Grid.ind2frac(n1,n2,n3);
            a = obj.Basis.a*(f1(2)-f1(1));
            b = obj.Basis.b*(f2(2)-f2(1));
            c = obj.Basis.c*(f3(2)-f3(1));
            
            [o1,o2,o3] = obj.Grid.frac2ind(0,0,0,false);
            o = 1-[o1,o2,o3];
            
            h = io.map.initHeader();
            h.nc = obj.Grid.N(h.mapc);
            h.nr = obj.Grid.N(h.mapr);
            h.ns = obj.Grid.N(h.maps);
            h.nx = obj.Grid.N(1);
            h.ny = obj.Grid.N(2);
            h.nz = obj.Grid.N(3);
            h.x_length = a;
            h.y_length = b;
            h.z_length = c;
            h.alpha = obj.Basis.alpha;
            h.beta = obj.Basis.beta;
            h.gamma = obj.Basis.gamma;
            
            h.ncstart = o(h.mapc);
            h.nrstart = o(h.mapr);
            h.nsstart = o(h.maps);
            
            h.ispg = spaceGroupNumber;%0; % -1 <--- NOTE: space group is set to P1 by default... I've not figured out symmetry operators in map files yet
            h.nsymbt = numel(symops)*80;
            h.sym = char(join(pad(symops,80)',2));
            
            if nargin > 1          
                io.map.write(mrcfilename,h,M);
            end

        end
        
        function M = export(obj,h5out,mapname,varargin)
            ndsets = -1 + find(~cellfun(@isnumeric,varargin),1,'first');
            if isempty(ndsets), ndsets = numel(varargin); end
            kwopts = varargin((ndsets+1):end);
            
            % kw args alter opts, or get passed to M.new_dataset
            opts = struct('append',false,'newmap',false);
            dset_ops = {};
            for j=2:2:numel(kwopts)
                if isfield(opts,kwopts{j-1})
                    opts.(kwopts{j-1}) = kwopts{j};
                else
                    dset_ops = [dset_ops,kwopts((j-1):j)];
                end
            end
            
            M = io.h5.MapFile(h5out);
            M.verbose = true;
            
            if opts.append
                assert(M.file_exists);
                assert(M.group_exists('/crystal'));
                assert(M.group_exists('/maps'));
                warning('Appending. Beware: crystal and grid parameters are not checked for consistency.');
            elseif opts.newmap
                assert(M.file_exists);
                assert(M.group_exists('/crystal'));
                if ~M.group_exists('/maps')
                    M.addgroup('/maps');
                end
                warning('Creating new map in existing file. Beware: crystal is not checked for consistency.');
            
            else
                
                if M.file_exists, warning('file %s exists and is being overwritten',h5out); delete(h5out); end
                
                switch obj.type
                    case {'density','patterson'}
                        M.init(obj.Basis,obj.SpaceGroup);
                    case {'structurefactor','intensity'}
                        M.init(obj.Basis.invert,obj.SpaceGroup);
                end
                
                M.addgroup('/maps');
            end
            
            mpath = fullfile('/maps',mapname);
            
            if ~M.group_exists(mpath)
                M.new_grid(mpath,obj.Grid,'mdx_class',class(obj),'type',obj.type,'isPeriodic',uint8(obj.isPeriodic));
            end
            
            for j=1:ndsets
                varname = inputname(3 + j);
                if isempty(varname), varname = sprintf('data_%d',j); end
                M.new_dataset(fullfile(mpath,varname),varargin{j},dset_ops{:});
            end
            
        end
        
        %         function [Tasu,ASUGrid] = symreduce(obj,T,afun)
        %
        %             if nargin < 3 || isempty(afun)
        %                 afun = @mean;
        %             end
        %
        %             % 1) symexpand each index and return those that are within the
        %             % asu
        %
        %             f = table2array(T(:,1:3));
        %
        %             Ops = obj.Symmetry.generalPositions;
        %
        %             nOps = numel(Ops);
        %             fasu = f;
        %             for n=1:nOps
        %                 fs = (Ops(n)*f')';
        %                 isasu = obj.isASU(fs(:,1),fs(:,2),fs(:,3));
        %                 fasu(isasu,:) = fs(isasu,:);
        %             end
        %             [n1,n2,n3] = obj.Grid.frac2ind(fasu(:,1),fasu(:,2),fasu(:,3),false);
        %             [~,ia,ib] = unique([n1,n2,n3],'rows');
        %
        %             Tasu = array2table([fasu(ia,:), accumarray(ib,T.(4),[numel(ia),1],afun)],...
        %                 'VariableNames',T.Properties.VariableNames(1:4));
        %
        %             [fmin,fmax] = bounds(table2array(Tasu(:,1:3)));
        %             ASUGrid = obj.resize2bounds(obj.Grid,fmin,fmax);
        %
        %             %
        %             %             asumask = obj.isASU;
        %             %             MTasu = obj;
        %             %
        %             %             if resizegrid
        %             %                 [f1,f2,f3] = obj.Grid.ind2frac([1;obj.Grid.N(1)],[1;obj.Grid.N(2)],[1;obj.Grid.N(3)]);
        %             %                 Ops = obj.Symmetry.generalPositions;
        %             %                 [fmin,fmax] = bounds(obj.applysymops(Ops,[f1,f2,f3]));
        %             %                 MTfull = obj;
        %             %                 MTfull.Grid = MTfull.resize2bounds(MTfull.Grid,fmin,fmax);
        %             %                 [NewGrid,asumask] = MTfull.shrink2mask(MTfull.Grid,MTfull.isASU);
        %             %                 MTasu.Grid = NewGrid;
        %             %             end
        %             %
        %             %             T = MTasu.array2table(reshape(1:numel(asumask),MTasu.Grid.N),asumask);
        %             %             asuInd = obj.symexpand(T);
        %             %             isIncl = ~isnan(asuInd);
        %             %             asuInd = asuInd(isIncl);
        %             %
        %             %             [n1,n2,n3] = ind2sub(MTasu.Grid.N,asuInd);
        %             %
        %             %             A = accumarray([n1,n2,n3],A0(isIncl),MTasu.Grid.N,afun,NaN);
        %
        %         end
        
        %         function [Tfull,FullGrid] = symexpand(obj,T)
        %
        %             f = table2array(T(:,1:3));
        %
        %             Ops = obj.Symmetry.generalPositions;
        %             nOps = numel(Ops);
        %
        %             Tfull = cell(nOps,1);
        %             for n=1:numel(Ops)
        %                 fn = (Ops(n)*f')';
        %                 Tfull{n} = T;
        %                 Tfull{n}.(1) = fn(:,1);
        %                 Tfull{n}.(2) = fn(:,2);
        %                 Tfull{n}.(3) = fn(:,3);
        %             end
        %             Tfull = cat(1,Tfull{:});
        %
        %             [n1,n2,n3] = obj.Grid.frac2ind(Tfull.(1),Tfull.(2),Tfull.(3),false);
        %             [~,ia] = unique([n1,n2,n3],'rows');
        %             Tfull = Tfull(ia,:);
        %
        %             [fmin,fmax] = bounds(table2array(Tfull(:,1:3)));
        %             FullGrid = obj.resize2bounds(obj.Grid,fmin,fmax);
        %
        %         end
        
        function [gridIndex,isOnGrid] = frac2ind(obj,f1,f2,f3)
            [n1,n2,n3] = obj.Grid.frac2ind(f1,f2,f3,obj.isPeriodic);
            isOnGrid = n1 >= 1 & n2 >= 1 & n3 >= 1 & ...
                n1 <= obj.Grid.N(1) & n2 <= obj.Grid.N(2) & n3 <= obj.Grid.N(3);
            gridIndex = sub2ind(obj.Grid.N,n1(isOnGrid),n2(isOnGrid),n3(isOnGrid));
        end
        
        function varargout = table2array(obj,T,mode,mode2)
            
            if nargin < 3 || isempty(mode)
                mode = 'direct'; % or 'symexpand'
            end
            
            if nargin < 4 || isempty(mode2)
                mode2 = 'replace'; % or 'mean' or 'sum'
            end
            
            fprintf(1,'mode = %s, mode2 = %s\n',mode,mode2);
            
            ncols = size(T,2) - 3; % number of columns to convert
            assert(ncols > 0);
            varargout = cell(1,ncols);
            
            % initialize outputs
            for n=1:ncols
                varargout{n} = zeros(obj.Grid.N);
            end
            
            switch lower(mode)
                case 'direct'
                    [ind,isIncl] = obj.frac2ind(T.(1),T.(2),T.(3));
                    for n=1:ncols
                        varargout{n}(ind) = T.(3 + n)(isIncl);
                    end
                case 'symexpand'
                    
                    switch obj.type
                        case {'density','patterson'}
                            symfun = @(op,f) (op*f')';
                        case 'intensity'
                            symfun = @(op,f) (f*op);
                        case 'structurefactor'
                            symfun = @(op,f) ([f,zeros(size(f,1),1)]*op);
                            % returns a 4-vector
                    end
                    
                    Ops = obj.Symmetry.generalPositions;
                    
                    x0 = table2array(T(:,1:3));
                    
                    N = zeros(obj.Grid.N);
                    
                    for j=1:numel(Ops)
                        xj = symfun(Ops(j),x0);
                        [ind,isIncl] = obj.frac2ind(xj(:,1),xj(:,2),xj(:,3));
                        if size(xj,1)==4
                            phaseFactor = exp(2i*pi*xj(isIncl,4));
                        else
                            phaseFactor = 1;
                        end
                        
                        N(ind) = N(ind) + 1;
                        
                        for n=1:ncols
                            
                            switch mode2
                                case 'replace'
                                    varargout{n}(ind) = phaseFactor.*T.(3 + n)(isIncl);
                                case 'mean'
                                    varargout{n}(ind) = varargout{n}(ind) + phaseFactor.*T.(3 + n)(isIncl);
                                otherwise
                                    error('mode2 not recognized');
                            end
                        end
                    end
                    
                    % post-processing
                    switch mode2
                        case 'mean'
                            for n=1:ncols
                                varargout{n} = varargout{n}./N;
                            end
                        case 'replace'
                            for n=1:ncols
                                varargout{n}(~logical(N)) = NaN;
                            end
                    end
                    
                otherwise
                    error('mode not recognized');
            end
            
        end
        
        function T = array2table(obj,A,mask)
            if isempty(A)
                A = zeros(obj.Grid.N);
            end
            if nargin < 3 || isempty(mask)
                mask = ~isnan(A);
            end
            [nmin,nmax] = mask_extents(mask);
            [G0,cropfun] = crop_grid(obj.Grid,nmin,nmax);
            m = cropfun(mask);
            A = cropfun(A);
            [f1,f2,f3] = G0.grid(); % <-- slow?
            
            switch obj.type
                case 'density'
                    colnames = {'x','y','z','rho'};
                case 'structurefactor'
                    colnames = {'h','k','l','F'};
                case 'patterson'
                    colnames = {'x','y','z','P'};
                case 'intensity'
                    colnames = {'h','k','l','I'};
            end
            
            T = table(f1(m),f2(m),f3(m),A(m),'VariableNames',colnames);
        end
        
        function varargout = read_data(obj,h5in,mapname,varargin)
            
           [MT,M,mpath] = obj.fromfile(h5in,mapname);
           
           try
               obj.check_compatibility(MT);
           catch EM
               fprintf(1,'error: grid in file is incompatible with MapTools object\n');
               rethrow(EM);
           end
           % need to check that maps are equivalent somehow
           
           % need to figure out where to start and stop reading!
           
           % first, find the bounds to read
           
           [nmin,nmax,resizefun] = resize_for_target_grid(MT.Grid,obj.Grid);
           
           varargout = cell(size(varargin));
           
           for j=1:numel(varargin)
                dset0 = M.read(fullfile(mpath,varargin{j}),nmin,1 + nmax - nmin);
                varargout{j} = resizefun(dset0);
           end
        end
        
        function check_compatibility(obj,MT)
            % check whether the two map pro.script.MapTools are compatible
            % with each other
            
            basistol = 1E-3;

            if ~isa(MT,class(obj))
                warning('class mismatch');
            end
            if ~strcmp(obj.type,MT.type)
                warning('type mismatch');
            end
            
            %assert(obj.isPeriodic == MT.isPeriodic);
            %assert(obj.SpaceGroup.number == MT.SpaceGroup.number);
            assert(abs(obj.Basis.a - MT.Basis.a)/obj.Basis.a < basistol);
            assert(abs(obj.Basis.b - MT.Basis.b)/obj.Basis.b < basistol);
            assert(abs(obj.Basis.c - MT.Basis.c)/obj.Basis.c < basistol);
            assert(abs(obj.Basis.alpha - MT.Basis.alpha)/obj.Basis.alpha < basistol);
            assert(abs(obj.Basis.beta - MT.Basis.beta)/obj.Basis.beta < basistol);
            assert(abs(obj.Basis.gamma - MT.Basis.gamma)/obj.Basis.gamma < basistol);
            
            % same grid delta
            assert(all(round(obj.Grid.N.*MT.Grid.P) == round(obj.Grid.P.*MT.Grid.N)));
            
        end
        
    end
    methods(Static)
        
        function [MT,M,mpath] = fromfile(h5in,mapname)
            M = io.h5.MapFile(h5in);
            
            if mapname(1)=='/'
                mpath = mapname;
            else
                mpath = fullfile('/maps',mapname);
            end
            
            [Basis,SpaceGroup] = M.read_crystal();
            
            Grid = M.read_grid(mpath);
            
            MT = proc.script.MapTools(...
                'Basis',Basis,...
                'SpaceGroup',SpaceGroup,...
                'Grid',Grid);
            
            a = M.readatt(mpath);
            if isfield(a,'mdx_class') && a.mdx_class == "proc.script.MapTools"
                MT.type = a.type;
                MT.isPeriodic = a.isPeriodic;
                switch lower(a.type)
                    case {'intensity','structurefactor'}
                        MT.Basis = MT.Basis.invert;
                end
            else
                warning('mdx_class is missing or did not match. assuming map represents electron density');
            end
            
        end
        
        function [MT,varargout] = import(h5in,mapname,varargin)
            
            [MT,M,mpath] = proc.script.MapTools.fromfile(h5in,mapname);
            
            varargout = cell(1,numel(varargin));
            for j=1:numel(varargin)
                varargout{j} = M.read(fullfile(mpath,varargin{j}));
            end
        end

    end
end



function fs = applysymops(Ops,f)
% assumes f is an array of column vectors [f1,f2,f3]
assert(size(f,2)==3);
nOps = numel(Ops);
fs = cell(nOps,1);
for n=1:nOps
    fs{n} = (Ops(n)*f')';
end
fs = cat(1,fs{:});
end

function [fmin,fmax] = grid_extents(G)

[f1,f2,f3] = G.ind2frac([1,G.N(1)],[1,G.N(2)],[1,G.N(3)]);
[fmin,fmax] = bounds([f1(:),f2(:),f3(:)]);

end

function [nmin,nmax,resizefun] = resize_for_target_grid(G_file,G_target)

% first, find the bounds to read
[fmin_target,fmax_target] = grid_extents(G_target);
[fmin_file,fmax_file] = grid_extents(G_file);

% get the bounds with respect to the file grid
nmin_target = f2n(G_file,fmin_target);
nmax_target = f2n(G_file,fmax_target);

nmin_file = f2n(G_file,fmin_file);
nmax_file = f2n(G_file,fmax_file);

if any(nmax_target > nmax_file | nmin_target < nmin_file)
    warning('Target grid has larger extent that data file. Will pad with NaN');
end

% now, figure out what portion of the file to read
nmin = max(nmin_target,nmin_file);
nmax = min(nmax_target,nmax_file);

G0 = resize_grid(G_file,nmin,nmax);

[~,resizefun] = resize_grid_to_bounds(G0,fmin_target,fmax_target);

    function n = f2n(G,f)
        % for convenience
        [n_1,n_2,n_3] = G.frac2ind(f(:,1),f(:,2),f(:,3),false);
        n = [n_1,n_2,n_3];
    end

end

function [NewGrid,resizefun] = resize_for_symexpand(G,Ops)

[n1,n2,n3] = ndgrid([1,G.N(1)],[1,G.N(2)],[1,G.N(3)]);
[f1,f2,f3] = G.ind2frac(n1,n2,n3);
f = [f1(:),f2(:),f3(:)];
fs = applysymops(Ops,f);
[fmin,fmax] = bounds(fs);

[NewGrid,resizefun] = resize_grid_to_bounds(G,fmin,fmax);

end

function [NewGrid,resizefun] = resize_grid_to_bounds(G,fmin,fmax)

[n1min,n2min,n3min] = G.frac2ind(fmin(1),fmin(2),fmin(3),false);
[n1max,n2max,n3max] = G.frac2ind(fmax(1),fmax(2),fmax(3),false);

[NewGrid,resizefun] = resize_grid(G,[n1min,n2min,n3min],[n1max,n2max,n3max]);

end

function [NewGrid,resizefun] = resize_grid(G,nmin,nmax)

N = 1 + nmax - nmin;
P = G.delta.*N;
[o1,o2,o3] = G.ind2frac(nmin(1),nmin(2),nmin(3));

NewGrid = latt.PeriodicGrid(N,[o1,o2,o3],P);

ncmin = max(nmin,[1,1,1]);
ncmax = min(nmax,G.N);

npadpre = ncmin - nmin;
npadpost = nmax - ncmax;
padval = NaN;

resizefun = @(m) ...
    padarray(...
    padarray(...
    m(ncmin(1):ncmax(1),ncmin(2):ncmax(2),ncmin(3):ncmax(3)),...
    npadpre,padval,'pre'),...
    npadpost,padval,'post');

end

function [NewGrid,cropfun] = crop_grid(G,nmin,nmax)

N = 1 + nmax - nmin;
P = G.delta.*N;
[o1,o2,o3] = G.ind2frac(nmin(1),nmin(2),nmin(3));

NewGrid = latt.PeriodicGrid(N,[o1,o2,o3],P);

cropfun = @(m) m(nmin(1):nmax(1),nmin(2):nmax(2),nmin(3):nmax(3));

end

function [nmin,nmax] = mask_extents(mask)

mask3 = squeeze(any(mask,[1,2]));
mask2 = squeeze(any(mask,[1,3]));
mask1 = squeeze(any(mask,[2,3]));
n1min = find(mask1,1,'first');
n1max = find(mask1,1,'last');
n2min = find(mask2,1,'first');
n2max = find(mask2,1,'last');
n3min = find(mask3,1,'first');
n3max = find(mask3,1,'last');

nmin = [n1min,n2min,n3min];
nmax = [n1max,n2max,n3max];

end



% adapted from proc.script.CoordinateTools
function [fmin,fmax] = bounding_box(B,ori,rmax)

ori = ori(:);

v3 = cross(B(:,1),B(:,2));
v3 = v3/norm(v3);

v1 = cross(B(:,2),B(:,3));
v1 = v1/norm(v1);

v2 = cross(B(:,3),B(:,1));
v2 = v2/norm(v2);

p1 = ori + [-1,1].*v1*rmax;
p2 = ori + [-1,1].*v2*rmax;
p3 = ori + [-1,1].*v3*rmax;

f1 = inv(B)*p1;
f2 = inv(B)*p2;
f3 = inv(B)*p3;

f1 = sort(f1(1,:),'ascend');
f2 = sort(f2(2,:),'ascend');
f3 = sort(f3(3,:),'ascend');

fmin = [f1(1),f2(1),f3(1)];
fmax = [f1(2),f2(2),f3(2)];

end