classdef LatticeDynamicsTools < util.propertyValueConstructor
    %LatticeDynamicsTools - fit lattice dynamics models to ADPs, diffuse
    %scattering, and perform simulations
    
    properties
        Cell
        M
        interp_osr = 4
        interp_dgrid = 0.5
        interp_Badd = 5
        interp_rmax = 4
        interp_mode = 'cubic'
        addLinearBackground = true;
        supercell = [1,1,1]
    end
    properties(Dependent)
        Basis
    end
    
    methods
        function obj = LatticeDynamicsTools(varargin)
            %LatticeDynamicsTools
            obj@util.propertyValueConstructor(varargin{:});
        end
        
        function B = get.Basis(obj)
            B = obj.Cell.Basis;
            if ~isa(B,'latt.OrientedBasis')
                B = B.orient;
            end
        end
        
        function [sAxes,ori,F,Fx,Fy,Fz] = map2sfgrid(obj,rho,MapGrid,MapBasis,tlsori)
            % EXPERIMENTAL FEATURE!
            assert(MapBasis.alpha==90 & MapBasis.beta==90 & MapBasis.gamma==90); % required for now
                        
            MT = proc.script.MapTools('Grid',MapGrid,'Basis',MapBasis,'type','density','isPeriodic',false);
            [x,y,z] = MT.Grid.grid();
            [x,y,z] = MT.Basis.frac2lab(x,y,z);
            x = x-tlsori(1);
            y = y-tlsori(2);
            z = z-tlsori(3);
            
            ori0 = MT.Grid.ori;
            MT.Grid = MT.Grid.invert.invert; % re-center the grid
            oshift = ori0 - MT.Grid.ori;
            [ox,oy,oz] = MT.Basis.frac2lab(oshift(1),oshift(2),oshift(3));
            ori = [ox,oy,oz];
            
            [MT,resizefun] = MT.resize('factor',obj.interp_osr);
            rho = resizefun(rho);
            x = resizefun(x);
            y = resizefun(y);
            z = resizefun(z);
            
            [F,MTr] = MT.fourier_transform(rho);
            Fx = MT.fourier_transform(rho.*x);
            Fy = MT.fourier_transform(rho.*y);
            Fz = MT.fourier_transform(rho.*z);
            
            [MTr,cropfun] = MTr.resize('radius',0.5/obj.interp_dgrid);
            F = cropfun(F);
            Fx = cropfun(Fx);
            Fy = cropfun(Fy);
            Fz = cropfun(Fz);
            
            [hv,kv,lv] = MTr.Grid.ind2frac(1:MTr.Grid.N(1),1:MTr.Grid.N(2),1:MTr.Grid.N(3));
            sxv = MTr.Basis.a*hv(:);
            syv = MTr.Basis.b*kv(:);
            szv = MTr.Basis.c*lv(:);
            
            sAxes = {sxv,syv,szv};
        end
        
        function [sAxes,ori,F,Fx,Fy,Fz] = atoms2sfgrid(obj,AtomFF)
            
            a = range(AtomFF.x)*obj.interp_osr;
            b = range(AtomFF.y)*obj.interp_osr;
            c = range(AtomFF.z)*obj.interp_osr;
            BigGrid = proc.script.GridDesigner(...
                'Basis',latt.Basis(a,b,c,90,90,90).orient,...
                'dgrid',obj.interp_dgrid,...
                'maxPrimeFactor',5).run();
            ori = mean(table2array(AtomFF(:,{'x','y','z'})),1);
            SF = proc.script.StructureFactor(...
                'LatticeGrid',BigGrid,...
                'Badd',obj.interp_Badd,...
                'rmax',obj.interp_rmax);
            x0 = AtomFF.x;
            y0 = AtomFF.y;
            z0 = AtomFF.z;
            
            AtomFF.x = x0 - ori(1);
            AtomFF.y = y0 - ori(2);
            AtomFF.z = z0 - ori(3);
            
            Ax = AtomFF; Ay = AtomFF; Az = AtomFF;
            Ax.mdxFormFactor = AtomFF.mdxFormFactor.rescale(x0);
            Ay.mdxFormFactor = AtomFF.mdxFormFactor.rescale(y0);
            Az.mdxFormFactor = AtomFF.mdxFormFactor.rescale(z0);
            
            F = SF.run(AtomFF);
            Fx = SF.run(Ax);
            Fy = SF.run(Ay);
            Fz = SF.run(Az);
            
            RGrid = BigGrid.invert;
            
            [hv,kv,lv] = RGrid.PeriodicGrid.ind2frac(1:BigGrid.N(1),1:BigGrid.N(2),1:BigGrid.N(3));
            sxv = RGrid.Basis.a*hv(:);
            syv = RGrid.Basis.b*kv(:);
            szv = RGrid.Basis.c*lv(:);
            
            sAxes = {sxv,syv,szv};
            
        end
        
        function sffun = calc1PSFInterpFromMap(obj,rho,MapGrid,MapBasis)
            % EXPERIMENTAL
            
            % for now, this only supports one rigid group. In future, 
            % will implement segmented maps
            nGroups = numel(obj.Cell.AsymmetricUnit);
            assert(nGroups==1);
            
            ori = cell(nGroups,1);
            GI = cell(nGroups,1);
            
            g = 1;
            
            tlsori = obj.Cell.AsymmetricUnit(g).ori;
            if isempty(tlsori)
                tlsori = [0,0,0];
            end
            
            fprintf(1,'computing interpolants for group %d of %d\n',g,nGroups);
                
            [sAxes,ori{g},F,Fx,Fy,Fz] = obj.map2sfgrid(rho,MapGrid,MapBasis,tlsori);
            
            GI{g} = griddedInterpolant(sAxes,...
                    cat(4,real(F),imag(F),real(Fx),imag(Fx),real(Fy),imag(Fy),real(Fz),imag(Fz)),...
                    obj.interp_mode,'none');
                
            sffun = obj.interpolant2sffun(GI,ori);
            
            fprintf(1,'done\n');
            
        end
        
        function sffun = calc1PSFInterp(obj,AtomFF,groupid)
        
            nGroups = numel(obj.Cell.AsymmetricUnit);
            assert(all(groupid <= nGroups & groupid >= 1));
            
            ori = cell(nGroups,1);
            GI = cell(nGroups,1);
            
            for g=1:nGroups
                
                A = AtomFF(groupid==g,:);
                tlsori = obj.Cell.AsymmetricUnit(g).ori;
                if isempty(tlsori)
                    tlsori = [0,0,0];
                end
                
                    A.x = A.x - tlsori(1);
                    A.y = A.y - tlsori(2);
                    A.z = A.z - tlsori(3);
                    
                fprintf(1,'computing interpolants for group %d of %d\n',g,nGroups);
                
                [sAxes,ori{g},F,Fx,Fy,Fz] = obj.atoms2sfgrid(A);
                
                ori{g} = ori{g} + tlsori(:)';
                
                GI{g} = griddedInterpolant(sAxes,...
                    cat(4,real(F),imag(F),real(Fx),imag(Fx),real(Fy),imag(Fy),real(Fz),imag(Fz)),...
                    obj.interp_mode,'none');
            end
            
            sffun = obj.interpolant2sffun(GI,ori);
            
            fprintf(1,'done\n');
            
        end
        
        function sffun = interpolant2sffun(obj,GI,ori)
            
            RBasis = obj.Basis.invert;
            Ops = obj.Cell.UnitCellOperators;
            
            sffun = @(h,k,l) mysffun(h,k,l,Ops,RBasis,GI,ori);
            
            fprintf(1,'done\n');
            
            function G = mysffun(hh,kk,ll,Ops,RBasis,GI,ori)
                nops = numel(Ops);
                ngroups = numel(ori);
                G = cell(nops,ngroups);
                for grp=1:ngroups
                    for o=1:nops
                        hkln = [hh,kk,ll,zeros(size(hh))]*Ops(o);
                        [sx,sy,sz] = RBasis.frac2lab(hkln(:,1),hkln(:,2),hkln(:,3));
                        pf = exp(2i*pi*(hkln(:,4) + sx*ori{grp}(1) + sy*ori{grp}(2) + sz*ori{grp}(3)));
                        Fa = GI{grp}(sx,sy,sz);
                        Fn = pf.*(Fa(:,:,:,1) + 1i*Fa(:,:,:,2));
                        Fxn = pf.*(Fa(:,:,:,3) + 1i*Fa(:,:,:,4));
                        Fyn = pf.*(Fa(:,:,:,5) + 1i*Fa(:,:,:,6));
                        Fzn = pf.*(Fa(:,:,:,7) + 1i*Fa(:,:,:,8));
                        G{grp,o} = 2*pi*[...
                            sx.*Fn,...
                            sy.*Fn,...
                            sz.*Fn,...
                            Fyn.*sz - Fzn.*sy,...
                            Fzn.*sx - Fxn.*sz,...
                            Fxn.*sy - Fyn.*sx];
                    end
                end
                G = cat(2,G{:});
            end
            
        end
        
        function [Gk,kspace_group] = precompute1PSFs(obj,sffun,h,k,l)

            [kspace_group,G_bz] = obj.kspace_groupings(h,k,l);
            
            Gk = cell(G_bz.N);
            
            nBZ = numel(Gk);
            
            for j=1:nBZ
                fprintf(1,'computing 1-phonon structure factor for k-vector %d of %d\n',j,nBZ);
                kg = kspace_group{j};
                if isempty(kg)
                    continue;
                end
                Gk{j} = sffun(h(kg),k(kg),l(kg));
            end
        end

        
        function [ind,G_bz] = kspace_groupings(obj,h,k,l)
            
            fprintf(1,'grouping h,k,l by k-vector\n');
            
            G_bz = latt.PeriodicGrid(obj.supercell,[0,0,0],[1,1,1]).invert.invert;
            
            [n1,n2,n3] = G_bz.frac2ind(h(:),k(:),l(:));
            
            ind = accumarray([n1,n2,n3],1:numel(h),G_bz.N,@(v) {v});
        end

        
        function I = calc1PIntensityUnitcell(obj,V,sffun,h,k,l)
            
            [kspace_group,G_bz] = obj.kspace_groupings(h,k,l);
            
            LD = nm.LatticeDynamics('V',V,'supercell',obj.supercell,'M',obj.M);
            C = LD.cov_unitcell;
            
            Ik = cell(G_bz.N);
            nBZ = numel(Ik);
            
            for j=1:nBZ
                kg = kspace_group{j};
                if isempty(kg)
                    continue;
                end
                fprintf(1,'computing 1-phonon structure factor for k-vector %d of %d\n',j,nBZ);
            
                Gk = sffun(h(kg),k(kg),l(kg));
                Ik{j} = real(dot(Gk,Gk*C,2));
            end
            
            I = zeros(size(h));
            for j=1:numel(Ik)
                I(kspace_group{j}) = Ik{j};
            end
        end
        
        function I = calc1PIntensity(obj,V,sffun,h,k,l)

            [kspace_group,G_bz] = obj.kspace_groupings(h,k,l);
            
            LD = nm.LatticeDynamics('V',V,'supercell',obj.supercell,'M',obj.M);
            Kinv = LD.Kinv;
            %C = real(mean(Kinv,[1,2,3]));
            %C = shiftdim(C,3);
            Kinv = shiftdim(Kinv,3);
            
            Ik = cell(G_bz.N);
            %I0k = cell(G_bz.N);
            nBZ = numel(Ik);
            
            for j=1:nBZ
                kg = kspace_group{j};
                if isempty(kg)
                    continue;
                end
                fprintf(1,'computing 1-phonon structure factor for k-vector %d of %d\n',j,nBZ);
            
                Gk = sffun(h(kg),k(kg),l(kg));
                Ik{j} = real(dot(Gk,Gk*conj(Kinv(:,:,j)),2));
                %I0k{j} = real(dot(Gk,Gk*C,2));
                %Ik{j} = dot(Gk,Gk*conj(Kinv(:,:,j)),2);
            end
            
            I = zeros(size(h));
            %I0 = zeros(size(h));
            for j=1:numel(Ik)
                I(kspace_group{j}) = Ik{j};
                %I0(kspace_group{j}) = I0k{j};
            end
            
        end
        
        
        function I = calc1PIntensityFrom1PSF(obj,V,Gk,ind,sz)
            
            LD = nm.LatticeDynamics('V',V,'supercell',obj.supercell,'M',obj.M);
            Kinv = LD.Kinv;
            Kinv = shiftdim(Kinv,3);
            
            Ik = cell(size(Gk));
            nBZ = numel(Gk);
            
            for j=1:nBZ
                if isempty(Gk{j})
                    continue;
                end
                Ik{j} = real(dot(Gk{j},Gk{j}*conj(Kinv(:,:,j)),2));
            end
            
            I = zeros(sz);
            for j=1:numel(Ik)
                I(ind{j}) = Ik{j};
            end
            
        end
        

        function [pfit,fitinfo,history] = fitHessianToHalos(obj,I,sigma,Gk,ind,Vfun,p0,pmin,pmax,varargin)
            
            assert(all(size(I,[2,3,4])==obj.supercell));
            
            sz = size(I);
            
            Icalc = @(p) obj.calc1PIntensityFrom1PSF(Vfun(p),Gk,ind,sz);
            
            optfun = @(p) haloFitFunction(I,sigma,Icalc(p),obj.addLinearBackground);
            
            [pfit,fitinfo,history] = run_lsqnonlin(optfun,p0,pmin,pmax,varargin{:});
        end
        
        function [P,massweights] = calcProjectionOperator(obj,toCOM)
            if nargin < 2 || isempty(toCOM)
                toCOM = true;
            end
            
            UC = obj.Cell;
            
            if toCOM
                coms = arrayfun(@(g) [g.com; sum(g.massvec)],UC.AsymmetricUnit,'Uni',0);
                UC.AsymmetricUnit = cellfun(@(v) nm.Group(v(1),v(2),v(3),v(4)),coms);
            end
            
            P = UC.tl2uxyz;
            
            massweights = repmat([UC.AsymmetricUnit.massvec],1,numel(UC.SpaceGroup.generalPositions));

        end

        function V = calcLatticeCovfromHessian(obj,Hessian,ProjectionOperator,massweights)
            ngroups = size(ProjectionOperator,1)/3;
            w = ones([1,1,ngroups]);
            
            if nargin < 4 || isempty(massweights)
                w(:) = w(:).*massweights(:);
            end
            
            LD = nm.LatticeDynamics('V',Hessian,'supercell',obj.supercell,'M',obj.M);
            C = LD.cov_supercell();
            C = shiftdim(mat2cell(C,size(C,1),size(C,2),ones(size(C,3),1),ones(size(C,4),1),ones(size(C,5),1)),2);
            V = cellfun(@(c) sum(w.^2.*get_diag_projection(c,ProjectionOperator),3)/sum(w.^2),C,'Uni',0);
        end
        
        function [pfit,fitinfo,history] = fitHessianToLatticeCov(obj,Vobs,weights,avgsub,Vfun,p0,pmin,pmax,varargin)
            
            % by default, use the center of mass of the ASU groups to
            % compute this guy?
            [P,mw] = obj.calcProjectionOperator(true);
            
            Vcalc = @(p) obj.calcLatticeCovfromHessian(Vfun(p),P,mw);
            
            optfun = @(p) covFitFun(Vcalc(p),Vobs,weights,avgsub);
            
            [pfit,fitinfo,history] = run_lsqnonlin(optfun,p0,pmin,pmax,varargin{:});
        end
        
        function [pfit,fitinfo,history] = fitHessianToCov(obj,Cobs,weights,Vfun,p0,pmin,pmax,varargin)

            Vcalc = @(p) nm.LatticeDynamics('V',Vfun(p),'supercell',obj.supercell,'M',obj.M).cov_supercell;

            optfun = @(p) (Vcalc(p) - Cobs).*weights;
            
            [pfit,fitinfo,history] = run_lsqnonlin(optfun,p0,pmin,pmax,varargin{:});
        end
        
        function [pfit,fitinfo,history] = fitHessianToADPs(obj,Uobs,Vfun,p0,pmin,pmax,varargin)
            % Uobs is a cell array of 3x3 matrices (U)
            Uobs = cat(3,Uobs{:}); % <- turn into a 3x3xN
            P = obj.Cell.AsymmetricUnit.tl2uxyz;
                        
            optfun = @(p) Uobs - calcUfromHessian(Vfun(p),obj.M,P,obj.supercell);
            
            [pfit,fitinfo,history] = run_lsqnonlin(optfun,p0,pmin,pmax,varargin{:});
            
            %H = Vfun(pfit);
            %Ucalc = calcUfromHessian(Vfun(pfit),P,obj.supercell);
            %Ucalc = mat2cell(Ucalc,3,3,ones(size(Ucalc,3),1));
        end
        
        function [pfit,fitinfo,history] = fitHessianToADPs_weighted(obj,Uobs,Vfun,p0,pmin,pmax,varargin)
            % Uobs is a cell array of 3x3 matrices (U)
            w = 1./cellfun(@trace,Uobs(:));
            w = shiftdim(w,-2);
            
            Uobs = cat(3,Uobs{:}); % <- turn into a 3x3xN
            P = obj.Cell.AsymmetricUnit.tl2uxyz;
                        
            optfun = @(p) w.*(Uobs - calcUfromHessian(Vfun(p),obj.M,P,obj.supercell));
            
            [pfit,fitinfo,history] = run_lsqnonlin(optfun,p0,pmin,pmax,varargin{:});
            
            %H = Vfun(pfit);
            %Ucalc = calcUfromHessian(Vfun(pfit),P,obj.supercell);
            %Ucalc = mat2cell(Ucalc,3,3,ones(size(Ucalc,3),1));
        end

        
    end
    
    methods(Static)
        
        function [hklRef,hklTable] = findStrongBraggPeaks(mtzFileName,numRef,resMin,resMax)
            
            [hklTable,Basis] = proc.script.ImportMTZ().run(mtzFileName);
            
            % filter by reflection range
            [sx,sy,sz] = Basis.invert.frac2lab(hklTable.h,hklTable.k,hklTable.l);
            s = sqrt(sx.^2 + sy.^2 + sz.^2);
            isIncl = 1./s >= resMin & 1./s <= resMax & ~isnan(hklTable.Fobs);
            hklTable = hklTable(isIncl,:);
            
            % sort by Fobs and choose the most intense entries
            [~,ixorder] = sort(hklTable.Fobs,'descend');
            hklTable = hklTable(ixorder(1:numRef),:);
            
            hklRef = table2array(hklTable(:,{'h','k','l'}));
            
        end
        
        function [hklGrid,I,sigma] = loadBrillouinZones(h5In,dset,hklRef)
            
            M = io.h5.MapFile(h5In);
            G = M.read_grid(dset);
            ndiv = 1./G.delta;
            G_bz = latt.PeriodicGrid(ndiv,[0,0,0],[1,1,1]).invert.invert;
            
            hklGrid = repmat(G_bz,size(hklRef,1),1);
            for j=1:size(hklRef,1)
                hklGrid(j).ori = hklGrid(j).ori + hklRef(j,:);
            end
            
            I = NaN*ones([size(hklRef,1),ndiv]);
            sigma = NaN*ones([size(hklRef,1),ndiv]);
            
            for j=1:size(hklRef,1)
                ori = hklGrid(j).ori;
                [n1,n2,n3] = G.frac2ind(ori(1),ori(2),ori(3));
                try
                    Ij = M.read(fullfile(dset,'I'),[n1,n2,n3],hklGrid(j).N);
                    sigmaj = M.read(fullfile(dset,'sigma'),[n1,n2,n3],hklGrid(j).N);
                catch EM
                    warning('error thrown by h5read for h,k,l=%d,%d,%d:\n message = ''%''',...
                        hklRef(j,:),EM.message);
                    continue;
                end
                I(j,:,:,:) = Ij;
                sigma(j,:,:,:) = sigmaj;
            end
            
            sigma(isnan(I) | isnan(sigma)) = Inf;
            I(isinf(sigma)) = 0;
            
        end
        
        
        
    end
end



function [xfit,fitinfo,history] = run_lsqnonlin(objfun,x0,xmin,xmax,varargin)

% for example:
%varargin = {'Algorithm','trust-region-reflective',...
%    'MaxFunctionEvaluations',1000,'UseParallel',false};

history = struct('iteration',{},'x',{},'resnorm',{});

opts = optimoptions(@lsqnonlin,varargin{:},'OutputFcn',@outfun);

[xfit,~,~,~,fitinfo] = lsqnonlin(objfun,x0,xmin,xmax,opts);

    function stop = outfun(x,optimValues,state)
        stop = false;
        
        switch state
            case 'init'
                % do nothing
            case 'iter'
                % Concatenate current point and objective function
                history = [history;...
                    struct('iteration',optimValues.iteration,...
                    'x',x,'resnorm',optimValues.resnorm)];
            case 'done'
                % do nothing
            otherwise
        end
    end

end

function U = calcUfromHessian(H,M,P,supercell)
% returns a 3x3xN matrix
C = nm.LatticeDynamics('V',H,'supercell',supercell,'M',M).cov_unitcell;
nc = size(P,2);

U = get_diag_projection(C(1:nc,1:nc),P);

%P = mat2cell(P,3*ones(size(P,1)/3,1),nc);
%U = cellfun(@(p) p*C(1:nc,1:nc)*p',P,'Uni',0);
end


function resid = haloFitFunction(I,sigma,Icalc,addLinearBackground)

if nargin < 4 || isempty(addLinearBackground)
    addLinearBackground = true;
end

supercell = size(I,[2,3,4]);

% fit background (constant + linear)
G_bz = latt.PeriodicGrid(supercell,[0,0,0],[1,1,1]).invert.invert;
[dh,dk,dl] = G_bz.grid;

if addLinearBackground
    A = ones(prod(supercell),4);
    A(:,2) = dh(:);
    A(:,3) = dk(:);
    A(:,4) = dl(:);
else
    A = ones(prod(supercell),1);
end

B = zeros(size(I));

for jj=1:size(Icalc,1)
    x = (A.*repmat(1./sigma(jj,:)',1,size(A,2)))\((I(jj,:) - Icalc(jj,:))./sigma(jj,:))';
    B(jj,:) = A*x;
end

% residual
resid = (I - Icalc - B)./sigma;

end




function Md = get_diag_projection(C,Pcom)

npts = size(Pcom,1)/3;

Md = zeros(3,3,npts);

for n=1:npts
   n1start = (n-1)*3 + 1;
   n1stop = n1start + 3 - 1;
   Pn = Pcom(n1start:n1stop,:);
   Md(:,:,n) = Pn*C*Pn';
end

end


function resid = covFitFun(Vcalc,Vobs,weights,avgsub)

if nargin < 4 || isempty(avgsub)
    avgsub = true;
end

idx = find(triu(ones(3,3)));

Vdiff = cellfun(@(vcalc,vobs) vcalc(idx)-vobs(idx),Vcalc,Vobs,'Uni',0);
resid = cat(2,Vdiff{:});

if avgsub
    resid = resid - mean(resid,2); % remove offset (i.e. global aniso scaling)
end
resid = resid.*(weights(:)');

end
