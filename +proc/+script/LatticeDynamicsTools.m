classdef LatticeDynamicsTools < util.propertyValueConstructor
    %LatticeDynamicsTools - fit lattice dynamics models to ADPs, diffuse
    %scattering, and perform simulations
    
    properties
        Cell
        interp_osr = 4
        interp_dgrid = 0.5
        interp_Badd = 5
        interp_rmax = 4
        interp_mode = 'cubic'
        supercell = [1,1,1]
    end
    
    methods
        function obj = LatticeDynamicsTools(varargin)
            %LatticeDynamicsSimulation
            obj@util.propertyValueConstructor(varargin{:});
        end
        
        function [sAxes,ori,F,Fx,Fy,Fz] = calcStructureFactorInterpolants(obj,AtomFF)
            
            Basis = obj.Cell.Basis;
            if ~isa(Basis,'latt.OrientedBasis')
                Basis = Basis.orient;
            end
            
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
            
            %phx = exp(2i*pi*sxv*ori(1));
            %phy = exp(2i*pi*syv*ori(2));
            %phz = exp(2i*pi*szv*ori(3));
            
            sAxes = {sxv,syv,szv};
            
        end
        
        function [Gk,ind] = onePhononStructureFactors(obj,AtomFF,groupid,h,k,l)
            
            h = h(:);
            k = k(:);
            l = l(:);
            
            Basis = obj.Cell.Basis;
            Ops = obj.Cell.UnitCellOperators;
            if ~isa(Basis,'latt.OrientedBasis')
                Basis = Basis.orient;
            end
            
            nGroups = numel(obj.Cell.AsymmetricUnit);
            assert(all(groupid <= nGroups & groupid >= 1));
            
            nOps = numel(Ops);
            G = cell(nGroups,nOps);
            
            %[s0x,s0y,s0z] = Basis.invert.frac2lab(h,k,l);
            
            for g=1:nGroups
                A = AtomFF(groupid==g,:);
                [sAxes,ori,F,Fx,Fy,Fz] = obj.calcStructureFactorInterpolants(A);
                %pf0 = exp(2i*pi*(s0x*ori(1) + s0y*ori(2) + s0z*ori(3)));
                GI = griddedInterpolant(sAxes,F,obj.interp_mode);
                GIx = griddedInterpolant(sAxes,Fx,obj.interp_mode);
                GIy = griddedInterpolant(sAxes,Fy,obj.interp_mode);
                GIz = griddedInterpolant(sAxes,Fz,obj.interp_mode);
                
                for n=1:nOps
                    Op = Ops(n);
                    
                    % apply unit cell operator
                    hkln = [h,k,l,zeros(size(h))]*Op;
                    hn = hkln(:,1);
                    kn = hkln(:,2);
                    ln = hkln(:,3);
                    pn = hkln(:,4);
                    
                    pf = exp(2i*pi*pn);
                    
                    [sx,sy,sz] = Basis.invert.frac2lab(hn,kn,ln);
                    
                    pf0 = exp(2i*pi*(sx*ori(1) + sy*ori(2) + sz*ori(3)));
                    
                    qx = 2*pi*sx;
                    qy = 2*pi*sy;
                    qz = 2*pi*sz;
                    
                    Fn  = pf.*pf0.*GI(sx,sy,sz);
                    Fxn = pf.*pf0.*GIx(sx,sy,sz);
                    Fyn = pf.*pf0.*GIy(sx,sy,sz);
                    Fzn = pf.*pf0.*GIz(sx,sy,sz);
                    
                    % G{1,n} = [ F*qx, F*qy, F*qz, Fy*qz - Fz*qy, Fz*qx - Fx*qz, Fx*qy - Fy*qx]
                    G{g,n} = [Fn.*qx, Fn.*qy, Fn.*qz,...
                        Fyn.*qz - Fzn.*qy,...
                        Fzn.*qx - Fxn.*qz,...
                        Fxn.*qy - Fyn.*qx];
                    
                end
            end
            G = cat(2,G{:}); % flatten
            [Gk,ind] = obj.grouprowsbyk(h,k,l,G);
        end
        
        function I = simulateOnePhononScattering(obj,V,AtomFF,groupid,h,k,l)
            tic
            [Gk,ind] = obj.onePhononStructureFactors(AtomFF,groupid,h,k,l);
            toc
            tic
            Ik = obj.onePhononScattering(V,Gk);
            toc
            tic
            I = obj.kungroup(Ik,ind,size(h));
            toc
        end
        
        function [Gk,ind] = grouprowsbyk(h,k,l,G)
            G_bz = latt.PeriodicGrid(obj.supercell,[0,0,0],[1,1,1]).invert.invert;
            [n1,n2,n3] = G_bz.frac2ind(h(:),k(:),l(:));
            ind = accumarray([n1,n2,n3],1:numel(h),G_bz.N,@(v) {v});
            Gk = cell(obj.supercell);
            for j=1:numel(ind)
                Gk{j} = G(ind{j},:);
            end
        end
        
        function Ik = onePhononScattering(obj,V,Gk)
            LD = nm.LatticeDynamics('V',V,'supercell',obj.supercell);
            Kinv = LD.Kinv;
            Kinv = shiftdim(Kinv,3);
            Ik = cell(size(Gk));
            for j=1:numel(Ik)
                if isempty(Gk{j}), continue; end
                Ik{j} = real(dot(Gk{j},Gk{j}*conj(Kinv(:,:,j)),2));
            end
        end
        
    end
    
    methods(Static)
        
        
        function I = kungroup(Ik,ind,sz)
            I = zeros(sz);
            for j=1:numel(Ik)
                I(ind{j}) = Ik{j};
            end
        end
        

    end
    
end