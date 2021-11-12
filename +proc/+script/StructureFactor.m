classdef StructureFactor < util.propertyValueConstructor
    %StructureFactor - compute using fast FFT-based method
    
    properties
        LatticeGrid
        SpaceGroup
        Badd
        rmax = 3.0
        hklTable
    end
    
    methods
        function obj = StructureFactor(varargin)
            %StructureFactor
            obj@util.propertyValueConstructor(varargin{:});
        end
        
        function F = run(obj,AtomFF)
            rho = obj.calc_electron_density(AtomFF);
            F = obj.calc_structure_factor(rho);
            if ~isempty(obj.SpaceGroup)
                F = obj.apply_symmetry(F);
            end
            if ~isempty(obj.hklTable)
            	F = obj.get_at_hkl(F);
            end
        end
        
        function rho = calc_electron_density(obj,AtomFF)
            ff = AtomFF.mdxFormFactor;
            if ~isempty(obj.Badd)
                ff = ff.addB(obj.Badd);
            end
            [k1,k2,k3] = obj.LatticeGrid.sphkernel(obj.rmax);
            
            % calculate electron density (splat algorithm)
            rhoFun = @(X,Y,Z,B) electronDensity(B,X,Y,Z);
            
            rho = obj.LatticeGrid.splat(k1,k2,k3,rhoFun,@plus,0,AtomFF.x,AtomFF.y,AtomFF.z,ff);
        end
        
        function F = calc_structure_factor(obj,rho)
            
            F = obj.LatticeGrid.ifft(rho);
            
            if ~isempty(obj.Badd) % need to sharpen
                F = F.*obj.calc_sharp_factor();
            end
            
        end
        
        function F = calc_sharp_factor(obj)
            
            Bsharp = latt.Blob(1,0).addB(-obj.Badd);
            [sx,sy,sz] = obj.LatticeGrid.invert().grid();
            smax = proc.script.GridDesigner.slimit(obj.LatticeGrid);
            F = Bsharp.scatteringAmplitude(sx,sy,sz);
            F(sx.^2 + sy.^2 + sz.^2 > smax^2) = 0;
            
        end
        
        function hklTable = create_hkl_table(obj,smax)
            R = obj.LatticeGrid.invert;
            [h,k,l] = R.PeriodicGrid.grid();
            [sx,sy,sz] = R.Basis.frac2lab(h,k,l);
            if nargin < 2 || isempty(smax)
                h = h(:); k = k(:); l = l(:);
            else
                isIncl = sx.^2 + sy.^2 + sz.^2 <= smax^2;
                h = h(isIncl); k = k(isIncl); l = l(isIncl);
            end
            hklTable = table(h,k,l);
        end
        
        function Fout = get_at_hkl(obj,F)
            hkl = table2array(obj.hklTable(:,{'h','k','l'}));
            R = obj.LatticeGrid.invert();
            [n1,n2,n3] = R.PeriodicGrid.frac2ind(hkl(:,1),hkl(:,2),hkl(:,3),false);
            isIncl = n1>0 & n1<=R.N(1) & n2>0 & n2<=R.N(2) & n3>0 & n3<=R.N(3);
            ind = sub2ind(R.N,n1(isIncl),n2(isIncl),n3(isIncl));
            Fout = NaN*ones(size(hkl,1),1);
            Fout(isIncl) = F(ind);
        end
        
        function Fcell = apply_symmetry(obj,F)
            % expand to unit cell
            R = obj.LatticeGrid.invert();
        
            [h,k,l] = R.PeriodicGrid.grid();
            Ops = obj.SpaceGroup.generalPositions;
            hklp0 = [h(:),k(:),l(:),zeros(numel(h),1)];
            Fcell = zeros(R.N);
            for n=1:numel(Ops)
                hklp = hklp0*Ops(n);
                phaseFactor = exp(2i*pi*hklp(:,4));
                [n1,n2,n3] = R.PeriodicGrid.frac2ind(hklp(:,1),hklp(:,2),hklp(:,3),false);
                isIncl = n1>0 & n1<=R.N(1) & n2>0 & n2<=R.N(2) & n3>0 & n3<=R.N(3);
                ind = sub2ind(R.N,n1(isIncl),n2(isIncl),n3(isIncl));
                Fcell(isIncl) = Fcell(isIncl) + phaseFactor(isIncl).*F(ind);
            end
        end
        
        
        
    end
    
end

