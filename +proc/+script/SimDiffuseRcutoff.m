classdef SimDiffuseRcutoff < util.propertyValueConstructor
    %SimDiffuseRcutoff - compute diffuse patterson by "splat" method
    
    properties
        L
        Basis
        SpaceGroup
        Ops
        maxpairdist = 20
        Badd = 5
        dgrid = 0.3
        supercell = [1,1,1]
    end
    
    properties(Dependent) % for convenience
        B
    end
    
    methods
        function obj = SimDiffuseRcutoff(varargin)
            %SimDiffuseRcutoff
            obj@util.propertyValueConstructor(varargin{:});
        end
        function val = get.B(obj)
            val = obj.Basis.orthogonalizationMatrix;
        end
        
        function [SF,I_all,I_asu_intra,I_asu_inter] = run(obj,A,S)
            
            Ulatt = obj.calc_ADPs(A);

            [A.mdxFormFactor,T] = obj.factor_out_ADPs(A.mdxFormFactor,Ulatt);
            [S.mdxFormFactor] = obj.factor_out_ADPs(S.mdxFormFactor,Ulatt);

            [isym,iasu,covsym,covasu,An,Tn] = obj.precompute_neighbor_covariances(A,T);

            SF = obj.design_grid();

            [I_asu_intra,I_asu_inter] = obj.splat(SF,A,S,T,isym,iasu,covsym,covasu,An,Tn);
            
            % symmetrize
            fprintf(1,'SimDiffuseRcutoff.splat: symmetrizing...\n');
            
            SF2 = SF;
            SF2.SpaceGroup = SF2.SpaceGroup.PointGroup;
            I_all = real(SF2.apply_symmetry(I_asu_intra + I_asu_inter));
            
        end
        
        function [Ijw,Ijk] = splat(obj,SF,A,S,T,isym,iasu,covsym,covasu,An,Tn)
            
            nAtoms = size(A,1);
            Ijk = zeros(SF.LatticeGrid.N);
            Ijw = zeros(SF.LatticeGrid.N);
            
            P = obj.calc_projection_operator(A);

            for j = 1:nAtoms
                
                fprintf(1,'SimDiffuseRcutoff.splat: %d of %d\n',j,nAtoms);
                
                k = isym{j};
                
                if isempty(k)
                    Fk = 0;
                else
                    
                    Ank = An(k,:);
                    [fk_plus,fk_minus] = prep_deltapdf_formfactors(...
                        Ank.mdxFormFactor,P{j},T(j).U,{Tn(k).U},covsym(k));
                    
                    A1k = Ank;
                    A2k = Ank;
                    A1k.mdxFormFactor = fk_plus;
                    A2k.mdxFormFactor = fk_minus;
                    
                    % don't bother with neighbor solvent correlations, it's not a big contributor
                    rhok = SF.calc_electron_density([A1k;A2k]);
                    Fk = SF.LatticeGrid.ifft(rhok);
                end
                
                w = iasu{j};
                
                if isempty(w)
                    Fw = 0;
                else
                    
                    Aw = A(w,:);
                    [fw_plus,fw_minus] = prep_deltapdf_formfactors(...
                        Aw.mdxFormFactor,P{j},T(j).U,{T(w).U},covasu(w));
                    
                    A1w = Aw;
                    A2w = Aw;
                    A1w.mdxFormFactor = fw_plus;
                    A2w.mdxFormFactor = fw_minus;
                    
                    rhow = SF.calc_electron_density([A1w;A2w]);
                    Fw = SF.LatticeGrid.ifft(rhow);
                end
                
                rhoj = SF.calc_electron_density([A(j,:);S(j,:)]);
                Fj = SF.LatticeGrid.ifft(rhoj);
                
                Ijk = Ijk + conj(Fj).*Fk;
                Ijw = Ijw + conj(Fj).*Fw;
                
                assert(~any(isnan(Ijk(:))),'nan value!');
                assert(~any(isnan(Ijw(:))),'nan value!');
                
            end
            
            % sharpen if required
            if ~isempty(obj.Badd) && obj.Badd ~= 0
                fprintf(1,'SimDiffuseRcutoff.splat: sharpening\n');
            
                Fsharp = SF.calc_sharp_factor();
            
                Ijk = Ijk.*Fsharp.^2;
                Ijk = Ijk.*Fsharp.^2;
                
            end
            
            
        end
        
        function [isym,iasu,covsym,covasu,An,Tn] = precompute_neighbor_covariances(obj,A,T)
            % lattice neighbors
            [isym,Tsn] = proc.script.CoordinateTools.symmetry_neighbor_search(...
                A.x,A.y,A.z,obj.B,obj.Ops,obj.maxpairdist);

            % asu neighbors
            iasu = proc.script.CoordinateTools.asu_neighbor_search(...
                A.x,A.y,A.z,obj.maxpairdist,A.mdxAltGroup);

            P = obj.calc_projection_operator(A);
            covsym = calc_symmetry_neighbor_covariances(obj.L,P,Tsn,obj.Ops,obj.B);
            covasu = calc_asu_neighbor_covariances(obj.L,P);

            % Transform neighbors
            An = transform_symmetry_neighbors(A,Tsn,obj.Ops,obj.B);
            %Sn = transform_symmetry_neighbors(S,Tsn,obj.Ops,obj.B);

            % quick hack to transform the debye waller factor
            tmp = A;
            tmp.mdxFormFactor = T';
            tmp = transform_symmetry_neighbors(tmp,Tsn,obj.Ops,obj.B);
            Tn = tmp.mdxFormFactor;
            
        end
        
        function Ulatt = calc_ADPs(obj,A)
            
            P = obj.calc_projection_operator(A);
            
            covMat = obj.L.covarianceMatrix(obj.L.normalModes(),[0,0,0]);
            Ulatt = cellfun(@(p) p*covMat(1:6,1:6)*p',P,'Uni',0);

        end
        
        function P = calc_projection_operator(obj,A)
            nAtoms = size(A,1);
            
            P = nm.Group(A.x,A.y,A.z,1).tl2uxyz;
            P = mat2cell(P,3*ones(nAtoms,1),6);
        end
        
        function SF = design_grid(obj)
            
            G = proc.script.GridDesigner(...
                'Basis',obj.Basis,...
                'SpaceGroup',obj.SpaceGroup,...
                'dgrid',obj.dgrid).run();
            
            G.PeriodicGrid.N = G.PeriodicGrid.N.*obj.supercell;
            G.PeriodicGrid.P = obj.supercell;
            
            SF =proc.script.StructureFactor(...
                'LatticeGrid',G,...
                'SpaceGroup',obj.SpaceGroup,...
                'Badd',obj.Badd); % <-- play with this?
            
        end
    end
    
    methods (Static)
        function [FF,T] = factor_out_ADPs(FF,U)
            % subtract ADP from atomic scattering factor
            FF = FF.addU(cellfun(@(u) -u,U,'Uni',0));
            FF = FF.partitionUB();
            
            % if B is negative, convert to U isoU then try again
            isValid = [FF.B] >= 0;
            FF(~isValid) = FF(~isValid).toIsoU().partitionUB();
            
            % if B is still negative, just set B to zero
            isValid = [FF.B] >= 0;
            FF(~isValid) = FF(~isValid).addB(-[FF(~isValid).B]);
            
            % at this point, everybody should be OK
            isValid = [FF.B] >= 0;
            assert(all(isValid)); % hmmmm
            
            % create Debye-Waller factor from U
            T = latt.Blob.debyeWallerFactor(U);
            
            % to do: calculate residuals caused by the above corrections
        end
        


    end
end



function [fk_plus,fk_minus] = prep_deltapdf_formfactors(fk,Pj,Uj,Uk,covMatNk)
Vjk = cellfun(@(c) Pj*c,covMatNk,'Uni',0);

U1jk = cellfun(@(uk,vjk) uk + Uj - (vjk + vjk'),Uk(:),Vjk(:),'Uni',0);
U2jk = cellfun(@(uk) uk + Uj,Uk(:),'Uni',0);

fk_plus = fk.addU(U1jk);
fk_minus = fk.addU(U2jk).rescale(-1);
end


function covMatN = calc_asu_neighbor_covariances(L,P)

covMat = L.covarianceMatrix(L.normalModes,[0,0,0]);
covMat = covMat(1:6,1:6);

covMatN = cellfun(@(p) covMat*p',P,'Uni',0);
end


function covMatN = calc_symmetry_neighbor_covariances(L,P,Tsn,Ops,B)

[o,~,Tsn.ilattop] = unique(Tsn(:,{'xshift','yshift','zshift'}),'rows');
lattOps = rowfun(@(dx,dy,dz) symm.SymmetryOperator(eye(3),[dx,dy,dz]),o,'OutputFormat','uniform');

covMat = L.covarianceMatrix(L.normalModes,-[lattOps.t]');

nLattOps = numel(lattOps);
nOps = numel(Ops);

covMat = mat2cell(covMat,6*ones(nOps,1),6*ones(nOps,1),ones(nLattOps,1));

[ix,syms] = group_by_columns(Tsn,{'isymop','ilattop'});

Pn = P(Tsn.iatom); % transform the projection operator
covMatN = cell(size(Pn));

Ba = symm.AffineTransformation(B);
for k=1:size(ix)
    SymOp = lattOps(syms.ilattop(k))*Ops(syms.isymop(k));
    SymA = Ba*SymOp*inv(Ba);
    covMatk = covMat{1,syms.isymop(k),syms.ilattop(k)};
    covMatN(ix{k}) = cellfun(@(p) covMatk*(SymA.r*p)',Pn(ix{k}),'Uni',0);
end

end



function An = transform_symmetry_neighbors(A,Tsn,Ops,B)

[o,~,Tsn.ilattop] = unique(Tsn(:,{'xshift','yshift','zshift'}),'rows');
lattOps = rowfun(@(dx,dy,dz) symm.SymmetryOperator(eye(3),[dx,dy,dz]),o,'OutputFormat','uniform');

[ix,syms] = group_by_columns(Tsn,{'isymop','ilattop'});

An = A(Tsn.iatom,{'x','y','z','mdxFormFactor'});

for k=1:size(ix)
    SymOp = lattOps(syms.ilattop(k))*Ops(syms.isymop(k));
    An(ix{k},:) = proc.script.CoordinateTools.transform_atoms(...
        An(ix{k},:),SymOp,B);
end

end



function [ix,r] = group_by_columns(T,cols)

[r,~,ind] = unique(T(:,cols),'rows');
ix = accumarray(ind,1:size(T,1),[size(r,1),1],@(v) {v});

end
