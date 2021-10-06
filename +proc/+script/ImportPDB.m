classdef ImportPDB < util.propertyValueConstructor
    %IMPORTMODEL for reading PDB files
    
    properties % control what happens when you "run"
        addAltGroups = true
        addElementData = true
        includeHeteroAtoms = true
        includeAnisoTemps = true
        ignoreZeroOccupancyAtoms = true
        addResidueCategory = true
        occupancyTol = 0.01 + eps; % maximum difference in occupancy - used in alt group assessments
    end
    
    properties(Constant=true)
        % some values from NIST
        atomicMasses = [1.0079 4.0026 6.94 9.0122 10.811 12.0107 14.0067 15.9994 18.9984 20.18 22.9898 24.3051 26.9815 28.0855 30.9738 32.0648 35.4529 39.9478 39.0983 40.078 44.9559 47.8667 50.9415 51.9961 54.938 55.8451 58.9332 58.6933 63.546 65.3778 69.7231 72.6276 74.9216 78.9594 79.9035 83.798 85.4677 87.6166 88.9058 91.2236 92.9064 95.9598 98 101.0649 102.9055 106.4153 107.8681 112.4116 114.8181 118.7101 121.7598 127.6031 126.9045 131.2928 132.9055 137.3269 138.9055 140.1157 140.9077 144.2416 145 150.3664 151.9644 157.2521 158.9254 162.4995 164.9303 167.2591 168.9342 173.0542 174.9668 178.485 180.9479 183.8418 186.2067 190.2249 192.2161 195.0845 196.9666 200.5992 204.3834 207.2169 208.9804 209 210 222 223 226 227 232.0381 231.0359 238.0289 237 244 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
        atomicSymbols = {'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'};
        vdwRadii = [1.2 1.4 1.82 0.63 1.75 1.775 1.5 1.45 1.47 1.54 2.27 1.73 1.5 2.1 1.9 1.8 1.75 1.88 2.75 1.95 1.32 1.95 1.06 1.13 1.19 1.26 1.13 1.63 1.4 1.39 1.87 1.48 0.83 1.9 1.85 2.02 2.65 2.02 1.61 1.42 1.33 1.75 2 1.2 1.22 1.63 1.72 1.58 1.93 2.17 1.12 1.26 1.98 2.16 3.01 2.41 1.83 1.86 1.62 1.79 1.76 1.74 1.96 1.69 1.66 1.63 1.61 1.59 1.57 1.54 1.53 1.4 1.22 1.26 1.3 1.58 1.22 1.72 1.66 1.55 1.96 2.02 1.73 1.21 1.12 2.3 3.24 2.57 2.12 1.84 1.6 1.75 1.71 1.67 1.66 1.65 1.64 1.63 1.62 1.61 1.6 1.59 1.58 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
        
        % These atoms are usually ions when observed in PDB files... except
        % when part of metal clusters. Such cases may need to be dealt with
        % manually.
        commonIons = {'Na','Mg','Al','Cl','K','Ca','Cr','Mn','Fe','Ni','Cu','Zn','Br','Rb','Sr','Cs','Ba','I'};
        
        % https://proteopedia.org/wiki/index.php/Standard_Residues
        standardAminoAcids = {...
            'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU',...
            'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR',...
            'MSE', 'ASX', 'GLX', 'UNK', 'SEC', 'PYL'}; % <-- last row is less standard
        
        standardNucleicAcids = {...
            'A', 'C', 'G', 'I', 'U','DA', 'DC', 'DG', 'DI', 'DT', 'DU'};
        
        waterNames = {'HOH','WAT'};
        
    end
    
    methods
        function obj = ImportPDB(varargin)
            %IMPORTMODEL
            obj@util.propertyValueConstructor(varargin{:});
        end
        
        function [Atoms,Basis,SpaceGroup] = run(obj,fileName)
            
            [Crystal,Atoms,HetAtoms,AnisoTemps,headerLines]  = io.pdb.read(fileName);
            Atoms.isHet = false(size(Atoms,1),1);
            
            if obj.includeHeteroAtoms
                HetAtoms.isHet = true(size(HetAtoms,1),1);
                Atoms = [Atoms;HetAtoms];
            end
            
            if obj.includeAnisoTemps
                [ia,locb] = ismember(AnisoTemps.serial,Atoms.serial);
                assert(all(ia),'some atom records are missing?');
                
                Atoms = [Atoms, array2table(NaN*ones(size(Atoms,1),6),'VariableNames',{'u11','u22','u33','u12','u13','u23'})];
                Atoms.u11(locb) = AnisoTemps.u00*1E-4;
                Atoms.u22(locb) = AnisoTemps.u11*1E-4;
                Atoms.u33(locb) = AnisoTemps.u22*1E-4;
                Atoms.u12(locb) = AnisoTemps.u01*1E-4;
                Atoms.u13(locb) = AnisoTemps.u02*1E-4;
                Atoms.u23(locb) = AnisoTemps.u12*1E-4;
            end
            
            if obj.ignoreZeroOccupancyAtoms
                Atoms = Atoms(Atoms.occupancy~=0,:);
            end
            
            % include unit cell info in CustomProperties, in case they are
            % needed later
            Atoms = addprop(Atoms,{'Crystal','Remarks'},{'table','table'});
            Atoms.Properties.CustomProperties.Crystal = Crystal;
            Atoms.Properties.CustomProperties.Remarks = strjoin(headerLines,'\n');
            
            % if requested, go ahead and generage the mdx-lib objects to
            % represent the unit cell
            if nargout >= 2
                Basis = latt.Basis(Crystal.a,Crystal.b,Crystal.c,...
                    Crystal.alpha,Crystal.beta,Crystal.gamma);
            end
            if nargout >= 3
                SpaceGroup = symm.SpaceGroup(strrep(Crystal.sGroup,' ',''));
            end
            
            if obj.addElementData
                Atoms = obj.assign_element_properties(Atoms);
            end
            
            if obj.addAltGroups
                Atoms = obj.assign_alternate_groups(Atoms);
            end
            
            if obj.addResidueCategory
                Atoms = obj.assign_residue_category(Atoms);
            end
            
        end
        
        
        function Atoms = assign_element_properties(obj,Atoms)
            
            % go row by row... inefficient but only need to do it once
            
            mdxAtomicSymbol = cell(size(Atoms,1),1);
            mdxVdwRadius = zeros(size(Atoms,1),1);
            mdxAtomicMass = zeros(size(Atoms,1),1);
            mdxIsIon = false(size(Atoms,1),1);
            
            for j=1:size(Atoms,1)
                ind = find(strcmpi(Atoms.element(j),obj.atomicSymbols),1,'first');
                assert(~isempty(ind),'element not found!');
                
                mdxAtomicSymbol{j} = obj.atomicSymbols{ind};
                mdxVdwRadius(j) = obj.vdwRadii(ind);
                mdxAtomicMass(j) = obj.atomicMasses(ind);
                if ~isempty(Atoms.charge{j}) || ismember(obj.atomicSymbols(ind),obj.commonIons)
                    mdxIsIon(j) = true;
                end
            end
            
            Atoms.mdxAtomicSymbol = mdxAtomicSymbol;
            Atoms.mdxAtomicMass = mdxAtomicMass;
            Atoms.mdxVdwRadius = mdxVdwRadius;
            Atoms.mdxIsIon = mdxIsIon;
            
        end
        
        function Atoms = assign_alternate_groups(obj,Atoms)
            
            % an "atom" can have multiple records corresopnding to different alternate
            % conformers. I think (?) that that chainID, resSeq, iCode, and name are
            % sufficient to identify a unique atom.
            [~,~,Atoms.mdxAtomID] = unique(Atoms(:,{'chainID','resSeq','iCode','name'}),'rows');
            
            % for convenience, convert string altLoc into a numeric index
            [altCode,~,altID] = unique(Atoms.altLoc);
            
            % in PDB files, atoms with no alt conformers have a blank altLoc "". This
            % should always be sorted to the first position by unique(). But good to
            % check:
            assert(isempty(altCode{1}),'assume altCode = "" is first column');
            altCode = altCode(2:end); % drop the first alt code
            
            % to analyze alt confs, let's make an array where each row is a particular
            % atom, and each column contains the row index of that atom's alt
            % conformers in the original Atoms table (or zero)
            altgroup_ind = accumarray([Atoms.mdxAtomID,altID],1:size(Atoms,1),[],@max);
            
            % an array with the same format, but with occupancies instead
            altgroup_occ = accumarray([Atoms.mdxAtomID,altID],Atoms.occupancy,[],@sum);
            
            % There are two kinds of atoms to consider. For polymers (protein) we can
            % assume that all (atoms should have occupancy of 1 if you sum over the
            % alternative conformers. For heteroatoms (water, ions) the occupancy might
            % be less than one, which is fine.
            %
            % Let's check the first scenario:
            
            isHet = accumarray(Atoms.mdxAtomID,Atoms.isHet,[],@any);
            totalOccupancy = sum(altgroup_occ,2);
            numConfs = sum(logical(altgroup_ind),2);
            
            assert(all(abs(totalOccupancy(~isHet)-1) <= obj.occupancyTol),...
                'occupancies should add to 1 for non-hetero atoms');
            
            % we also need to make sure that, if the "blank" alternate conformer is
            % present, there are no "A","B" etc.
            %hasAlternates = numConfs > 1;
            %assert(~hasAlternates & logical(altgroup_ind(1,:)),'this should never happen');
            
            % now, assign each atom record to one or more alternate conformer groups
            numAlts = numel(altCode);
            mdxAltGroup = false([size(Atoms,1),numAlts]);
            
            for j=1:size(altgroup_ind,1)
                K = find(altgroup_ind(j,:));
                rowIndices = altgroup_ind(j,K); % row indices in the Atoms table
                hasDefaultConformer = logical(altgroup_ind(j,1));
                
                if isHet(j)
                    if hasDefaultConformer
                        mdxAltGroup(rowIndices,:) = true;
                    else
                        ind = sub2ind(size(mdxAltGroup),rowIndices,K-1);
                        mdxAltGroup(ind) = true;
                    end
                else % not isHet(j)
                    if hasDefaultConformer
                        mdxAltGroup(rowIndices,:) = true;
                    elseif numel(K) == numAlts % all alts are assigned
                        ind = sub2ind(size(mdxAltGroup),rowIndices,K-1);
                        mdxAltGroup(ind) = true;
                    else % ambiguous situation
                        % say there are 3 alt confs A, B, C.
                        % atom 1 is assigned A and B, but not C.
                        % what did the modeler intend to happen with atom 1 in state C?
                        % Probably it does not disappear. Probably it exists in either
                        % A or B, or both A and B. For simplicity, here we choose
                        % whichever one has higher occupancy. It would be BETTER to
                        % check for clashes and therefore infer the modeler's intent,
                        % but that's beyond the scope here. Just issue a warning.
                        
                        [~,default_conf] = max(altgroup_occ(j,2:end));
                        for k=1:numAlts
                            if logical(altgroup_ind(j,k+1))
                                rowIndex = altgroup_ind(j,k+1);
                                mdxAltGroup(rowIndex,k) = true;
                            else
                                rowIndex = altgroup_ind(j,default_conf + 1);
                                warning(strjoin({...
                                    'Ambiguous conformer assignment for atom %d',...
                                    '  - assuming state %d is occupied by altLoc %s'...
                                    },'\n'),rowIndex,k,altCode{default_conf});
                                mdxAltGroup(rowIndex,k) = true;
                            end
                        end
                    end
                end
                
            end
            
            Atoms.mdxAltGroup = mdxAltGroup;
            
        end
        
        function A = assign_residue_category(obj,A)
            isAA = ~A.isHet & ismember(A.resName,obj.standardAminoAcids);
            isNA = ~A.isHet & ismember(A.resName,obj.standardNucleicAcids);
            
            isWater = A.isHet & ismember(A.resName,obj.waterNames);
            isLigand = A.isHet & ~isWater;
            
            aminoAcidChains = unique(A.chainID(isAA));
            nucleicAcidChains = unique(A.chainID(isNA));
            
            isAA = ~A.isHet & ismember(A.chainID,aminoAcidChains);
            isNA = ~A.isHet & ismember(A.chainID,nucleicAcidChains);
            
            %isOther = ~isAA & ~isNA & ~isLigand & ~isWater; % ideally this should be nothing
            
            % make sure each atom has only one assignment
            assert(all(sum([isAA,isNA,isWater,isLigand],2)<=1));
            
            valueset = {'protein' 'nucleic acid','ligand','water'};
            A.mdxResidueCategory = categorical(repmat({''},size(A,1),1),valueset);
            A.mdxResidueCategory(isAA) = {'protein'};
            A.mdxResidueCategory(isNA) = {'nucleic acid'};
            A.mdxResidueCategory(isLigand) = {'ligand'};
            A.mdxResidueCategory(isWater) = {'water'};
            
            % now, assign backbone / sidechain for protein atoms
            isBackbone = ismember(A.name,{'C','N','O','CA'}) & A.mdxResidueCategory == 'protein';
            isSidechain = ~isBackbone & A.mdxResidueCategory == 'protein';
            
            valueset = {'backbone' 'sidechain','nucleobase'}; % haven't implemented base / sugar for DNA/RNA yet
            A.mdxChemicalGroup = categorical(repmat({''},size(A,1),1),valueset);
            A.mdxChemicalGroup(isBackbone) = {'backbone'};
            A.mdxChemicalGroup(isSidechain) = {'sidechain'};
            
        end
        
        
    end
    
    methods(Static = true)
        function AltAtoms = splitAltGroups(Atoms)
            % split into alt groups. if an entry is present in multiple alt
            % groups, the occupancy column is split evenly
            
            numAlt = size(Atoms.mdxAltGroup,2);
            
            AltAtoms = cell(1,numAlt);
            
            Atoms.occupancy = rowfun(@(o,ag) o/nnz(ag),Atoms(:,{'occupancy','mdxAltGroup'}),'OutputFormat','Uniform');
            
            for j=1:numAlt
                AltAtoms{j} = removevars(Atoms(Atoms.mdxAltGroup(:,j),:),'mdxAltGroup');
            end
        end
        
        function referenceID = assign_by_proximity(Atoms,referenceID)
            % referenceID is a number labeling each group
            % and it is zero if the atom is unassigned.
            %
            % assign_by_proximity finds the index of nearest atom in the
            % reference group, then assigns the atom its reference ID.
            %
            
            rowIndices = find(~logical(referenceID));
            %ind = find_nearest_neighbors(Atoms,unassignedRows);

            A = Atoms(:,{'x','y','z','mdxAltGroup'});
            A.isIncluded = false(size(A,1),1);
            A.isIncluded(rowIndices) = true;
            A.row = (1:size(A,1))';
            
            numAlt = size(A.mdxAltGroup,2);
            
            ind = zeros(size(A,1),1);

            for n=1:numAlt
                An = A(A.mdxAltGroup(:,n),:);
                AnIncl = An(An.isIncluded,:);
                AnRem = An(~An.isIncluded,:);
                xyzIncl = table2array(AnIncl(:,{'x','y','z'}));
                xyzRem = table2array(AnRem(:,{'x','y','z'}));
                ixRem = knnsearch(xyzRem,xyzIncl);
                ind(AnIncl.row) = AnRem.row(ixRem);
            end
            
            ind = ind(rowIndices);
            
            referenceID(rowIndices) = referenceID(ind);
        end
        
        function Scatterers = to_xray_structure(Atoms)
            % Convert each line in the Atoms table to an X-ray scattering
            % factor.
            %
            % Fields x, y, z, and occupancy are carried through.
            %
            % If u11, u22, ... exist then tempFactor is ignored
            %
            % Otherwise, an isotropic B-factor is assigned using tempFactor
            %
            % A scattering factor object (latt.Blob) is created for each row.
            %
            % If the column "mdxAtomicSymbol" exists, it is used for
            % scattering factor lookup. Otherwise, the "element" column is
            % used. If neither column exists, the code checks if a "v"
            % column exists, which would indicate a blobbified solvent
            % mask. An ~uniform ellipsoidal form factor is used in this
            % case.
                        
            if ismember('mdxAtomicSymbol',Atoms.Properties.VariableNames)
                % use mdxAtomicSymbol for form factor lookup
                structureFactorType = 'mdxAtomicSymbol';
            elseif ismember('element',Atoms.Properties.VariableNames)
                % use element for form factor lookup
                structureFactorType = 'element';
            else
                error('should not happen');
            end
            
            hasUaniso = all(ismember({'u11','u22','u33','u12','u13','u23'},Atoms.Properties.VariableNames));
            hasBiso = ismember('tempFactor',Atoms.Properties.VariableNames);
            hasOccupancy = ismember('occupancy',Atoms.Properties.VariableNames);
            
            assert(hasOccupancy)
            
            if hasUaniso
                params = {structureFactorType,'occupancy','u11','u22','u33','u12','u13','u23'};
                fffun = @(el,o,u11,u22,u33,u12,u13,u23) latt.Blob.atomicScatteringFactor(el).rescale(o).addU([u11,u12,u13;u12,u22,u23;u13,u23,u33]);
            elseif hasBiso
                params = {structureFactorType,'occupancy','tempFactor'};
                fffun = @(el,o,B) latt.Blob.atomicScatteringFactor(el).rescale(o).addB(B);
            else
                error('should not happen');
            end
            
            varsToKeep = intersect({'mdxAtomID','mdxAltGroup','x','y','z'},Atoms.Properties.VariableNames,'stable');
            Scatterers = Atoms(:,varsToKeep);
            
            Scatterers.mdxFormFactor = rowfun(fffun,Atoms(:,params),'OutputFormat','uniform');
        end
        
    end
    
end

