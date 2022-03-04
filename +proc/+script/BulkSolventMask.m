classdef BulkSolventMask < util.propertyValueConstructor
    %BULKSOLVENTMASK
    
    properties
        vdwProb = 1.2;
        ionProb = 0.8;
        rShrink = 0.8;
        rMax = []; % if empty, make a guess based on prob radius etc.
        LatticeGrid = latt.LatticeGrid.empty(); % lattice grid object
        SpaceGroup = symm.SpaceGroup.empty();
        Atoms % table of atoms. Required fields are x,y,z,mdxVdwRadius,mdxIsIon,mdxAltGroup
    end
    
    methods
        function obj = BulkSolventMask(varargin)
            %BULKSOLVENTMASK
            obj@util.propertyValueConstructor(varargin{:});
        end
        
        function [indMaps] = run(obj)
            % number of alternate conformers
            numAlt = size(obj.Atoms.mdxAltGroup,2);
            
            indMaps = cell(1,numAlt);
            solvExclMask = cell(1,numAlt);
            
            for n=1:numAlt
                isIncl = obj.Atoms.mdxAltGroup(:,n);
                [indMaps{n},probDist] = obj.distance_map(isIncl);
                
                % compute solvent accessibility (slow because of slow symmetry
                % expansion step for some reason...)
                solvExclMask{n} = obj.solvent_excluded_region(probDist);
                
                indMaps{n}(~solvExclMask{n}) = 0;
            end
            
            isConsensus = obj.consensus_map(indMaps);
            indMapConsensus = indMaps{1};
            indMapConsensus(~isConsensus) = 0;
            probDistConsensus = probDist;
            probDistConsensus(~isConsensus) = Inf;
            
            isConsensusAtom = ismember((1:size(obj.Atoms,1))',indMapConsensus);
            
            for n=1:numAlt
                isIncl = obj.Atoms.mdxAltGroup(:,n) & ~isConsensusAtom;
                if ~any(isIncl)
                    continue;
                end
                indMaps{n} = obj.distance_map_update(isIncl,indMapConsensus,probDistConsensus);
                indMaps{n}(~solvExclMask{n}) = 0;
            end
            
        end
%         
        function [Sol,bloblist] = blobify(obj,indMaps,binH)
            
            if nargin < 3 || isempty(binH)
                binH = false;
            end
            
            bloblist = obj.solvent_mask_blobs(indMaps);
            
            if binH
                hInd =find(strcmp(obj.Atoms.mdxAtomicSymbol,'H'));
                hMap = obj.find_nearest_neighbors(hInd);
                
                for j=1:numel(hInd)
                    bloblist{hMap(j)} = [bloblist{hMap(j)}; bloblist{hInd(j)}];
                    bloblist{hInd(j)} = [];
                end

            end
            
            % for each item in bloblist, compute various parameters:
            % x,y,z,u11,u22,u33,u12,u13,u23,v
            %
            % also include (from Atoms):
            % mdxAtomID, mdxAltGroup
            
            hasBlob = cellfun(@(r) size(r,1) > 0,bloblist);
            
            v = det(obj.LatticeGrid.delta)*cellfun(@(r) size(r,1),bloblist(hasBlob));
            cen = cellfun(@(r) mean(r,1),bloblist(hasBlob),'Uni',0);
            U = cellfun(@(r,c) (r-c)'*(r-c)/size(r,1),bloblist(hasBlob),cen,'Uni',0);
            
            u11 = cellfun(@(u) u(1,1),U);
            u22 = cellfun(@(u) u(2,2),U);
            u33 = cellfun(@(u) u(3,3),U);
            u12 = cellfun(@(u) u(1,2),U);
            u13 = cellfun(@(u) u(1,3),U);
            u23 = cellfun(@(u) u(1,3),U);
            
            x = cellfun(@(c) c(1),cen);
            y = cellfun(@(c) c(2),cen);
            z = cellfun(@(c) c(3),cen);
            
            Sol = [table(x,y,z,v,u11,u22,u33,u12,u13,u23),obj.Atoms(hasBlob,{'mdxAtomID','mdxAltGroup'})];
        end
        
        function S = to_xray_structure(obj,Sol)
            warning('deprecated'); % <-- will move to another class
           
            %S = Sol(:,{'x','y','z'});
            
            params = {'v','u11','u22','u33','u12','u13','u23'};
            Sol.mdxFormFactor = rowfun(...
                @(v,u11,u22,u33,u12,u13,u23) latt.Blob.betterEllipsoid([u11,u12,u13;u12,u22,u23;u13,u23,u33]).rescale(v),...
                Sol(:,params),...
                'OutputFormat','uniform');
            Sol = removevars(Sol,params);
            
            % apply occupancies and Uaniso from parent atoms to solvent
            % mask blobs

            hasUaniso = all(ismember({'u11','u22','u33','u12','u13','u23'},obj.Atoms.Properties.VariableNames));
            hasBiso = ismember('tempFactor',obj.Atoms.Properties.VariableNames);
            
            uid_keys = {'mdxAtomID','mdxAltGroup'};
            
            if ~hasUaniso
                assert(hasBiso)
                params = {'occupancy','tempFactor'};
                fffun = @(ff,o,b) ff.rescale(o).addB(b);
            else
                params = {'occupancy','u11','u22','u33','u12','u13','u23'};
                fffun = @(ff,o,u11,u22,u33,u12,u13,u23) ff.rescale(o).addU([u11,u12,u13;u12,u22,u23;u13,u23,u33]);
            end
            
            % look up occupancy etc from parent structure
            [~,ia] = ismember(Sol(:,uid_keys),obj.Atoms(:,uid_keys),'rows');
            Sol = [Sol,obj.Atoms(ia,params)];

            Sol.mdxFormFactor = rowfun(fffun,...
                Sol(:,[{'mdxFormFactor'},params]),...
                'OutputFormat','uniform');
            
            S = Sol(:,[uid_keys,{'x','y','z','mdxFormFactor'}]);
        end
        
    end
    methods(Hidden = true)
        
        function ind = find_nearest_neighbors(obj,rowIndices)
            A = obj.Atoms(:,{'x','y','z','mdxAltGroup'});
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
        end
        
        function [indMapOut,probDistOut] = distance_map_update(obj,isIncluded,indMap,probDist)
            [a,b] = obj.distance_map(isIncluded);
            isAltered = b<probDist;
            indMapOut = indMap;
            probDistOut = probDist;
            probDistOut(isAltered) = b(isAltered);
            indMapOut(isAltered) = a(isAltered);
        end
        
        function [indMapOut,probDist] = distance_map(obj,isIncluded)
            
            if nargin < 2 || isempty(isIncluded)
                isIncluded = true(size(obj.Atoms,1),1);
            end
            
            A0 = obj.Atoms(isIncluded,{'x','y','z','mdxVdwRadius','mdxIsIon'});
            
            % add probe radius
            A0.r = A0.mdxVdwRadius + (~A0.mdxIsIon)*obj.vdwProb + (A0.mdxIsIon)*obj.ionProb;
            
            if isempty(obj.rMax)
                rKer = max(A0.r) + max(sqrt(sum(obj.LatticeGrid.delta.^2,1))); % maximum radius to calculate
            else
                rKer = obj.rMax;
            end
            
            % define the distance function
            distfun = @(xx,yy,zz,rr) sqrt(xx.*xx + yy.*yy + zz.*zz) - rr;
            
            [k1,k2,k3] = obj.LatticeGrid.sphkernel(rKer);
            
            [probDist,indMap] = obj.LatticeGrid.splat(k1,k2,k3,distfun,@le,Inf,A0.x,A0.y,A0.z,A0.r);
            
            % take care of atom selection - look up indices of original input array
            indMapOut = zeros(size(indMap));
            indOut = find(isIncluded);
            indMapOut(logical(indMap)) = indOut(indMap(logical(indMap)));
            
        end
        
        function [isSolventExcluded,isAccessible,isASU] = solvent_excluded_region(obj,probDistASU)
            
            neighborOps = setdiff(obj.SpaceGroup.generalPositions,symm.SymmetryOperator());
            
            tmp = obj.symmetrize_probdist(probDistASU,neighborOps);
            
            % compute
            isAccessible = obj.shrink_mask(tmp > 0 & probDistASU > 0);
            isASU = probDistASU < tmp;
            
            isSolventExcluded = ~isAccessible & isASU;
            
        end
        
        function isAccessible = shrink_mask(obj,isAccessible)
            
            % design a kernel for dilating the map
            [k1,k2,k3] = obj.LatticeGrid.sphkernel(obj.rShrink);
            n1 = max(k1); n2 = max(k2); n3 = max(k3);
            nhood = false(n1*2 + 1,n2*2 + 1,n3*2 + 1);
            ind = sub2ind(size(nhood),k1 + n1 + 1,k2 + n2 + 1,k3 + n3 + 1);
            nhood(ind) = true;
            
            % dilate the map
            isAccessible = imdilate(isAccessible,nhood);
            
        end

        
        function newMap = symmetrize_probdist(obj,oldMap,SymOps)
            
            newMap = inf*ones(size(oldMap));
            
            for j=1:numel(SymOps)
                tmp = obj.transform_map(oldMap,SymOps(j));
                isAltered = tmp < newMap;
                newMap(isAltered) = tmp(isAltered);
            end

        end
        
        function newMap = transform_map(obj,oldMap,SymOp)
            
            [f1,f2,f3] = obj.LatticeGrid.PeriodicGrid.grid();
            f123t = SymOp*[f1(:),f2(:),f3(:)]';
            [n1,n2,n3] = obj.LatticeGrid.PeriodicGrid.frac2ind(f123t(1,:),f123t(2,:),f123t(3,:));
            ind = reshape(sub2ind(obj.LatticeGrid.N,n1,n2,n3),obj.LatticeGrid.N);
            
            newMap = oldMap(ind);
            
        end
        
        function bloblist = solvent_mask_blobs(obj,indMaps)
            
            nAtoms = size(obj.Atoms,1);
            
            gridIndex = 1:numel(indMaps{1});
            bloblist = repmat({[]},[nAtoms,1]);
            
            for k=1:numel(indMaps)
                isIncl = logical(indMaps{k});
                blobn = accumarray(indMaps{k}(isIncl),gridIndex(isIncl),[nAtoms,1],@(d) {d});
                for n=1:nAtoms
                    bloblist{n} = unique([bloblist{n};blobn{n}]);
                end
            end
            
            % now, convert blobs to coordinates
            for j=1:numel(bloblist)
                [n1,n2,n3] = ind2sub(obj.LatticeGrid.N,bloblist{j});
                [f1,f2,f3] = obj.LatticeGrid.PeriodicGrid.ind2frac(n1,n2,n3);
                [fa1,fa2,fa3] = obj.LatticeGrid.Basis.lab2frac(obj.Atoms.x(j),obj.Atoms.y(j),obj.Atoms.z(j));
                f1 = mod(f1-fa1 + 0.5,1) + fa1 - 0.5;
                f2 = mod(f2-fa2 + 0.5,1) + fa2 - 0.5;
                f3 = mod(f3-fa3 + 0.5,1) + fa3 - 0.5;
                [x,y,z] = obj.LatticeGrid.Basis.frac2lab(f1,f2,f3);
                bloblist{j} = [x,y,z];
            end
            
        end

    end
    methods(Static = true, Hidden = true)
        
                
        function consensusMask = consensus_map(indMaps)
            
            consensusMask = true(size(indMaps{1}));
            indMap = indMaps{1};
            
            for j=2:numel(indMaps)
                consensusMask = consensusMask & indMap == indMaps{j};
                indMap(~consensusMask) = 0;
            end
            
        end
    end
end
