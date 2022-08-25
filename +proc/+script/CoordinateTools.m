classdef CoordinateTools
    %COORDINATETOOLS
    
    properties
    end
    
    %     methods
    %         function obj = CoordinateTools()
    %             %COORDINATETOOLS
    %         end
    %     end
    
    methods(Static = true)
        
        
        function T = transform_atoms(T,A,B)
            
            if nargin == 3
                % assume that A is a symmetry operator and B is either a latt.Basis
                % object or an orthogonalization matrix
                assert(isa(A,'symm.SymmetryOperator'));
                if isa(B,'latt.Basis')
                    B = B.orthogonalizationMatrix;
                end
                B = symm.AffineTransformation(B,[0;0;0]);
                A = B*A*inv(B);
            end
            
            % next, transform the adp components: U <- R*U*R'
            if isnumeric(A) % assume it's a 3x3 matrix
                R = A;
            else
                R = A.r; % thake the rotational part of the symm.AffineTransformation
            end
            
            % first, transform the coordinates (x,y,z)
            xyz = table2array(T(:,{'x','y','z'}));
            xyz = (A*(xyz'))';
            T.x = xyz(:,1);
            T.y = xyz(:,2);
            T.z = xyz(:,3);
            
            if ismember('u11',T.Properties.VariableNames)
                % transform adps also
                % fprintf(1,'transforming u11, u22, ...\n');
                
                [T.u11,T.u22,T.u33,T.u12,T.u13,T.u23] = transform_adp(R,...
                    T.u11,T.u22,T.u33,T.u12,T.u13,T.u23);
                
            end
            
            if ismember('mdxFormFactor',T.Properties.VariableNames)
                % transform ADP component of form factor also
                fprintf(1,'transforming mdxFormFactor\n');
                T.mdxFormFactor = T.mdxFormFactor.transform(R);
            end
            
        end
        
        
        function [AtomFF,SolFF] = to_xray_structure(Atoms,Sol)
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
            % used.
            
            if ismember('mdxAtomicSymbol',Atoms.Properties.VariableNames)
                % use mdxAtomicSymbol for form factor lookup
                structureFactorType = 'mdxAtomicSymbol';
            elseif ismember('element',Atoms.Properties.VariableNames)
                % use element for form factor lookup
                structureFactorType = 'element';
            else
                error('should not happen');
            end
            
            uid_keys = {'mdxAtomID','mdxAltGroup'}; % <- required for now
            assert(all(ismember(uid_keys,Atoms.Properties.VariableNames)));
            
            hasUaniso = all(ismember({'u11','u22','u33','u12','u13','u23'},Atoms.Properties.VariableNames));
            hasBiso = ismember('tempFactor',Atoms.Properties.VariableNames);
            hasOccupancy = ismember('occupancy',Atoms.Properties.VariableNames);
            
            assert(hasOccupancy) % required for now
            
            if hasUaniso
                params = {structureFactorType,'occupancy','u11','u22','u33','u12','u13','u23'};
                fffun = @(el,o,u11,u22,u33,u12,u13,u23) latt.Blob.atomicScatteringFactor(el).rescale(o).addU([u11,u12,u13;u12,u22,u23;u13,u23,u33]);
            elseif hasBiso
                params = {structureFactorType,'occupancy','tempFactor'};
                fffun = @(el,o,B) latt.Blob.atomicScatteringFactor(el).rescale(o).addB(B);
            else
                error('should not happen');
            end
            
            AtomFF = Atoms(:,[uid_keys,{'x','y','z'}]);
            
            AtomFF.mdxFormFactor = rowfun(fffun,Atoms(:,params),'OutputFormat','uniform');
            
            
            if nargin < 2
                return; % done
            end
            
            % now, deal with Sol
            params = {'v','u11','u22','u33','u12','u13','u23'};
            
            required_columns = [uid_keys,{'x','y','z'},params];
            assert(all(ismember(required_columns,Sol.Properties.VariableNames)));
            
            SolFF = Sol(:,[uid_keys,{'x','y','z'}]);
            
            SolFF.mdxFormFactor = rowfun(...
                @(v,u11,u22,u33,u12,u13,u23) latt.Blob.betterEllipsoid([u11,u12,u13;u12,u22,u23;u13,u23,u33]).rescale(v),...
                Sol(:,params),...
                'OutputFormat','uniform');
            
            if hasUaniso
                params = {'occupancy','u11','u22','u33','u12','u13','u23'};
                fffun = @(ff,o,u11,u22,u33,u12,u13,u23) ff.rescale(o).addU([u11,u12,u13;u12,u22,u23;u13,u23,u33]);
            elseif hasBiso
                params = {'occupancy','tempFactor'};
                fffun = @(ff,o,b) ff.rescale(o).addB(b);
            else
                warning('this should not happen');
            end
            
            % look up occupancy etc from parent structure
            [~,ia] = ismember(SolFF(:,uid_keys),Atoms(:,uid_keys),'rows');
            
            SolFF.mdxFormFactor = rowfun(fffun,...
                [SolFF(:,'mdxFormFactor'),Atoms(ia,params)],...
                'OutputFormat','uniform');
            
        end
        
        
        function [nlist,nexcl] = asu_neighbor_search(x,y,z,maxpairdist,groupSele)
            
            xyz0 = [x(:),y(:),z(:)];
            nlist = rangesearch(xyz0,xyz0,maxpairdist);
            
            if nargin < 5
                groupSele = true(size(xyz,1),1); % only one selection, and they're all in it
            end
            
            % now, filter to remove self-references and deal with alt conformers
            nAtoms = size(xyz0,1);
            
            nexcl = cell(size(nlist));
            
            for j=1:nAtoms
                isIncl = true(size(nlist{j}));
                isIncl = isIncl & nlist{j} ~= j; % remove self references
                
                % make sure the atoms are present in the same group at least once
                isSameGroup = any(and(groupSele(j,:),groupSele(nlist{j},:)),2);
                isIncl = isIncl & isSameGroup';
                
                nexcl{j} = nlist{j}(~isIncl);
                nlist{j} = nlist{j}(isIncl);
                
            end
            
        end
        
        
        
        function [nlist,Aout] = symmetry_neighbor_search(x,y,z,B,Ops,maxpairdist)
            
            xyz0 = [x(:),y(:),z(:)];
            
            com = mean(xyz0,1);
            rmax = sqrt(max(sum((xyz0-com).^2,2)));
            
            % convert to fractional coordinates
            maxcomdist = 2*rmax + maxpairdist;
            
            T = neighbor_ops(com,B,Ops,maxcomdist);
            
            % now, symmetry expand and filter by distance using rangesearch
            
            nsymop = size(T,1);
            
            nlist = cell(nsymop,1);
            Aout = cell(nsymop,1);
            
            f0 = inv(B)*xyz0';
            
            for j=1:nsymop
                fj = Ops(T.isymop(j))*f0 + [T.xshift(j);T.yshift(j);T.zshift(j)];
                xyzj = (B*fj)';
                n = rangesearch(xyzj,xyz0,maxpairdist);
                if ~all(cellfun(@isempty,n))
                    [iatom,~,c] = unique(cell2mat(n')');
                    n = mat2cell(c,cellfun(@numel,n));
                    nlist{j} = n;
                    Aout{j} = [table(iatom),repmat(T(j,:),numel(iatom),1)];
                end
            end
            
            % filter
            isIncl = ~cellfun(@isempty,nlist);
            
            nlist = nlist(isIncl);
            Aout = Aout(isIncl);
            
            % flatten
            nn = cellfun(@(t) size(t,1),Aout); % number of neighbors
            nadd = cumsum([0;nn(1:(end-1))]); % index to add
            
            nsymop = numel(nlist);
            
            for j=1:nsymop
                nlist{j} = cellfun(@(n) n + nadd(j),nlist{j},'Uni',0);
            end
            
            % flatten Aout
            Aout = vertcat(Aout{:});
            
            % flatten nlist
            nlist = horzcat(nlist{:});
            tmp = cell(size(nlist,1),1);
            for j=1:size(nlist,1) % loop over atoms
                tmp{j} = vertcat(nlist{j,:});
            end
            nlist = tmp;
            
        end
        
    end
end



function [ur11,ur22,ur33,ur12,ur13,ur23] = transform_adp(R,u11,u22,u33,u12,u13,u23)

if numel(u11) > 1 % vectorize
    [ur11,ur22,ur33,ur12,ur13,ur23] = arrayfun(...
        @(u11,u22,u33,u12,u13,u23) transform_adp(R,u11,u22,u33,u12,u13,u23),...
        u11,u22,u33,u12,u13,u23);
else
    U = [u11,u12,u13; u12, u22, u23; u13, u23, u33];
    Ur = R*U*R';
    ur11 = Ur(1,1);
    ur22 = Ur(2,2);
    ur33 = Ur(3,3);
    ur12 = Ur(1,2);
    ur13 = Ur(1,3);
    ur23 = Ur(2,3);
end

end



function T = neighbor_ops(xcen,B,Ops,dmax)
% return a table of symmetry operators and shifted coordinates

[b1,b2,b3] = bounding_box(B,xcen,dmax);

f0 = inv(B)*xcen(:);
nOps = numel(Ops);
T = cell(nOps,1);

for j=1:nOps
    fj = Ops(j)*f0;
    sr1 = [ceil(b1(1)-fj(1)),floor(b1(2)-fj(1))];
    sr2 = [ceil(b2(1)-fj(2)),floor(b2(2)-fj(2))];
    sr3 = [ceil(b3(1)-fj(3)),floor(b3(2)-fj(3))];
    
    [isymop,xshift,yshift,zshift] = ndgrid(j,sr1(1):sr1(2),sr2(1):sr2(2),sr3(1):sr3(2));
    fs = fj + [xshift(:),yshift(:),zshift(:)]';
    isIncl = sum((B*(fs-f0)).^2,1) < dmax^2;
    
    % get rid of self op?
    if Ops(j)==symm.SymmetryOperator() % if Ops(j) is the identity operator
        isIncl(xshift==0 & yshift==0 & zshift==0) = false;
    end
    
    isymop = isymop(isIncl);
    xshift = xshift(isIncl);
    yshift = yshift(isIncl);
    zshift = zshift(isIncl);
    
    isymop = isymop(:);
    xshift = xshift(:);
    yshift = yshift(:);
    zshift = zshift(:);
    
    T{j} = table(isymop,xshift,yshift,zshift);
end

T = vertcat(T{:});

end


function [f1,f2,f3] = bounding_box(B,ori,rmax)

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

end