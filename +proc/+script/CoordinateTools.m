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
        
        function [] = neighbor_search(A,Basis,SpaceGroup,maxDist)
            % A = Atoms() table
            
            a = 15; % box size
            [A.xbin,A.ybin,A.zbin,xcen,ycen,zcen] = boxem(A.x,A.y,A.z,a);
            [b,~,ind] = unique(A(:,{'xbin','ybin','zbin'}),'rows');
            boxinfo = rowfun(@(n1,n2,n3) deal(xcen(n1),ycen(n2),zcen(n3)),b,...
                'OutputVariableNames',{'xcen','ycen','zcen'});
            ix = accumarray(ind,1:numel(ind),[],@(v) {v});
            
            Ops = SpaceGroup.generalPositions;
            B = Basis.orthogonalizationMatrix;
            xyz = table2array(A(:,{'x','y','z'}));
            
            T = table();
            
            for j=1:numel(ix)
                Aj = A(ix{j},:);
                rcen = table2array(boxinfo(j,{'xcen','ycen','zcen'}));
                
                [x,y,z,xshift,yshift,zshift,isymop,irow] = symmetry_neighbor(...
                    rcen(:),xyz',Ops,B);
                
                isIncl = abs(x-rcen(1)) < (a/2 + dmax) & abs(y-rcen(2)) < (a/2 + dmax) & abs(z-rcen(3)) < (a/2 + dmax); % box
                x = x(isIncl);
                y = y(isIncl);
                z = z(isIncl);
                xshift = xshift(isIncl);
                yshift = yshift(isIncl);
                zshift = zshift(isIncl);
                isymop = isymop(isIncl);
                irow = irow(isIncl);
                
                Tj = table(irow,isymop,xshift,yshift,zshift,x,y,z);
                
                T = [T;Tj];
                
            end
            
        end
        
        
        
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
                fprintf(1,'transforming u11, u22, ...\n');
                
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




function [d,xshift,yshift,zshift,isymop,ir1] = symmdist(r0,r1,Ops,B)

% convert to fractional coordinates
M = B'*B;
f0 = inv(B)*r0;
f1 = inv(B)*r1;

nOps = numel(Ops);
nPoints = size(r1,2);

isymop = repmat((1:nOps)',1,nPoints); % for reference
ir1 = repmat((1:nPoints),nOps,1); % for reference

d = zeros(nOps,nPoints);

xshift = zeros(nOps,nPoints);
yshift = zeros(nOps,nPoints);
zshift = zeros(nOps,nPoints);

for j=1:numel(Ops)
    f1j = Ops(j)*f1;
    df01 = f1j - f0;
    shift = - floor(df01 + 1/2);
    df01 = df01 + shift;
    d(j,:) = sqrt(dot(df01,M*df01));
    xshift(j,:) = shift(1,:);
    yshift(j,:) = shift(2,:);
    zshift(j,:) = shift(3,:);
end

end


function [xbin,ybin,zbin,xcen,ycen,zcen] = boxem(x,y,z,a)

[xmin,xmax] = bounds(x);
[ymin,ymax] = bounds(y);
[zmin,zmax] = bounds(z);

latt.Basis(a,a,a,90,90,90);
N = [ceil((xmax-xmin)/a),ceil((ymax-ymin)/a),ceil((zmax-zmin)/a)];
xEdges = (xmin+xmax)/2 + a*((0:N(1)) - N(1)/2);
yEdges = (ymin+ymax)/2 + a*((0:N(2)) - N(2)/2);
zEdges = (zmin+zmax)/2 + a*((0:N(3)) - N(3)/2);

xbin = discretize(x,xEdges);
ybin = discretize(y,yEdges);
zbin = discretize(z,zEdges);

xcen = xEdges(1:(end-1))*.5 + xEdges(2:end)*.5;
ycen = yEdges(1:(end-1))*.5 + yEdges(2:end)*.5;
zcen = zEdges(1:(end-1))*.5 + zEdges(2:end)*.5;

end