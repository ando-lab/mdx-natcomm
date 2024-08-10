classdef PointGroup
    %POINTGROUP 
    %
    % note: this is an old piece of code that's used mainly by proc.
    % It will be replaced by symm.PointGroup
    %
    % For now, keep the legacy interface but back it using
    % symm.SpaceGroupInfo
    %
    % Dropped support for XYZ format of operators -- is this used anywhere?
    
    properties (SetAccess=immutable)
        name
        LaueGroup
        %generalPositionsXYZ 
        generalPositionsMat
        %generatorsXYZ
        generatorsMat
        multiplicity
        testASU
    end

    methods
        function obj = PointGroup(searchID)
            sginfo = symm.SpaceGroupInfo(searchID);

            if ~sginfo.isLaueGroup && sginfo.isPointGroup % it's a point group, add Laue group assignment
                obj.LaueGroup = geom.lattice.PointGroup(sginfo.laueGroupNumber);
            elseif sginfo.isLaueGroup
                obj.testASU = symm.LaueGroup(searchID).testASU;
            else
                error('did not find point group or Laue group')
            end
            
            obj.name = sginfo.name;
            obj.generalPositionsMat = cellfun(@(v) v(1:3,1:3), ...
                sginfo.generalPositions,'Uni',0);
            obj.generatorsMat = cellfun(@(v) v(1:3,1:3), ...
                sginfo.generators,'Uni',0);
            obj.multiplicity = length(obj.generalPositionsMat);
        end
    end

    methods(Static,Hidden = true)

        function opList = xyz2op(xyzList)

            if iscell(xyzList)
                opList = cell(size(xyzList));
                for j=1:numel(xyzList)
                    opList{j} = xyzconv(xyzList{j});
                end
            else
                opList = xyzconv(xyzList);
            end

        end
        function xyzList = op2xyz(opList)
            if iscell(opList)
                xyzList = cell(size(opList));
                for j=1:numel(opList)
                    xyzList{j} = matconv(opList{j});
                end
            else
                xyzList = matconv(opList);
            end
        end
        function opList = op2general(opsIn)
            opList = opsIn;
            nOps = 0; % so that the while loop executes at least once
            nNewOps = length(opList);
            while nNewOps > nOps
                nOps = nNewOps;
                opList = binaryProducts(opList);
                opList = uniqueOperators(opList);
                nNewOps = length(opList);
            end
        end

    end

end

function op = xyzconv(xyz)
syms x y z
c = eval(['[',xyz,']']);

op = zeros(3,3);

op(:,1) = subs(c,[x,y,z],[1,0,0]);
op(:,2) = subs(c,[x,y,z],[0,1,0]);
op(:,3) = subs(c,[x,y,z],[0,0,1]);
end

function xyz = matconv(op)
syms x y z
fullstr = char(transpose(op*[x;y;z]));
matchstr = regexp(fullstr,'matrix\(\[\[(.*)\]\]\)','tokens');
xyz = matchstr{1}{1};
end

function opMat = binaryProducts(opList)
nOps = length(opList);
opMat = cell(nOps,nOps);

for i = 1:nOps
    for j=1:nOps
        opMat{i,j} = opList{i}*opList{j};
    end
end
end
function opList = uniqueOperators(opList)
opsAsRows = reshape(cat(3,opList{:}),9,[])';
[~,ia,~] = unique(opsAsRows,'rows','stable');
opList = opList(ia);
end
