classdef PointGroup
    %POINTGROUP 
    %
    % note: this is an old piece of code that's used mainly by proc.
    % It will be replaced by symm.PointGroup
    
    properties (SetAccess=immutable)
        name
        LaueGroup
        generalPositionsXYZ
        generalPositionsMat
        generatorsXYZ
        generatorsMat
        multiplicity
        testASU
    end

    methods
        function obj = PointGroup(searchID)
            if nargin==0
                searchID = '1';
            end

            if ischar(searchID)
                tblidx = find(...
                    strcmpi(searchID,obj.pointGroupNames),1,'first');
                if tblidx % it's a point group?
                    ispg = true;
                else % it's a laue group?
                    ispg = false;
                    tblidx = find(...
                        strcmpi(searchID,obj.laueGroupNames),1,'first');
                end
            else
                tblidx = [];
            end
            assert(~isempty(tblidx),'did not find point or laue group');

            if ispg % it's a point group
                obj.name = obj.pointGroupNames{tblidx};
                laueGroupName = obj.laueGroupNames{tblidx};
                obj.LaueGroup = geom.lattice.PointGroup(laueGroupName);
                obj.generalPositionsXYZ = obj.allxyz{tblidx};
                obj.generalPositionsMat = obj.allmat{tblidx};
                obj.generatorsXYZ = obj.genxyz{tblidx};
                obj.generatorsMat = obj.genmat{tblidx};
                obj.multiplicity = length(obj.generalPositionsMat);
            else
                obj.name = obj.laueGroupNames{tblidx};
                obj.generalPositionsXYZ = obj.allxyzlaue{tblidx};
                obj.generalPositionsMat = obj.allmatlaue{tblidx};
                obj.generatorsXYZ = obj.genxyzlaue{tblidx};
                obj.generatorsMat = obj.genmatlaue{tblidx};
                obj.multiplicity = length(obj.generalPositionsMat);
                obj.testASU = obj.testReciprocalASU{tblidx};
            end
        end
    end

    properties (Constant=true, Hidden=true)
        pointGroupNames = {'1', '3', '312', '321', '2', '4', '422', '6', '622', '23', '432', '222'};

        laueGroupNames = {'-1','-3','-3 1 m','-3 m 1','2/m','4/m','4/m m m','6/m','6/m m m','m -3','m -3 m','m m m'};

testReciprocalASU = {@(h,k,l)(0<=l)&~(l==0&h<0)&~(h==0&l==0&k<0);
@(h,k,l)(0<=h&0<=k)&~(k==0&l<0)&~(h==0&l<=0);
@(h,k,l)(0<=h&h<=k)&~(h==0&l<0);
@(h,k,l)(0<=h&0<=k&0<=l)&~(l==0&h<k);
@(h,k,l)(0<=k&0<=l)&~(l==0&h<0);
@(h,k,l)(0<=h&0<=k&0<=l)&~(h==0&k>0);
@(h,k,l)(0<=h&h<=k&0<=l);
@(h,k,l)(0<=h&0<=k&0<=l)&~(h==0&k>0);
@(h,k,l)(0<=h&h<=k&0<=l);
@(h,k,l)(0<=h&h<=k&h<=l)&~(h==l&h<k);
@(h,k,l)(0<=h&h<=k&k<=l);
@(h,k,l)(0<=h&0<=k&0<=l)};

        genxyz = {{' x, y, z'};
            {' x, y, z', '-y,x-y, z'};
            {' x, y, z', '-y,x-y, z', '-y,-x,-z'};
            {' x, y, z', '-y,x-y, z', ' y, x,-z'};
            {' x, y, z', '-x,  y,-z'};
            {' x, y, z', '-x, -y, z', '-y, x, z'};
            {' x, y, z', '-x, -y, z', '-y, x, z', '-x, y,-z'};
            {' x, y, z', '-y,x-y, z', '-x,-y, z'};
            {' x, y, z', '-y,x-y, z', '-x,-y, z', ' y, x,-z'};
            {' x, y, z', '-x, -y, z', '-x, y,-z', ' z, x, y'};
            {' x, y, z', '-x, -y, z', '-x, y,-z', ' z, x, y', ' y, x,-z'};
            {' x, y, z', '-x, -y, z', '-x, y,-z'}};

        allxyz = {{'x, y, z'};
            {'x, y, z', '-y, x - y, z', 'y - x, -x, z'};
            {'x, y, z', '-y, x - y, z', '-y, -x, -z', 'y - x, -x, z', 'y - x, y, -z', 'x, x - y, -z'};
            {'x, y, z', '-y, x - y, z', 'y, x, -z', 'y - x, -x, z', 'x - y, -y, -z', '-x, y - x, -z'};
            {'x, y, z', '-x, y, -z'};
            {'x, y, z', '-x, -y, z', '-y, x, z', 'y, -x, z'};
            {'x, y, z', '-x, -y, z', '-y, x, z', '-x, y, -z', 'y, -x, z', 'x, -y, -z', 'y, x, -z', '-y, -x, -z'};
            {'x, y, z', '-y, x - y, z', '-x, -y, z', 'y - x, -x, z', 'y, y - x, z', 'x - y, x, z'};
            {'x, y, z', '-y, x - y, z', '-x, -y, z', 'y, x, -z', 'y - x, -x, z', 'y, y - x, z', 'x - y, -y, -z', '-y, -x, -z', '-x, y - x, -z', 'x - y, x, z', 'y - x, y, -z', 'x, x - y, -z'};
            {'x, y, z', '-x, -y, z', '-x, y, -z', 'z, x, y', 'x, -y, -z', 'z, -x, -y', '-z, -x, y', '-z, x, -y', 'y, z, x', '-y, z, -x', 'y, -z, -x', '-y, -z, x'};
            {'x, y, z', '-x, -y, z', '-x, y, -z', 'z, x, y', 'y, x, -z', 'x, -y, -z', 'z, -x, -y', '-y, -x, -z', '-z, -x, y', 'y, -x, z', '-z, x, -y', 'y, z, x', 'x, z, -y', '-y, x, z', '-z, y, x', '-y, z, -x', '-x, z, y', '-z, -y, -x', 'y, -z, -x', '-x, -z, -y', 'z, y, -x', '-y, -z, x', 'x, -z, y', 'z, -y, x'};
            {'x, y, z', '-x, -y, z', '-x, y, -z', 'x, -y, -z'}};

        genxyzlaue = {{' x, y, z', '-x, -y, -z'};
            {' x, y, z', '-y,x-y, z', '-x, -y, -z'};
            {' x, y, z', '-y,x-y, z', '-y,-x,-z', '-x, -y, -z'};
            {' x, y, z', '-y,x-y, z', ' y, x,-z', '-x, -y, -z'};
            {' x, y, z', '-x,  y,-z', '-x, -y, -z'};
            {' x, y, z', '-x, -y, z', '-y, x, z', '-x, -y, -z'};
            {' x, y, z', '-x, -y, z', '-y, x, z', '-x, y,-z', '-x, -y, -z'};
            {' x, y, z', '-y,x-y, z', '-x,-y, z', '-x, -y, -z'};
            {' x, y, z', '-y,x-y, z', '-x,-y, z', ' y, x,-z', '-x, -y, -z'};
            {' x, y, z', '-x, -y, z', '-x, y,-z', ' z, x, y', '-x, -y, -z'};
            {' x, y, z', '-x, -y, z', '-x, y,-z', ' z, x, y', ' y, x,-z', '-x, -y, -z'};
            {' x, y, z', '-x, -y, z', '-x, y,-z', '-x, -y, -z'}};

        allxyzlaue = {{'x, y, z', '-x, -y, -z'};
            {'x, y, z', '-y, x - y, z', '-x, -y, -z', 'y - x, -x, z', 'y, y - x, -z', 'x - y, x, -z'};
            {'x, y, z', '-y, x - y, z', '-y, -x, -z', '-x, -y, -z', 'y - x, -x, z', 'y - x, y, -z', 'y, y - x, -z', 'x, x - y, -z', 'y, x, z', 'x - y, x, -z', 'x - y, -y, z', '-x, y - x, z'};
            {'x, y, z', '-y, x - y, z', 'y, x, -z', '-x, -y, -z', 'y - x, -x, z', 'x - y, -y, -z', 'y, y - x, -z', '-x, y - x, -z', '-y, -x, z', 'x - y, x, -z', 'y - x, y, z', 'x, x - y, z'};
            {'x, y, z', '-x, y, -z', '-x, -y, -z', 'x, -y, z'};
            {'x, y, z', '-x, -y, z', '-y, x, z', '-x, -y, -z', 'y, -x, z', 'x, y, -z', 'y, -x, -z', '-y, x, -z'};
            {'x, y, z', '-x, -y, z', '-y, x, z', '-x, y, -z', '-x, -y, -z', 'y, -x, z', 'x, -y, -z', 'x, y, -z', 'y, x, -z', 'y, -x, -z', '-y, -x, -z', 'x, -y, z', '-y, x, -z', '-x, y, z', '-y, -x, z', 'y, x, z'};
            {'x, y, z', '-y, x - y, z', '-x, -y, z', '-x, -y, -z', 'y - x, -x, z', 'y, y - x, z', 'y, y - x, -z', 'x, y, -z', 'x - y, x, z', 'x - y, x, -z', '-y, x - y, -z', 'y - x, -x, -z'};
            {'x, y, z', '-y, x - y, z', '-x, -y, z', 'y, x, -z', '-x, -y, -z', 'y - x, -x, z', 'y, y - x, z', 'x - y, -y, -z', 'y, y - x, -z', '-y, -x, -z', 'x, y, -z', '-x, y - x, -z', '-y, -x, z', 'x - y, x, z', 'x - y, x, -z', 'y - x, y, -z', '-y, x - y, -z', 'y - x, y, z', 'x, x - y, -z', 'y, x, z', 'x, x - y, z', 'y - x, -x, -z', 'x - y, -y, z', '-x, y - x, z'};
            {'x, y, z', '-x, -y, z', '-x, y, -z', 'z, x, y', '-x, -y, -z', 'x, -y, -z', 'z, -x, -y', 'x, y, -z', '-z, -x, y', 'x, -y, z', '-z, x, -y', 'y, z, x', '-z, -x, -y', '-x, y, z', '-y, z, -x', '-z, x, y', 'y, -z, -x', 'z, x, -y', '-y, -z, x', 'z, -x, y', '-y, -z, -x', 'y, -z, x', '-y, z, x', 'y, z, -x'};
            {'x, y, z', '-x, -y, z', '-x, y, -z', 'z, x, y', 'y, x, -z', '-x, -y, -z', 'x, -y, -z', 'z, -x, -y', '-y, -x, -z', 'x, y, -z', '-z, -x, y', 'y, -x, z', 'x, -y, z', '-z, x, -y', 'y, z, x', 'x, z, -y', '-z, -x, -y', '-y, x, z', '-z, y, x', '-y, -x, z', '-x, y, z', '-y, z, -x', '-x, z, y', '-z, x, y', '-z, -y, -x', 'y, x, z', 'y, -z, -x', '-x, -z, -y', 'z, x, -y', 'z, y, -x', '-y, x, -z', '-y, -z, x', 'x, -z, y', 'z, -x, y', '-y, -z, -x', '-x, -z, y', 'z, -y, x', 'y, -x, -z', 'z, -y, -x', 'y, -z, x', 'x, -z, -y', 'z, y, x', '-y, z, x', 'x, z, y', '-z, -y, x', 'y, z, -x', '-x, z, -y', '-z, y, -x'};
            {'x, y, z', '-x, -y, z', '-x, y, -z', '-x, -y, -z', 'x, -y, -z', 'x, y, -z', 'x, -y, z', '-x, y, z'}};

        genmat = {{[1 0 0;0 1 0;0 0 1]};
            {[1 0 0;0 1 0;0 0 1], [0 -1 0;1 -1 0;0 0 1]};
            {[1 0 0;0 1 0;0 0 1], [0 -1 0;1 -1 0;0 0 1], [0 -1 0;-1 0 0;0 0 -1]};
            {[1 0 0;0 1 0;0 0 1], [0 -1 0;1 -1 0;0 0 1], [0 1 0;1 0 0;0 0 -1]};
            {[1 0 0;0 1 0;0 0 1], [-1 0 0;0 1 0;0 0 -1]};
            {[1 0 0;0 1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [0 -1 0;1 0 0;0 0 1]};
            {[1 0 0;0 1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [0 -1 0;1 0 0;0 0 1], [-1 0 0;0 1 0;0 0 -1]};
            {[1 0 0;0 1 0;0 0 1], [0 -1 0;1 -1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1]};
            {[1 0 0;0 1 0;0 0 1], [0 -1 0;1 -1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [0 1 0;1 0 0;0 0 -1]};
            {[1 0 0;0 1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [-1 0 0;0 1 0;0 0 -1], [0 0 1;1 0 0;0 1 0]};
            {[1 0 0;0 1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [-1 0 0;0 1 0;0 0 -1], [0 0 1;1 0 0;0 1 0], [0 1 0;1 0 0;0 0 -1]};
            {[1 0 0;0 1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [-1 0 0;0 1 0;0 0 -1]}};

        allmat = {{[1 0 0;0 1 0;0 0 1]};
            {[1 0 0;0 1 0;0 0 1], [0 -1 0;1 -1 0;0 0 1], [-1 1 0;-1 0 0;0 0 1]};
            {[1 0 0;0 1 0;0 0 1], [0 -1 0;1 -1 0;0 0 1], [0 -1 0;-1 0 0;0 0 -1], [-1 1 0;-1 0 0;0 0 1], [-1 1 0;0 1 0;0 0 -1], [1 0 0;1 -1 0;0 0 -1]};
            {[1 0 0;0 1 0;0 0 1], [0 -1 0;1 -1 0;0 0 1], [0 1 0;1 0 0;0 0 -1], [-1 1 0;-1 0 0;0 0 1], [1 -1 0;0 -1 0;0 0 -1], [-1 0 0;-1 1 0;0 0 -1]};
            {[1 0 0;0 1 0;0 0 1], [-1 0 0;0 1 0;0 0 -1]};
            {[1 0 0;0 1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [0 -1 0;1 0 0;0 0 1], [0 1 0;-1 0 0;0 0 1]};
            {[1 0 0;0 1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [0 -1 0;1 0 0;0 0 1], [-1 0 0;0 1 0;0 0 -1], [0 1 0;-1 0 0;0 0 1], [1 0 0;0 -1 0;0 0 -1], [0 1 0;1 0 0;0 0 -1], [0 -1 0;-1 0 0;0 0 -1]};
            {[1 0 0;0 1 0;0 0 1], [0 -1 0;1 -1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [-1 1 0;-1 0 0;0 0 1], [0 1 0;-1 1 0;0 0 1], [1 -1 0;1 0 0;0 0 1]};
            {[1 0 0;0 1 0;0 0 1], [0 -1 0;1 -1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [0 1 0;1 0 0;0 0 -1], [-1 1 0;-1 0 0;0 0 1], [0 1 0;-1 1 0;0 0 1], [1 -1 0;0 -1 0;0 0 -1], [0 -1 0;-1 0 0;0 0 -1], [-1 0 0;-1 1 0;0 0 -1], [1 -1 0;1 0 0;0 0 1], [-1 1 0;0 1 0;0 0 -1], [1 0 0;1 -1 0;0 0 -1]};
            {[1 0 0;0 1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [-1 0 0;0 1 0;0 0 -1], [0 0 1;1 0 0;0 1 0], [1 0 0;0 -1 0;0 0 -1], [0 0 1;-1 0 0;0 -1 0], [0 0 -1;-1 0 0;0 1 0], [0 0 -1;1 0 0;0 -1 0], [0 1 0;0 0 1;1 0 0], [0 -1 0;0 0 1;-1 0 0], [0 1 0;0 0 -1;-1 0 0], [0 -1 0;0 0 -1;1 0 0]};
            {[1 0 0;0 1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [-1 0 0;0 1 0;0 0 -1], [0 0 1;1 0 0;0 1 0], [0 1 0;1 0 0;0 0 -1], [1 0 0;0 -1 0;0 0 -1], [0 0 1;-1 0 0;0 -1 0], [0 -1 0;-1 0 0;0 0 -1], [0 0 -1;-1 0 0;0 1 0], [0 1 0;-1 0 0;0 0 1], [0 0 -1;1 0 0;0 -1 0], [0 1 0;0 0 1;1 0 0], [1 0 0;0 0 1;0 -1 0], [0 -1 0;1 0 0;0 0 1], [0 0 -1;0 1 0;1 0 0], [0 -1 0;0 0 1;-1 0 0], [-1 0 0;0 0 1;0 1 0], [0 0 -1;0 -1 0;-1 0 0], [0 1 0;0 0 -1;-1 0 0], [-1 0 0;0 0 -1;0 -1 0], [0 0 1;0 1 0;-1 0 0], [0 -1 0;0 0 -1;1 0 0], [1 0 0;0 0 -1;0 1 0], [0 0 1;0 -1 0;1 0 0]};
            {[1 0 0;0 1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [-1 0 0;0 1 0;0 0 -1], [1 0 0;0 -1 0;0 0 -1]}};

        genmatlaue = {{[1 0 0;0 1 0;0 0 1], [-1 0 0;0 -1 0;0 0 -1]};
            {[1 0 0;0 1 0;0 0 1], [0 -1 0;1 -1 0;0 0 1], [-1 0 0;0 -1 0;0 0 -1]};
            {[1 0 0;0 1 0;0 0 1], [0 -1 0;1 -1 0;0 0 1], [0 -1 0;-1 0 0;0 0 -1], [-1 0 0;0 -1 0;0 0 -1]};
            {[1 0 0;0 1 0;0 0 1], [0 -1 0;1 -1 0;0 0 1], [0 1 0;1 0 0;0 0 -1], [-1 0 0;0 -1 0;0 0 -1]};
            {[1 0 0;0 1 0;0 0 1], [-1 0 0;0 1 0;0 0 -1], [-1 0 0;0 -1 0;0 0 -1]};
            {[1 0 0;0 1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [0 -1 0;1 0 0;0 0 1], [-1 0 0;0 -1 0;0 0 -1]};
            {[1 0 0;0 1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [0 -1 0;1 0 0;0 0 1], [-1 0 0;0 1 0;0 0 -1], [-1 0 0;0 -1 0;0 0 -1]};
            {[1 0 0;0 1 0;0 0 1], [0 -1 0;1 -1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [-1 0 0;0 -1 0;0 0 -1]};
            {[1 0 0;0 1 0;0 0 1], [0 -1 0;1 -1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [0 1 0;1 0 0;0 0 -1], [-1 0 0;0 -1 0;0 0 -1]};
            {[1 0 0;0 1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [-1 0 0;0 1 0;0 0 -1], [0 0 1;1 0 0;0 1 0], [-1 0 0;0 -1 0;0 0 -1]};
            {[1 0 0;0 1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [-1 0 0;0 1 0;0 0 -1], [0 0 1;1 0 0;0 1 0], [0 1 0;1 0 0;0 0 -1], [-1 0 0;0 -1 0;0 0 -1]};
            {[1 0 0;0 1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [-1 0 0;0 1 0;0 0 -1], [-1 0 0;0 -1 0;0 0 -1]}};

        allmatlaue = {{[1 0 0;0 1 0;0 0 1], [-1 0 0;0 -1 0;0 0 -1]};
            {[1 0 0;0 1 0;0 0 1], [0 -1 0;1 -1 0;0 0 1], [-1 0 0;0 -1 0;0 0 -1], [-1 1 0;-1 0 0;0 0 1], [0 1 0;-1 1 0;0 0 -1], [1 -1 0;1 0 0;0 0 -1]};
            {[1 0 0;0 1 0;0 0 1], [0 -1 0;1 -1 0;0 0 1], [0 -1 0;-1 0 0;0 0 -1], [-1 0 0;0 -1 0;0 0 -1], [-1 1 0;-1 0 0;0 0 1], [-1 1 0;0 1 0;0 0 -1], [0 1 0;-1 1 0;0 0 -1], [1 0 0;1 -1 0;0 0 -1], [0 1 0;1 0 0;0 0 1], [1 -1 0;1 0 0;0 0 -1], [1 -1 0;0 -1 0;0 0 1], [-1 0 0;-1 1 0;0 0 1]};
            {[1 0 0;0 1 0;0 0 1], [0 -1 0;1 -1 0;0 0 1], [0 1 0;1 0 0;0 0 -1], [-1 0 0;0 -1 0;0 0 -1], [-1 1 0;-1 0 0;0 0 1], [1 -1 0;0 -1 0;0 0 -1], [0 1 0;-1 1 0;0 0 -1], [-1 0 0;-1 1 0;0 0 -1], [0 -1 0;-1 0 0;0 0 1], [1 -1 0;1 0 0;0 0 -1], [-1 1 0;0 1 0;0 0 1], [1 0 0;1 -1 0;0 0 1]};
            {[1 0 0;0 1 0;0 0 1], [-1 0 0;0 1 0;0 0 -1], [-1 0 0;0 -1 0;0 0 -1], [1 0 0;0 -1 0;0 0 1]};
            {[1 0 0;0 1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [0 -1 0;1 0 0;0 0 1], [-1 0 0;0 -1 0;0 0 -1], [0 1 0;-1 0 0;0 0 1], [1 0 0;0 1 0;0 0 -1], [0 1 0;-1 0 0;0 0 -1], [0 -1 0;1 0 0;0 0 -1]};
            {[1 0 0;0 1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [0 -1 0;1 0 0;0 0 1], [-1 0 0;0 1 0;0 0 -1], [-1 0 0;0 -1 0;0 0 -1], [0 1 0;-1 0 0;0 0 1], [1 0 0;0 -1 0;0 0 -1], [1 0 0;0 1 0;0 0 -1], [0 1 0;1 0 0;0 0 -1], [0 1 0;-1 0 0;0 0 -1], [0 -1 0;-1 0 0;0 0 -1], [1 0 0;0 -1 0;0 0 1], [0 -1 0;1 0 0;0 0 -1], [-1 0 0;0 1 0;0 0 1], [0 -1 0;-1 0 0;0 0 1], [0 1 0;1 0 0;0 0 1]};
            {[1 0 0;0 1 0;0 0 1], [0 -1 0;1 -1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [-1 0 0;0 -1 0;0 0 -1], [-1 1 0;-1 0 0;0 0 1], [0 1 0;-1 1 0;0 0 1], [0 1 0;-1 1 0;0 0 -1], [1 0 0;0 1 0;0 0 -1], [1 -1 0;1 0 0;0 0 1], [1 -1 0;1 0 0;0 0 -1], [0 -1 0;1 -1 0;0 0 -1], [-1 1 0;-1 0 0;0 0 -1]};
            {[1 0 0;0 1 0;0 0 1], [0 -1 0;1 -1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [0 1 0;1 0 0;0 0 -1], [-1 0 0;0 -1 0;0 0 -1], [-1 1 0;-1 0 0;0 0 1], [0 1 0;-1 1 0;0 0 1], [1 -1 0;0 -1 0;0 0 -1], [0 1 0;-1 1 0;0 0 -1], [0 -1 0;-1 0 0;0 0 -1], [1 0 0;0 1 0;0 0 -1], [-1 0 0;-1 1 0;0 0 -1], [0 -1 0;-1 0 0;0 0 1], [1 -1 0;1 0 0;0 0 1], [1 -1 0;1 0 0;0 0 -1], [-1 1 0;0 1 0;0 0 -1], [0 -1 0;1 -1 0;0 0 -1], [-1 1 0;0 1 0;0 0 1], [1 0 0;1 -1 0;0 0 -1], [0 1 0;1 0 0;0 0 1], [1 0 0;1 -1 0;0 0 1], [-1 1 0;-1 0 0;0 0 -1], [1 -1 0;0 -1 0;0 0 1], [-1 0 0;-1 1 0;0 0 1]};
            {[1 0 0;0 1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [-1 0 0;0 1 0;0 0 -1], [0 0 1;1 0 0;0 1 0], [-1 0 0;0 -1 0;0 0 -1], [1 0 0;0 -1 0;0 0 -1], [0 0 1;-1 0 0;0 -1 0], [1 0 0;0 1 0;0 0 -1], [0 0 -1;-1 0 0;0 1 0], [1 0 0;0 -1 0;0 0 1], [0 0 -1;1 0 0;0 -1 0], [0 1 0;0 0 1;1 0 0], [0 0 -1;-1 0 0;0 -1 0], [-1 0 0;0 1 0;0 0 1], [0 -1 0;0 0 1;-1 0 0], [0 0 -1;1 0 0;0 1 0], [0 1 0;0 0 -1;-1 0 0], [0 0 1;1 0 0;0 -1 0], [0 -1 0;0 0 -1;1 0 0], [0 0 1;-1 0 0;0 1 0], [0 -1 0;0 0 -1;-1 0 0], [0 1 0;0 0 -1;1 0 0], [0 -1 0;0 0 1;1 0 0], [0 1 0;0 0 1;-1 0 0]};
            {[1 0 0;0 1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [-1 0 0;0 1 0;0 0 -1], [0 0 1;1 0 0;0 1 0], [0 1 0;1 0 0;0 0 -1], [-1 0 0;0 -1 0;0 0 -1], [1 0 0;0 -1 0;0 0 -1], [0 0 1;-1 0 0;0 -1 0], [0 -1 0;-1 0 0;0 0 -1], [1 0 0;0 1 0;0 0 -1], [0 0 -1;-1 0 0;0 1 0], [0 1 0;-1 0 0;0 0 1], [1 0 0;0 -1 0;0 0 1], [0 0 -1;1 0 0;0 -1 0], [0 1 0;0 0 1;1 0 0], [1 0 0;0 0 1;0 -1 0], [0 0 -1;-1 0 0;0 -1 0], [0 -1 0;1 0 0;0 0 1], [0 0 -1;0 1 0;1 0 0], [0 -1 0;-1 0 0;0 0 1], [-1 0 0;0 1 0;0 0 1], [0 -1 0;0 0 1;-1 0 0], [-1 0 0;0 0 1;0 1 0], [0 0 -1;1 0 0;0 1 0], [0 0 -1;0 -1 0;-1 0 0], [0 1 0;1 0 0;0 0 1], [0 1 0;0 0 -1;-1 0 0], [-1 0 0;0 0 -1;0 -1 0], [0 0 1;1 0 0;0 -1 0], [0 0 1;0 1 0;-1 0 0], [0 -1 0;1 0 0;0 0 -1], [0 -1 0;0 0 -1;1 0 0], [1 0 0;0 0 -1;0 1 0], [0 0 1;-1 0 0;0 1 0], [0 -1 0;0 0 -1;-1 0 0], [-1 0 0;0 0 -1;0 1 0], [0 0 1;0 -1 0;1 0 0], [0 1 0;-1 0 0;0 0 -1], [0 0 1;0 -1 0;-1 0 0], [0 1 0;0 0 -1;1 0 0], [1 0 0;0 0 -1;0 -1 0], [0 0 1;0 1 0;1 0 0], [0 -1 0;0 0 1;1 0 0], [1 0 0;0 0 1;0 1 0], [0 0 -1;0 -1 0;1 0 0], [0 1 0;0 0 1;-1 0 0], [-1 0 0;0 0 1;0 -1 0], [0 0 -1;0 1 0;-1 0 0]};
            {[1 0 0;0 1 0;0 0 1], [-1 0 0;0 -1 0;0 0 1], [-1 0 0;0 1 0;0 0 -1], [-1 0 0;0 -1 0;0 0 -1], [1 0 0;0 -1 0;0 0 -1], [1 0 0;0 1 0;0 0 -1], [1 0 0;0 -1 0;0 0 1], [-1 0 0;0 1 0;0 0 1]}};

        % this table is used to generate all of the variables above (see
        % static hidden functions, below)
        masterTable =...
            {'1',       {' x, y, z'},@(h,k,l) (0<=l) & ~(l==0 & h<0) & ~(h==0 & l==0 & k<0) ;... % 1
            '3',        {' x, y, z','-y,x-y, z'},@(h,k,l) (0<=h & 0<=k) & ~(k==0 & l<0) & ~(h==0 & l<=0);... % 3
            '312',      {' x, y, z','-y,x-y, z','-y,-x,-z'},@(h,k,l) (0<=h & h<=k) & ~(h==0 & l<0);... % 312, note: sometimes called just 32
            '321',      {' x, y, z','-y,x-y, z',' y, x,-z'},@(h,k,l) (0<=h & 0<=k & 0<=l) & ~(l==0 & h<k);... % 321, note: sometimes called just 32
            '2',        {' x, y, z','-x,  y,-z'},@(h,k,l) (0<=k & 0<=l) & ~(l==0 & h<0);... %2
            '4',        {' x, y, z','-x, -y, z','-y, x, z'},@(h,k,l) (0<=h & 0<=k & 0<=l) & ~(h==0 & k>0);... % 4
            '422',      {' x, y, z','-x, -y, z','-y, x, z','-x, y,-z'},@(h,k,l) (0<=h & h<=k & 0<=l);... % 422
            '6',        {' x, y, z','-y,x-y, z','-x,-y, z'},@(h,k,l) (0<=h & 0<=k & 0<=l) & ~(h==0 & k>0);... % 6
            '622',      {' x, y, z','-y,x-y, z','-x,-y, z',' y, x,-z'},@(h,k,l) (0<=h & h<=k & 0<=l);... % 622
            '23',       {' x, y, z','-x, -y, z','-x, y,-z',' z, x, y'},@(h,k,l) (0<=h & h<=k & h<=l) & ~(h==l & h<k);... % 23
            '432',      {' x, y, z','-x, -y, z','-x, y,-z',' z, x, y',' y, x,-z'},@(h,k,l) (0<=h & h<=k & k<=l);... % 432
            '222',      {' x, y, z','-x, -y, z','-x, y,-z'},@(h,k,l) (0<=h & 0<=k & 0<=l)... % 222
            };
    end

    methods(Static,Hidden = true)

        % methods for re-generating operator database from master table
        function printPointGroupNames()
            t = geom.lattice.PointGroup.masterTable;
            fprintf(1,'pointGroupNames = {');
            fprintf(1,'''%s'', ',t{:,1});
            fprintf(1,'\b\b};\n');
        end

        function printTestReciprocalASU()
            t = geom.lattice.PointGroup.masterTable;
            fprintf(1,'testReciprocalASU = {');
            c = cellfun(@(v) char(v),t(:,3),'uniformOutput',0);
            fprintf(1,'%s;\n',c{:});
            fprintf(1,'\b\b};\n');
        end

        function printgenxyz()
            t = geom.lattice.PointGroup.masterTable;
            fprintf(1,'genxyz = {');
            for j=1:size(t,1)
                fprintf(1,'{');
                fprintf(1,'''%s'', ',t{j,2}{:});
                fprintf(1,'\b\b};\n');
            end
            fprintf(1,'\b\b};\n');
        end

        function printgenxyzlaue()
            t = geom.lattice.PointGroup.masterTable;
            fprintf(1,'genxyzlaue = {');
            for j=1:size(t,1)
                fprintf(1,'{');
                fprintf(1,'''%s'', ',t{j,2}{:});
                fprintf(1,'''%s'', ','-x, -y, -z');
                fprintf(1,'\b\b};\n');
            end
            fprintf(1,'\b\b};\n');
        end

        function printgenmat()
            t = geom.lattice.PointGroup.masterTable;
            fprintf(1,'genmat = {');
            for j=1:size(t,1)
                op = geom.lattice.PointGroup.xyz2op(t{j,2});
                thisopstr = cellfun(@(m) mat2str(m),op,'uniformOutput',0);
                fprintf(1,'{');
                fprintf(1,'%s, ',thisopstr{:});
                fprintf(1,'\b\b};\n');
            end
            fprintf(1,'\b\b};\n');
        end

        function printgenmatlaue()
            t = geom.lattice.PointGroup.masterTable;
            fprintf(1,'genmatlaue = {');
            for j=1:size(t,1)
                op = geom.lattice.PointGroup.xyz2op(t{j,2});
                op = [op,-eye(3)]; % add inversion operator
                thisopstr = cellfun(@(m) mat2str(m),op,'uniformOutput',0);
                fprintf(1,'{');
                fprintf(1,'%s, ',thisopstr{:});
                fprintf(1,'\b\b};\n');
            end
            fprintf(1,'\b\b};\n');
        end

        function printallmat()
            t = geom.lattice.PointGroup.masterTable;
            fprintf(1,'allmat = {');
            for j=1:size(t,1)
                op = geom.lattice.PointGroup.xyz2op(t{j,2});
                op = geom.lattice.PointGroup.op2general(op);
                thisopstr = cellfun(@(m) mat2str(m),op,'uniformOutput',0);
                fprintf(1,'{');
                fprintf(1,'%s, ',thisopstr{:});
                fprintf(1,'\b\b};\n');
            end
            fprintf(1,'\b\b};\n');
        end

        function printallmatlaue() % add the inversion operator
            t = geom.lattice.PointGroup.masterTable;
            fprintf(1,'allmatlaue = {');
            for j=1:size(t,1)
                op = geom.lattice.PointGroup.xyz2op(t{j,2});
                op = [op,-eye(3)]; % add inversion operator
                op = geom.lattice.PointGroup.op2general(op);
                thisopstr = cellfun(@(m) mat2str(m),op,'uniformOutput',0);
                fprintf(1,'{');
                fprintf(1,'%s, ',thisopstr{:});
                fprintf(1,'\b\b};\n');
            end
            fprintf(1,'\b\b};\n');
        end

        function printallxyzlaue()
            t = geom.lattice.PointGroup.masterTable;
            fprintf(1,'allxyzlaue = {');
            for j=1:size(t,1)
                op = geom.lattice.PointGroup.xyz2op(t{j,2});
                op = [op,-eye(3)]; % add inversion operator
                op = geom.lattice.PointGroup.op2general(op);
                xyz = geom.lattice.PointGroup.op2xyz(op);
                fprintf(1,'{');
                fprintf(1,'''%s'', ',xyz{:});
                fprintf(1,'\b\b};\n');
            end
            fprintf(1,'\b\b};\n');
        end

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
