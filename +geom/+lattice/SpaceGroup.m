classdef SpaceGroup
    % SpaceGroup
    %
    % note: this is an old piece of code that's used mainly by proc.
    % It has been replaced by symm.SpaceGroup
    %
    % For now, maintain the legacy interface but back it using symm.SpaceGroupInfo
    properties (SetAccess=immutable)
        spaceGroupNumber
        inputName
        fullName
        crystalSystem
        PointGroup
        bravaisType
        checkProperties % function handle replacing original member function
    end

    methods
        function obj = SpaceGroup(searchID)
            
            sginfo = symm.SpaceGroupInfo(searchID);
            obj.spaceGroupNumber = sginfo.number;
            obj.inputName = sginfo.name;
            obj.fullName = sginfo.formattedName;
            obj.crystalSystem = sginfo.crystalSystem;
            obj.bravaisType = sginfo.centeringType;
            obj.PointGroup = geom.lattice.PointGroup(sginfo.pointGroupNumber);
            obj.checkProperties = sginfo.testLattice;
            
        end 
    end
end

