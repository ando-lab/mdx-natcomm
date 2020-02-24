classdef SpaceGroup
    % SpaceGroup
    %
    % note: this is an old piece of code that's used mainly by proc.
    % It will be replaced by symm.SpaceGroup
    properties (SetAccess=immutable)
        spaceGroupNumber
        inputName
        fullName
        crystalSystem
        PointGroup
        bravaisType
    end
    
    properties (Constant=true, Hidden=true)
spaceGroupNumbers = [1, 3, 4, 5, 16, 17, 18, 19, 20, 21, 22, 23, 24, 75, 76, 77, 78, 79, 80, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 143, 144, 145, 146, 149, 150, 151, 152, 153, 154, 155, 168, 169, 170, 171, 172, 173, 177, 178, 179, 180, 181, 182, 195, 196, 197, 198, 199, 207, 208, 209, 210, 211, 212, 213, 214];
fullNames = {'P 1', 'P 2', 'P 2_1', 'C 2', 'P 2 2 2', 'P 2 2 2_1', 'P 2_1 2_1 2', 'P 2_1 2_1 2_1', 'C 2 2 2_1', 'C 2 2 2', 'F 2 2 2', 'I 2 2 2', 'I 2_1 2_1 2_1', 'P 4', 'P 4_1', 'P 4_2', 'P 4_3', 'I 4', 'I 4_1', 'P 4 2 2', 'P 4 2_1 2', 'P 4_1 2 2', 'P 4_1 2_1 2', 'P 4_2 2 2', 'P 4_2 2_1 2', 'P 4_3 2 2', 'P 4_3 2_1 2', 'I 4 2 2', 'I 4_1 2 2', 'P 3', 'P 3_1', 'P 3 2', 'R 3', 'P 3 1 2', 'P 3 2 1', 'P 3_1 1 2', 'P 3_1 2 1', 'P 3_2 1 2', 'P 3_2 2 1', 'R 3 2', 'P 6', 'P 6_1', 'P 6_5', 'P 6_2', 'P 6_4', 'P 6_3', 'P 6 2 2', 'P 6_1 2 2', 'P 6_5 2 2', 'P 6_2 2 2', 'P 6_4 2 2', 'P 6_3 2 2', 'P 2 3', 'F 2 3', 'I 2 3', 'P 2_1 3', 'I 2_1 3', 'P 4 3 2', 'P 4_2 3 2', 'F 4 3 2', 'F 4_1 3 2', 'I 4 3 2', 'P 4_3 3 2', 'P 4_1 3 2', 'I 4_1 3 2'};
inputNames = {'P1', 'P2', 'P21', 'C2', 'P222', 'P2221', 'P21212', 'P212121', 'C2221', 'C222', 'F222', 'I222', 'I212121', 'P4', 'P41', 'P42', 'P43', 'I4', 'I41', 'P422', 'P4212', 'P4122', 'P41212', 'P4222', 'P42212', 'P4322', 'P43212', 'I422', 'I4122', 'P3', 'P31', 'P32', 'R3', 'P312', 'P321', 'P3112', 'P3121', 'P3212', 'P3221', 'R32', 'P6', 'P61', 'P65', 'P62', 'P64', 'P63', 'P622', 'P6122', 'P6522', 'P6222', 'P6422', 'P6322', 'P23', 'F23', 'I23', 'P213', 'I213', 'P432', 'P4232', 'F432', 'F4132', 'I432', 'P4332', 'P4132', 'I4132'};
crystalSystems = {'Triclinic', 'Monoclinic', 'Monoclinic', 'Monoclinic', 'Orthorhombic', 'Orthorhombic', 'Orthorhombic', 'Orthorhombic', 'Orthorhombic', 'Orthorhombic', 'Orthorhombic', 'Orthorhombic', 'Orthorhombic', 'Tetragonal', 'Tetragonal', 'Tetragonal', 'Tetragonal', 'Tetragonal', 'Tetragonal', 'Tetragonal', 'Tetragonal', 'Tetragonal', 'Tetragonal', 'Tetragonal', 'Tetragonal', 'Tetragonal', 'Tetragonal', 'Tetragonal', 'Tetragonal', 'Trigonal', 'Trigonal', 'Trigonal', 'Trigonal', 'Trigonal', 'Trigonal', 'Trigonal', 'Trigonal', 'Trigonal', 'Trigonal', 'Trigonal', 'Hexagonal', 'Hexagonal', 'Hexagonal', 'Hexagonal', 'Hexagonal', 'Hexagonal', 'Hexagonal', 'Hexagonal', 'Hexagonal', 'Hexagonal', 'Hexagonal', 'Hexagonal', 'Cubic', 'Cubic', 'Cubic', 'Cubic', 'Cubic', 'Cubic', 'Cubic', 'Cubic', 'Cubic', 'Cubic', 'Cubic', 'Cubic', 'Cubic'};
pointGroups = {'1', '2', '2', '2', '222', '222', '222', '222', '222', '222', '222', '222', '222', '4', '4', '4', '4', '4', '4', '422', '422', '422', '422', '422', '422', '422', '422', '422', '422', '3', '3', '3', '3', '312', '321', '312', '321', '312', '321', '321', '6', '6', '6', '6', '6', '6', '622', '622', '622', '622', '622', '622', '23', '23', '23', '23', '23', '432', '432', '432', '432', '432', '432', '432', '432'};
bravaisTypes = {'P', 'P', 'P', 'C', 'P', 'P', 'P', 'P', 'C', 'C', 'F', 'I', 'I', 'P', 'P', 'P', 'P', 'I', 'I', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'I', 'I', 'P', 'P', 'P', 'R', 'P', 'P', 'P', 'P', 'P', 'P', 'R', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'F', 'I', 'P', 'I', 'P', 'P', 'F', 'F', 'I', 'P', 'P', 'I'};
    end

    methods
        function obj = SpaceGroup(searchID)
            
            if nargin==0
                tblidx = 1;
            elseif isnumeric(searchID)
                tblidx = find(...
                    obj.spaceGroupNumbers==searchID,1,'first');
            elseif ischar(searchID)
                tblidx = find(...
                    strcmpi(searchID,obj.inputNames),1,'first');
            else
                tblidx = [];
            end
            assert(~isempty(tblidx),'did not find space group');
            
            obj.spaceGroupNumber = obj.spaceGroupNumbers(tblidx);
            obj.inputName = obj.inputNames{tblidx};
            obj.fullName = obj.fullNames{tblidx};
            obj.crystalSystem = obj.crystalSystems{tblidx};
            obj.bravaisType = obj.bravaisTypes{tblidx};
            pointGroupName = obj.pointGroups{tblidx};
            obj.PointGroup = geom.lattice.PointGroup(pointGroupName);
            
        end
        
        function tf = checkProperties(obj,a,b,c,alpha,beta,gamma)
            % see table 6-6 on page 300 in Rupp book
            switch lower(obj.crystalSystem)
                case 'triclinic'
                    % nothing to check
                    tf = true;
                case 'monoclinic'
                    tf = alpha==90 & gamma==90;
                case 'orthorhombic'
                    tf = alpha==90 & beta==90 & gamma==90;
                case 'tetragonal'
                    tf = a==b & alpha==90 & beta==90 & gamma==90;
                case {'trigonal','hexagonal'}
                    tf = a==b & alpha==90 & beta==90 & gamma==120;
                case 'cubic'
                    tf = a==b & b==c & alpha==90 & beta==90 & gamma==90;
                otherwise
                    error('this really shouldn''t happen');
            end
            
        end
        
    end
end

