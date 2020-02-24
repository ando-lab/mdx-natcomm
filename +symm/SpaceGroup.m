classdef SpaceGroup
    % SpaceGroup - symmetry operators and other information related to a
    % crystallographic space group
    %
    % Implemented as a set of dynamic properties that wraps
    % around a SpaceGroupInfo object
    
    properties
        Info
    end
    
    properties(Dependent = true)
        name
        number
        generators
        generalPositions
        PointGroup
        LaueGroup
    end
    
    methods
        function obj = SpaceGroup(varargin)
            %SpaceGroup Construct an instance of this class
            obj.Info = symm.SpaceGroupInfo(varargin{:});
        end
        function val = get.generators(obj)
            val = cellfun(@symm.SymmetryOperator,obj.Info.generators);
        end
        function val = get.generalPositions(obj)
            val = cellfun(@symm.SymmetryOperator,obj.Info.generalPositions);
        end
        function val = get.name(obj)
            val = obj.Info.name;
        end
        function val = get.PointGroup(obj)
            val = symm.PointGroup(obj.Info.pointGroupNumber);
        end
        function val = get.LaueGroup(obj)
            val = symm.LaueGroup(obj.Info.laueGroupNumber);
        end
        function val = get.number(obj)
            val = obj.Info.number;
        end
        
        function tf = checkLattice(obj,a,b,c,alpha,beta,gamma)
            tf = obj.Info.testLattice(a,b,c,alpha,beta,gamma);
        end

    end
end

