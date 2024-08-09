classdef PointGroup
    %PointGroup - symmetry operators and other information related to a
    %crystallographic point group
    %
    % Implemented as a set of dynamic properties that wraps
    % around a SpaceGroupInfo object
    
    properties
        Info
    end
    
    properties(Dependent = true)
        name
        generators
        generalPositions
    end
    
    methods
        function obj = PointGroup(varargin)
            %PointGroup Construct an instance of this class
            if nargin==2 && strcmpi(varargin{1},'name')
                if ~strncmpi(varargin{2},'P',1)
                    varargin{2} = ['P' varargin{2}];
                end
            elseif nargin==1 && ~isnumeric(varargin{1})
                if ~any(strncmpi(varargin{1},{'P','R'},1))
                    varargin{1} = ['P' varargin{1}];
                end
            end
            obj.Info = symm.SpaceGroupInfo(varargin{:});
            
            assert(obj.Info.isPointGroup,...
                'requested space group is not a point group');
        end
        function val = get.generators(obj)
            val = cellfun(@symm.SymmetryOperator,obj.Info.generators);
        end
        function val = get.generalPositions(obj)
            val = cellfun(@symm.SymmetryOperator,obj.Info.generalPositions);
        end
        function val = get.name(obj)
            val = obj.Info.name(2:end);
        end 
    end
    
end

