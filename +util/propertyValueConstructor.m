classdef propertyValueConstructor
    properties(Dependent = true, Access = protected)
        propsToUpdate
    end
    methods
        function obj = propertyValueConstructor(varargin)
            if ~isempty(varargin)
                [thisobj,fn,fnExcl] = parse_inputs(obj.propsToUpdate,class(obj),varargin{:});
                if ~isempty(fnExcl)
                    error('in constructor, the following arguments were not assigned: %s',strjoin(fnExcl));
                end
            else
                fn = {}; % do nothing
            end
            for j=1:length(fn)
                obj.(fn{j}) = thisobj.(fn{j});
            end
        end
        
        function propNames = get.propsToUpdate(obj)
            mc = metaclass(obj);
            allProps = mc.PropertyList;
            ispublic = strcmp({allProps.SetAccess},'public') & ...
                not([allProps.Dependent]);
            propNames = {allProps(ispublic).Name};
        end
        
        function newobj = update(obj,varargin)
            [thisobj,fn] = parse_inputs(obj.propsToUpdate,class(obj),varargin{:});
            newobj = obj;
            
            for j=1:length(fn)
                if ~isempty(thisobj.(fn{j}))
                    % if the property is a class with an update method,
                    % call that instead of assigning it a new value
                    if isa(newobj.(fn{j}),'util.propertyValueConstructor')
                        newobj.(fn{j}) = update(newobj.(fn{j}),thisobj.(fn{j}));
                    else % just assign the value
                        newobj.(fn{j}) = thisobj.(fn{j});
                    end
                end
            end
        end
    end
end


function [thisobj,fn,fnExcl] = parse_inputs(propsToUpdate,myClass,varargin)
if length(varargin)>1 && ischar(varargin{1})
    % input is in the form of name, value pairs.
    thisobj = cell2struct(varargin(2:2:end)',varargin(1:2:(end-1)));
    fn = fieldnames(thisobj);
    isIncl = ismember(fn,propsToUpdate);
    fnExcl = fn(~isIncl);
    fn = fn(isIncl);
elseif length(varargin)==1 && isstruct(varargin{1})
    thisobj = varargin{1};
    fn = fieldnames(thisobj);
    isIncl = ismember(fn,propsToUpdate);
    fnExcl = fn(~isIncl);
    fn = fn(isIncl);
elseif length(varargin{1})==1 && isa(varargin{1},myClass)
    thisobj = varargin{1};
    fn = propsToUpdate;
    fnExcl = {};
elseif length(varargin{1})==1
    % check to see if class(varargin{1}) is a subclass of myClass
    mc = meta.class.fromName(myClass);
    % the following checks only the first in the list of superclasses
    % it may not work if there is more than one superclass (rare)
    if ~isempty(mc.SuperclassList) && isa(varargin{1},mc.SuperclassList(1).Name)
        thisobj = varargin{1};
        fn = thisobj.propsToUpdate;
        fnExcl = {};
    else
        error('%s is not a subclass of %s',mc.Name,myClass);
    end
else
    error('did not recognize input arguments');
end
end

