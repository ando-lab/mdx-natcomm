classdef SymmetryOperator < symm.AffineTransformation
    % SYMMETRYOPERATOR
    %
    % some info at: http://www.crystallography.fr/mathcryst/pdf/Gargnano/Aroyo_Gargnano_1.pdf
    
    properties(Dependent = true)
        % fourMatrix % four matrix
        xyzForm
        charForm
    end
    
    methods
        function obj = SymmetryOperator(varargin)
            %SYMMETRYOPERATOR Construct an instance of this class
            obj@symm.AffineTransformation(varargin{:});
            
        end
        
        function val = mtimes(obj1,obj2)
            val = mtimes@symm.AffineTransformation(obj1,obj2);
            if isa(obj1,'symm.SymmetryOperator') && isa(obj2,'symm.SymmetryOperator')
               val = symm.SymmetryOperator(val.r,val.t);
            end
        end
        
        function val = inv(obj)
            val = inv@symm.AffineTransformation(obj);
            val = symm.SymmetryOperator(val.r,val.t);
        end
        
        function disp(obj)
            for j=1:length(obj)
                disp(['   [' obj(j).xyzForm ']']);
            end
        end

        function isEq = eq(obj,obj2)
            
            % commented out because it slows things down:
            % assert(isa(obj,class(obj2)),'mixed types not supported');
            
            if numel(obj) == 1 && numel(obj2) == 1
                mat1 = obj.fourMatrix;
                mat2 = obj2.fourMatrix;
                
                % use common denominator (integer comparison)
                %[n,d] = rat([mat1,mat2]);
                %common_denominator = prod(d(:)); % don't bother finding lcd
                common_denominator = 48; % HACK ALERT: this should always be the case for crystallographic space groups?
                mat1 = round(mat1*common_denominator);
                mat2 = round(mat2*common_denominator);
                
                isEq = all(mat1(:) == mat2(:));
                
            elseif all(size(obj) == size(obj2))
                isEq = false(size(obj));
                for j=1:numel(isEq)
                    isEq(j) = eq(obj(j),obj2(j));
                end
            elseif numel(obj)==1 && ~isempty(obj2)
                isEq = false(size(obj2));
                for j=1:numel(isEq)
                    isEq(j) = eq(obj,obj2(j));
                end
            elseif numel(obj2)==1 && ~isempty(obj)
                isEq = eq(obj2,obj);
            else
                error('size mismatch');
            end
        end
        
        function isNotEq = ne(obj,obj2)
            isNotEq = ~eq(obj,obj2);
        end

        function [val,sortOrder] = sort(obj,varargin)
             % sort alphabetically using the xyzForm string
             [~,sortOrder] = sort({obj.xyzForm},varargin{:});
             val = obj(sortOrder);
        end
        
        function str = get.xyzForm(obj)
            str = mat2xyzstr(obj.r,obj.t);
        end
        
        function str = get.charForm(obj)
            coeffStr = cellfun(@(v) strtrim(rats(v)),...
                num2cell([obj.r,obj.t]),'UniformOutput',false);
            str = ['[' strjoin(coeffStr(1,:),',') ';' ...
                strjoin(coeffStr(2,:),',') ';' ...
                strjoin(coeffStr(3,:),',') ']'];
        end
    end
    
    methods(Static)
        function PG = space2point(SG)
            % convert a space group to a point group.
            % basically, just set t = [0;0;0] for all operators.
            PG = SG;
            for j=1:length(PG)
                PG(j).t = PG(j).t*0;
            end
        end    
            
        function [SymGroup] = generators2group(SymGen)
            maxDepth = 4; % don't multiply more than 4 operators together
            
            % first, make sure the SymGen are as expected
            for j=1:length(SymGen)
                % [numer,denom] = rat(SymGen(j).t); % a bit slow, replaced
                % with the code below
                numer = round(SymGen(j).t * 48);
                denom = 48; % HACK ALERT: the translational part is always an integer multiple of 1/48
                SymGen(j).t = mod(numer,denom)./denom; % map translational part back to [0,1) interval
            end
            
            SymGroup = SymGen;
            
            for d = 1:maxDepth
                newOps = symm.SymmetryOperator.empty();
                
                % try to make new operators by combining pairs of existing
                % ones
                for j=1:length(SymGroup)
                    for k=1:length(SymGroup)
                        opjk = SymGroup(j)*SymGroup(k);
                        
                        %[numer,denom] = rat(opjk.t); % a bit slow, replaced
                        % with the code below
                        numer = round(opjk.t * 48);
                        denom = 48; % HACK ALERT: the translational part is always an integer multiple of 1/48
                        
                        opjk.t = mod(numer,denom)./denom; % map translational part back to [0,1) interval

                        if ~ismember(opjk,SymGroup) && (isempty(newOps) || ~ismember(opjk,newOps))
                            newOps = [newOps,opjk];
                        end
                    end
                end
                if isempty(newOps) % no new operators were generated
                    %disp(d);
                    return;
                elseif d == maxDepth
                    error('max depth reached');
                end
                SymGroup = [SymGroup,newOps]; %add to group and try again                
            end

        end
    end
end

function str = mat2xyzstr(r,t)
coordStrings = {'x','y','z'};

outStr = {{'','','',''},{'','','',''},{'','','',''}};

for j=1:3
    for k=1:3
        ajk = r(j,k);
        if ajk
            ajkstr = strtrim(rats(ajk));
            outStr{j}{k} = [ajkstr '*' coordStrings{k}];
        end
    end
    aj = t(j,1);
    if aj
        ajstr = strtrim(rats(aj));
        outStr{j}{4} = ajstr;
    end
    outStr{j} = strjoin(outStr{j}(~cellfun(@isempty,outStr{j})),'+');
end

% clean it up
str = strrep(strjoin(outStr,','),'+-','-');
str = regexprep(str,'(?<!\d)1\*','');
end

