classdef CombineTables < util.propertyValueConstructor
    % CombineTables - 
    properties
        %hklRef = table(); 
        hklcols1 = {'h','k','l','I','sigma'}; % first 3 must be miller indices (h,k,l)
        hklcols2 = {'h','k','l','I','sigma'};  % first 3 must be miller indices (h,k,l)
        
        postfix1 = '';
        postfix2 = '2'; % string to append to the end of column names (hklcols2{4:end})
        
        map2asu = false; % re-map hkl values of the mtz to the ASU
        Crystal = geom.Crystal(); % only used if map2asu is true
        
        fid = 1; % where to print stats. 1 is command window. can also print to log file.
    end
    
    methods
        function obj = CombineTables(varargin)
            obj@util.propertyValueConstructor(varargin{:});
        end
        
        function [hklOut,ia,ib] = run(obj,hklTable1,hklTable2)
            
            hklTable1 = hklTable1(:,obj.hklcols1);
            hklTable2 = hklTable2(:,obj.hklcols2);

            
            if obj.map2asu
                [hklTable1.(1),hklTable1.(2),hklTable1.(3)] = obj.Crystal.hkl2asu(...
                    hklTable1.(1),hklTable1.(2),hklTable1.(3));
                [hklTable2.(1),hklTable2.(2),hklTable2.(3)] = obj.Crystal.hkl2asu(...
                    hklTable2.(1),hklTable2.(2),hklTable2.(3));
            end
            
            hklArray1 = table2array(hklTable1(:,1:3));
            hklArray2 = table2array(hklTable2(:,1:3));
            
            hklArray = unique([hklArray1;hklArray2],'rows');
            [~,ia] = ismember(hklArray1,hklArray,'rows');
            [~,ib] = ismember(hklArray2,hklArray,'rows');
            
            nrows = size(hklArray,1);
            
            tmp1 = NaN*ones(nrows,size(hklTable1,2)-3);
            tmp1(ia,:) = table2array(hklTable1(:,4:end));
            
            tmp2 = NaN*ones(nrows,size(hklTable2,2)-3);
            tmp2(ib,:) = table2array(hklTable2(:,4:end));
            
            % modify column names
            cols1 = cellfun(@(c) [c,obj.postfix1],...
                obj.hklcols1(4:end),'UniformOutput',false);
            cols2 = cellfun(@(c) [c,obj.postfix2],...
                obj.hklcols2(4:end),'UniformOutput',false);
                        
            hklOut = array2table([hklArray,tmp1,tmp2],'VariableNames',...
                [obj.hklcols1(1:3),cols1,cols2]);
            
        end
    end
end