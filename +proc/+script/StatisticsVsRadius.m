classdef StatisticsVsRadius < util.propertyValueConstructor
    % StatisticsVsRadius -
    properties
        edges
        Basis = latt.Basis.empty(); 
        PeriodicGrid = latt.PeriodicGrid.empty();
        
        fid = 1; % where to print stats. 1 is command window. can also print to log file.
    end
    
    methods
        function obj = StatisticsVsRadius(varargin)
            obj@util.propertyValueConstructor(varargin{:});
        end
        
        function T = tableimport(obj,T_in,cols)
            [x,y,z] = obj.Basis.frac2lab(T_in.(1),T_in.(2),T_in.(3));
            r = sqrt(x.^2 + y.^2 + z.^2);
            
            if ~isempty(obj.edges)
                isIncl = r > min(obj.edges) & r < max(obj.edges);
            else
                isIncl = true;
            end
            
            if nargin < 3 || isempty(cols)
                cols = T_in.Properties.VariableNames(4:end);
            end
            
            A = T_in.(cols{1});
            A = A(isIncl);
            r = r(isIncl);
            ind = ones(size(r));
            
            if numel(cols)==1
                T = table(r,ind,A);
            else %numel(cols)==2
                B = T_in.(cols{2});
                B = B(isIncl);
                T = table(r,ind,A,B);
            end
            
        end
        
        function T = array2table(obj,A,B)
            
            [x,y,z] = latt.LatticeGrid(obj.PeriodicGrid,obj.Basis).grid();
            r = sqrt(x.^2 + y.^2 + z.^2);
            isOneInput = nargin <= 2;
            
            if ~isempty(obj.edges)
                isIncl = r > min(obj.edges) & r < max(obj.edges);
            else
                isIncl = true;
            end
            
            if isOneInput
                isIncl = isIncl & ~isnan(A);
                r = r(isIncl);
                A = A(isIncl);
                ind = find(isIncl(:));
                T = table(r,ind,A);
            else
                isIncl = isIncl & ~isnan(A) & ~isnan(B);
                r = r(isIncl);
                A = A(isIncl);
                B = B(isIncl);
                ind = find(isIncl(:));
                T = table(r,ind,A,B);
            end
            
        end
        
        function [stats,ind] = table2stats(obj,T)
            
            isOneInput = size(T,2)==3;
            
            if isempty(obj.edges)
                rEdges = linspace(0,max(T.r),51)';
            else
                rEdges = obj.edges(:);
            end
            
            rCenters = rEdges(2:end)*0.5 + rEdges(1:(end-1))*0.5;
            
            d3r = (4*pi/3)*(rEdges(2:end).^3 - rEdges(1:(end-1)).^3); % spherical shell volume
            
            [~,ind] = histc(T.r,rEdges);
            isIncl = ind > 0 & ind < numel(rEdges);
            ind(~isIncl) = NaN;
            
            m = accumarray(ind(isIncl),1,[numel(rCenters),1],@sum);
            Aav = accumarray(ind(isIncl),T.A(isIncl),[numel(rCenters),1],@sum);
            Aav = Aav./m;
            
            Asub = T.A(isIncl) - Aav(ind(isIncl));
            Avar = accumarray(ind(isIncl),Asub.^2,[numel(rCenters),1],@sum);
            Avar = Avar./m;
            
            if isOneInput
                stats = table(rCenters,Aav,Avar,m,d3r,'VariableNames',{'r','av','var','n','d3r'});
            else
                Bav = accumarray(ind(isIncl),T.B(isIncl),[numel(rCenters),1],@sum);
                Bav = Bav./m;
                Bsub = T.B(isIncl) - Bav(ind(isIncl));
                Bvar = accumarray(ind(isIncl),Bsub.^2,[numel(rCenters),1],@sum);
                Bvar = Bvar./m;
                
                ABcov = accumarray(ind(isIncl),Bsub.*Asub,[numel(rCenters),1],@sum);
                ABcov = ABcov./m;
                
                cc = ABcov./sqrt(Avar.*Bvar);
                
                stats = table(rCenters,Aav,Avar,Bav,Bvar,ABcov,cc,m,d3r,'VariableNames',{'r','av1','var1','av2','var2','cov','cc','n','d3r'});
            end
            
        end
        
        function [stats,ind] = run(obj,varargin)
            % refactored (Feb 3 2022)
            T = obj.array2table(varargin{:});
            [stats,row_index] = obj.table2stats(T);
            ind = NaN*ones(size(varargin{1}));
            ind(T.ind(~isnan(row_index))) = row_index(~isnan(row_index));
        end
        
        function [stats,ind] = run_v0(obj,A,B)
            % old version, before refactoring
            [x,y,z] = latt.LatticeGrid(obj.PeriodicGrid,obj.Basis).grid();
            r = sqrt(x.^2 + y.^2 + z.^2);
            
            isOneInput = nargin <= 2;
            
            if isOneInput
                isIncl = ~isnan(A);
            else
                isIncl = ~isnan(A) & ~isnan(B);
            end
            
            if isempty(obj.edges)
                rEdges = linspace(0,max(r(isIncl)),51)';
            else
                rEdges = obj.edges(:);
            end
            rCenters = rEdges(2:end)*0.5 + rEdges(1:(end-1))*0.5;
            
            d3r = (4*pi/3)*(rEdges(2:end).^3 - rEdges(1:(end-1)).^3); % spherical shell volume
            
            %4*pi*diff(rEdges).*rCenters.^2; % approximate shell volume
            
            [~,ind] = histc(r,rEdges);
            isIncl = isIncl & ind > 0 & ind < numel(rEdges);
            
            m = accumarray(ind(isIncl),1,[numel(rCenters),1],@sum);
            Aav = accumarray(ind(isIncl),A(isIncl),[numel(rCenters),1],@sum);
            Aav = Aav./m;
            
            Asub = A(isIncl) - Aav(ind(isIncl));
            Avar = accumarray(ind(isIncl),Asub.^2,[numel(rCenters),1],@sum);
            Avar = Avar./m;
            
            if isOneInput
                stats = table(rCenters,Aav,Avar,m,d3r,'VariableNames',{'r','av','var','n','d3r'});
            else
                Bav = accumarray(ind(isIncl),B(isIncl),[numel(rCenters),1],@sum);
                Bav = Bav./m;
                Bsub = B(isIncl) - Bav(ind(isIncl));
                Bvar = accumarray(ind(isIncl),Bsub.^2,[numel(rCenters),1],@sum);
                Bvar = Bvar./m;
                
                ABcov = accumarray(ind(isIncl),Bsub.*Asub,[numel(rCenters),1],@sum);
                ABcov = ABcov./m;
                
                cc = ABcov./sqrt(Avar.*Bvar);
                
                stats = table(rCenters,Aav,Avar,Bav,Bvar,ABcov,cc,m,d3r,'VariableNames',{'r','av1','var1','av2','var2','cov','cc','n','d3r'});
            end
            
            
        end
        
    end
end