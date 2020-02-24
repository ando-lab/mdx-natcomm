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
        function [stats,ind] = run(obj,A,B)
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