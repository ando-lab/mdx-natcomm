classdef ImportMTZ < util.propertyValueConstructor
    %MTZIMPORTER  
    properties
        % some options
        convertPhaseToComplexF = false;
        translateColumnNames = true;
        
        % some standard column names from http://legacy.ccp4.ac.uk/html/mtzformat.html
        %
        % I omitted M/ISYM, FPH<n>, DP, DPH<n>, SIGDP, SIGFPH<n>,
        % SIGDPH<n>, HLA,HLB,HLC,HLD, etc.
        %
        % Here are some examples of colunm names that I've encountered:
        %
        % REFMAC4 (via ccp4i2): 
        %
        % H, K, L, FREER, FP, SIGFP, FC, PHIC, FC_ALL, PHIC_ALL, FWT, PHWT,
        % DELFWT, PHDELWT, FOM, PHCOMB, FAN, PHAN, DELFAN, PHDELAN,
        % HLACOMB, HLBCOMB, HLCCOMB, HLDCOMB, FC_ALL_LS, PHIC_ALL_LS
        % 
        % CTRUNCATE:
        %
        % H, K, L, F, SIGF, IMEAN, SIGIMEAN
        %
        % XDS --> POINTLESS --> AIMLESS
        %
        % H, K, L, M/ISYM, BATCH, I, SIGI, FRACTIONCALC, XDET, YDET, ROT, LP, FLAG
        %
        % Note: sometimes you see F(+), F(-), etc. This has not been
        % implemented yet.
        
        %standardColNames = {'H','K','L','BATCH','I','SIGI','FRACTIONCALC','IMEAN','SIGIMEAN','FP','FC','SIGFP','PHIC','PHIM','FOM','WT','FREE','FreeR_flag','FC_ALL','PHIC_ALL'};
        mtz2mdx = struct(...
            'H','h',...
            'K','k',...
            'L','l',...
            'FREE','isFree',...
            'FreeR_flag','isFree',...
            'FREER','isFree',...
            'BATCH','batch',...
            'I','Iobs',...
            'SIGI','sigmaIobs',...
            'IMEAN','Iobs',...
            'SIGIMEAN','sigmaIobs',...
            'FP','Fobs',...
            'SIGFP','sigmaFobs',...
            'F','Fobs',...
            'SIGF','sigmaFobs',...
            'FC','Fc',...
            'PHIC','phiFc',...
            'FC_ALL','FcAll',...
            'PHIC_ALL','phiFcAll');
        
    end
    
    methods
        function obj = ImportMTZ(varargin)
            %MTZIMPORTER
            obj@util.propertyValueConstructor(varargin{:});
        end
        
        function [hklTable,Basis,SpaceGroup] = run(obj,fileName)
            [mtz,colNames,c] = io.mtz.read(fileName);
            SpaceGroup = symm.SpaceGroup(c.spaceGroupNumber);
            Basis = latt.Basis(c.a,c.b,c.c,c.alpha,c.beta,c.gamma);
            
            hklTable = array2table(mtz);
            hklTable.Properties.VariableDescriptions = colNames;
            
            if obj.translateColumnNames
    
                isIncl = ismember(colNames,fieldnames(obj.mtz2mdx));
                newCols = cellfun(@(n) obj.mtz2mdx.(n),colNames(isIncl),'Uni',0);
                assert(numel(unique(newCols))==numel(newCols),'redundant column mapping?');
            
                hklTable = renamevars(hklTable,hklTable.Properties.VariableNames(isIncl),newCols);
            
            else
                hklTable = renamevars(hklTable,hklTable.Properties.VariableNames,colNames);
            end
            
            if obj.convertPhaseToComplexF
                error('not implemented yet');
            end
        end
    end
end

