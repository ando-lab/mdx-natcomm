classdef Integrate
    methods(Static)
        
        function [wedge,suggestedValues] = read(fileName)
            if nargin==0 || isempty(fileName)
                fileName = 'INTEGRATE.LP';
            elseif isdir(fileName)
                fileName = fullfile(fileName,'INTEGRATE.LP');
            end
            [wedge,suggestedValues] = read_integrate_lp(fileName);
        end
        
        function [wedge,suggestedValues] = profiles(fileName)
            if nargin==0 || isempty(fileName)
                fileName = 'INTEGRATE.LP';
            elseif isdir(fileName)
                fileName = fullfile(fileName,'INTEGRATE.LP');
            end
            [filePath,~,~] = fileparts(fileName);
            xdsInp = io.xds.Inp.read(filePath);
            nptsAlphaBeta = xdsInp.number_of_profile_grid_points_along_alpha_beta;
            nptsGamma = xdsInp.number_of_profile_grid_points_along_gamma;
            
            if isempty(nptsAlphaBeta)
                nptsAlphaBeta = 13; % older xds versions had default of 9, new is 13
            end
            if isempty(nptsGamma)
                nptsGamma = 9;
            end
            
            wedge = read_profiles(fileName,nptsAlphaBeta,nptsGamma);
        end
        
        function dataTable = hkl(fileName)
            if nargin==0 || isempty(fileName)
                fileName = 'INTEGRATE.HKL';
            elseif isdir(fileName)
                fileName = fullfile(fileName,'INTEGRATE.HKL');
            end
            [data,colNames] = read_integrate_hkl(fileName);
            dataTable = array2table(data,'VariableNames',lower(colNames));
        end
        
        function plot(wedge)
            p = [wedge.xparm];
            s = [wedge.stats];
            tab = [wedge.table];
            
            % calculate unit cell volume
            v = zeros(1,length(p));
            for j=1:length(p)
                v(j) = det([p(j).a_axis',p(j).b_axis',p(j).c_axis']);
            end
            
            fh = figure(gcf);clf; fh.Position(3) = 560; fh.Position(4) = 750;
            
            subplot(4,1,1);plot([tab.image],[tab.sigmar]);
            xlabel('Image Number');ylabel('Sigmar (degrees)');title('Mosaicity');
            
            subplot(4,1,2);plot(reshape([s.crystal_rotation]',3,[])','o-');
            xlabel('Wedge');ylabel('Angular deviation (deg)');
            legend({'X','Y','Z'},'location','best');title('Crystal Rotation');
            
            subplot(4,1,3);plot([tab.image],[tab.scale]);
            xlabel('Image Number');ylabel('Scale Factor');title('Scaling');
            
            subplot(4,1,4);plot(v,'o-');
            xlabel('Wedge');ylabel('Unit Cell Volume (Å^3)');
        end
    end
end