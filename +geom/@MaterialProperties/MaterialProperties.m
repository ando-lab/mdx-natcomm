classdef MaterialProperties
    methods(Static)
        function mu = calcMu(material,wavelength)
            % requires wavelength in angstrom
            % returns mu in mm^-1
            switch lower(material)
                case {'si','silicon'}
                    data = abs_si_data();
                case 'air'
                    data = abs_air_data();
                otherwise
                    error('did not recognize material');
            end
            
            mu = interp1q(data(:,1),data(:,2),wavelength);
        end
    end
end