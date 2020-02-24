classdef ScatteringFactor
    % SCATTERINGFACTOR - calculate spherical scattering factors of atoms,
    % including incoherent and coherent parts
    %
    % Note that H coefficients are those of bonded hydrogen, not atomic hydrogen 
    % (see private/gauss_fits.m)
    
    properties(SetAccess = private)
        atomName
        coeffs_coh
        coeffs_incoh
    end
    
    properties(Constant = true)
        nGauss = 4; % number of gaussians used in the approximation
    end
    
    properties(Constant = true, Access = private)
        %                    a1        b1        a2        b2        a3        b3        a4        b4         c
        coeffs_coh_H    = [ 0.15513  38.58412   0.45182   9.29027   0.28463  20.46901   0.10691   2.59639   0.00119];
        coeffs_incoh_H  = [-0.01666  78.34007  -0.18635  65.57034  -0.68394  27.06857  -0.11323   9.37823   0.99993];
        coeffs_coh_C    = [ 1.64910  38.79730   2.50350  14.18380   0.41149   1.10034   1.28055   0.44630   0.15407];
        coeffs_incoh_C  = [-0.96147  69.25535  -2.95605  19.30025  -0.83316   2.17020  -1.21382   0.83637   5.97046];
        coeffs_coh_N    = [ 1.69158  31.60101   3.17716  11.10210   0.33005   4.80316   1.56126   0.42304   0.23730];
        coeffs_incoh_N  = [-1.51297  49.91666  -3.06972  15.48366  -0.61218   3.46102  -1.73276   0.77338   6.93353];
        coeffs_coh_O    = [ 1.01195  31.47584   3.13890  12.44928   2.05913   5.46894   1.55362   0.31742   0.23571];
        coeffs_incoh_O  = [-1.31543  43.08746  -3.71940  13.05113  -0.98355   4.37218  -1.88096   0.64385   7.90290];
        coeffs_coh_P    = [ 1.68552  64.51201   3.99744  26.48583   6.53149   1.88896   1.68108   0.49789   1.10261];
        coeffs_incoh_P  = [-2.05335 100.80698  -3.26298  25.23137  -5.59975   2.58398  -3.15324   0.57943  14.08141];
        coeffs_coh_S    = [ 1.85079  52.67786   4.95031  21.63097   7.07400   1.44924   1.79104   0.13284   0.33190];
        coeffs_incoh_S  = [-2.12480  83.74962  -4.10117  21.49352  -5.35044   2.24776  -3.40558   0.55746  14.99166];
        coeffs_coh_Na   = [1.10611 130.67271   3.00509   9.11321   4.91188   3.37837   1.27678   0.32888   0.69248];
        coeffs_incoh_Na = [-1.03573 283.73739  -2.23973  14.34056  -5.19136   4.32041  -2.21449   0.48411  10.68622];
        coeffs_coh_Cl   = [1.68209  47.45921   6.20226  18.57122   7.04719   1.19817   1.58706   0.12449   0.48178];
        coeffs_incoh_Cl = [-2.20870  70.81993  -4.95948  17.86395  -5.23199   1.93036  -3.53565   0.50602  15.94432];
        coeffs_coh_Ar   = [1.63657  44.27479   7.34806  15.73812   1.05051   1.52953   6.54987   0.84947   1.41343 ];
        coeffs_incoh_Ar = [-2.25880  63.39193  -5.88349  14.26711  -5.35694   1.62014  -3.47233   0.41776  16.97830 ];
        coeffs_coh_K    = [1.00650 221.45871   8.28544  12.92495   0.83205  46.97673   7.44488   0.78492   1.43009];
        coeffs_incoh_K  = [-1.34306 299.57084  -4.60243  22.37030  -4.15582   6.33013  -7.20635   0.83210  17.34038];
        coeffs_coh_Ca   = [0.90971 185.87019   1.68553  89.83236   8.63106  10.51589   7.38916   0.66946   1.38414 ];
        coeffs_incoh_Ca = [-2.13660 273.83657  -3.14101  24.59351  -5.55582   6.86469  -7.40535   0.76950  18.24808];
        coeffs_coh_Fe   = [2.29909  78.37701   3.51246  15.36448  11.81670   4.78062   7.42070   0.30421   0.94061 ];
        coeffs_incoh_Fe = [-2.32463 161.37395  -3.91559  18.29992  -9.17871   3.64320  -8.38488   0.50098  23.81334];
        coeffs_coh_Co   = [2.32608  73.11073   3.94429  13.70573  12.39253   4.31555   7.42865   0.27554   0.89672];
        coeffs_incoh_Co = [-2.33070 151.45801  -4.17923  16.94664  -9.60074   3.34341  -8.59881   0.47693  24.71981 ];
        coeffs_coh_Cd   = [1.54149  97.38631   3.86245  27.30179  18.30986   7.17992  19.29216   0.60685   4.99263 ];
        coeffs_incoh_Cd = [-2.55464 155.48986  -5.13124  23.95549 -15.20921   3.41365 -17.67659   0.44209  40.58075];
        
    end
    
    methods
        function obj = ScatteringFactor(atomName)
            
            switch lower(atomName)
                case 'h'
                    obj.atomName = 'H';
                    obj.coeffs_coh = obj.coeffs_coh_H;
                    obj.coeffs_incoh = obj.coeffs_incoh_H;
                case 'c'
                    obj.atomName = 'C';
                    obj.coeffs_coh = obj.coeffs_coh_C;
                    obj.coeffs_incoh = obj.coeffs_incoh_C;
                case 'n'
                    obj.atomName = 'N';
                    obj.coeffs_coh = obj.coeffs_coh_N;
                    obj.coeffs_incoh = obj.coeffs_incoh_N;
                case 'o'
                    obj.atomName = 'O';
                    obj.coeffs_coh = obj.coeffs_coh_O;
                    obj.coeffs_incoh = obj.coeffs_incoh_O;
                case 'p'
                    obj.atomName = 'P';
                    obj.coeffs_coh = obj.coeffs_coh_P;
                    obj.coeffs_incoh = obj.coeffs_incoh_P;
                case 's'
                    obj.atomName = 'S';
                    obj.coeffs_coh = obj.coeffs_coh_S;
                    obj.coeffs_incoh = obj.coeffs_incoh_S;
                case 'na'
                    obj.atomName = 'Na';
                    obj.coeffs_coh = obj.coeffs_coh_Na;
                    obj.coeffs_incoh = obj.coeffs_incoh_Na;
                case 'cl'
                    obj.atomName = 'Cl';
                    obj.coeffs_coh = obj.coeffs_coh_Cl;
                    obj.coeffs_incoh = obj.coeffs_incoh_Cl;
                case 'ar'
                    obj.atomName = 'Ar';
                    obj.coeffs_coh = obj.coeffs_coh_Ar;
                    obj.coeffs_incoh = obj.coeffs_incoh_Ar;
                case 'k'
                    obj.atomName = 'K';
                    obj.coeffs_coh = obj.coeffs_coh_K;
                    obj.coeffs_incoh = obj.coeffs_incoh_K;
                case 'ca'
                    obj.atomName = 'Ca';
                    obj.coeffs_coh = obj.coeffs_coh_Ca;
                    obj.coeffs_incoh = obj.coeffs_incoh_Ca;
                case 'fe'
                    obj.atomName = 'Fe';
                    obj.coeffs_coh = obj.coeffs_coh_Fe;
                    obj.coeffs_incoh = obj.coeffs_incoh_Fe;
                case 'co'
                    obj.atomName = 'Co';
                    obj.coeffs_coh = obj.coeffs_coh_Co;
                    obj.coeffs_incoh = obj.coeffs_incoh_Co;
                case 'cd'
                    obj.atomName = 'Cd';
                    obj.coeffs_coh = obj.coeffs_coh_Cd;
                    obj.coeffs_incoh = obj.coeffs_incoh_Cd;
                otherwise
                    error('did not recognize atom name');
            end
        end
        
        function sf = f_coh(obj,s,b_add)
            % atomic scattering factor (coherent, elastic) in the kinematic
            % approximation (without anomalous corrections)
            %
            % default angle-variable is s = 2*sin(th)/lambda = q/(2*pi).
            %
            % optional second argument adds an isotropic B-factor
            
            coeffs = obj.coeffs_coh;
            if nargin==3 && ~isempty(b_add)
                coeffs = add_b_to_coeffs(coeffs,b_add);
            end
            
            % Gaussian fits were done using sin(th)/lambda, which is 0.5*s;
            sf = sumOfGaussians(0.5*s,coeffs);
        end
        function rho = rho_e(obj,r,b_add)
            
            if nargin<3 || isempty(b_add)
                b_add = 0;
            end
            
            sf_coeffs = add_b_to_coeffs(obj.coeffs_coh,b_add);
            % note: adding a b-factor of zero does nothing, unless one is
            % using a constant offset. Then, the number of coefficients
            % increased by 1, representing the constant offset as exp(-0)
            
            ai = sf_coeffs(1:2:(end-1));
            bi = sf_coeffs(2:2:(end));
            
            rho_coeffs = zeros(size(sf_coeffs));
            rho_coeffs(1:2:(end-1)) = ai.*(4*pi./bi).^(3/2);
            rho_coeffs(2:2:end) = 16*pi^2./bi;
            rho = sumOfGaussians(0.5*r,rho_coeffs);
        end
        function sf = I_incoh(obj,s)
            % incoherent scattering intensity
            %
            % default angle-variable is s = 2*sin(th)/lambda = q/(2*pi).
            
            % Gaussian fits were done using sin(th)/lambda, which is 0.5*s;
            sf = sumOfGaussians(0.5*s,obj.coeffs_incoh);
        end
        
    end
    
end

function coeffs = add_b_to_coeffs(coeffs_in,b_add)
ncoeffs = length(coeffs_in);
if floor(ncoeffs/2)~=ncoeffs/2 % odd, has constant offset
    coeffs = [coeffs_in,0]; % represent constant offset as gaussian with b=0
else
    coeffs = coeffs_in;
end
coeffs(2:2:end) = coeffs(2:2:end) + b_add; % add b-factor
end

