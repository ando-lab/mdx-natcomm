classdef GaussianAtom < latt.GaussianDensitySum
    %GAUSSIANATOM
    %
    % a,b are the coefficients of the N-gaussian structure factor
    %
    %    $sf = \sum_{i=1}^N a(i)*exp(-b(i)*s^2/4)$
    %
    % or, if N= length(b) and length(a) == N+1, then a(N+1) represent a
    % constant offset
    %
    %    $sf = \sum_{i=1}^N a(i)*exp(-b(i)*s^2/4) + a(N+1)$
    %
    % U is the 3x3 atomic displacement tensor, with Debye waller factor
    %
    %   $T = exp(-2*pi^2*s'*U*s)$
    %
    % or, U could be given as a scalar Uiso, such that
    %
    %   $U = Uiso*eye(3)$
    % In order to add an isotropic B-factor, it is necessary to convert using
    %
    %   $B = 8*pi^2*U$
    %
    % or
    %
    %   $U = B/(8*pi^2)$
    %
    % If the U is omitted, then U=0 is assumed.
    %
    % Note: it is also assumed that $U == U'$.
    %
    % The total structure factor is
    %
    %   $sf = \sum_{i=1}^N a(i)*exp(-b(i)*s^2/4)*exp(-2*pi^2*s'*U*s)$
    %   $sf = \sum_{i=1}^N a(i)*exp(- 2*pi^2*s'*V(j)*s / 4)$
    %
    % where $V(j) = U + (b(j)/(8*pi^2))*eye(3);%
    %
    % The corresponding electron density is
    %
    %   $rho = \sum{i=1}^N a(i)*(exp(-r'*inv(V(i))*r/2)/sqrt((2*pi)^3*det(V(i)))%
    
    properties(SetAccess = immutable)
        %a
        b
        U
    end
    
    properties(Dependent = true)
        Biso
        Uiso
        electrons
    end
    
    methods
        function obj = GaussianAtom(a,b,U)
            %GAUSSIANATOM
            
            assert(nargin>=2) % a and b are required, U is optional
            
            if nargin < 3
                U = 0; % default value of U
            end
            
            N = numel(b);
            assert(N == numel(a) | (N + 1) == numel(a));
            
            if numel(a) == N + 1
                N = N + 1;
                b(N) = 0;
            end
            
            assert(numel(U)==1 | all(size(U)==[3,3]));
            
            isIsotropic = numel(U) == 1;
            
            if isIsotropic
                V = zeros(1,N);
                for j=1:N
                    V(j) = U + b(j)/(8*pi^2);
                end
            else
                U = 0.5*(U + U'); % make sure its symmetric
                V = zeros(3,3,N);
                for j=1:N
                    V(:,:,j) = U + (b(j)/(8*pi^2))*eye(3);
                end
            end
            obj@latt.GaussianDensitySum(a,V);
            
            obj.U = U;
            obj.b = b;
            
        end
        
        
        function newGA = mtimes(obj1,obj2)
            % I havent implement matrix multiplication yet. scalars only
            assert(numel(obj1)==1 & numel(obj2)==1,'scalars only');
            
            [a1,a2] = ndgrid(obj1.a,obj2.a);
            [b1,b2] = ndgrid(obj1.b,obj2.b);
            anew = a1(:).*a2(:);
            bnew = b1(:) + b2(:);
            if (obj1.isIsotropic && obj2.isIsotropic) || ...
                    ~(obj1.isIsotropic && obj2.isIsotropic)
                Unew = obj1.U + obj2.U;
            elseif obj1.isIsotropic % then obj2 is not
                Unew = eye(3)*obj1.U + obj2.U;
            else % obj2 is isotropic, obj1 is not
                Unew = obj1.U + eye(3)*obj2.U;
            end
            newGA = latt.GaussianAtom(anew,bnew,Unew);
        end
        
        function val = get.Uiso(obj)
            if numel(obj.U)==1
                val = obj.U;
            else
                val = trace(obj.U)/3;
            end
        end
        
        function val = get.Biso(obj)
            val = 8*pi^2*obj.Uiso;
        end
        
        function val = get.electrons(obj)
            % this is equal to F(0) = sum(a(i))
            val = sum(obj.a);
        end
        
        function rho = electronDensity(obj,x,y,z)
            
            if all(obj.U(:)==0) && any(obj.b==0)
                error('FT of a constant is a delta function. Add a B-factor!');
            end
            
            rho = electronDensity@latt.GaussianDensitySum(obj,x,y,z);
            
        end
        
        function Uadd = fix(obj,dmin)
            Uadd = zeros(size(obj.U));

            if numel(obj.U)==1
                if obj.U < dmin
                    Uadd = dmin - obj.U;
                end
                
            else
                [v,d] = eig(obj.U,'vector');
                
                if any(d < dmin)
                   dp = lsqnonlin(@(x) [100*sum(x),x],...
                            [0,0,0],dmin-d,[Inf,Inf,Inf],...
                            optimset('Display','off'));
                   Uadd = v*diag(dp)*v';
                end
            end
        end
        
        function newobj = transform(obj,A)
            if isa(A,'symm.AffineTransformation')
                A = A.r; % get the rotation part
            end
            assert(all(size(A)==[3,3]));
            
            isRotInv = (abs(abs(det(A))-1) < 10*eps);
            
            if obj.isIsotropic && isRotInv
                % nothing needs to be done
                newobj = obj;
                return;
            end
            
            if ~isRotInv
                newobj = transform@latt.GaussianDensitySum(obj,A);
                warning('operator was not a rotation or inversion - returned generalized electron density');
            else
                % the b parts do not transform. only need to modify U
                newobj = obj.resetU(A*obj.U*A');
            end
            
        end
        
        function GA = resetU(obj,U)
            GA = latt.GaussianAtom(obj.a,obj.b,U);
        end
        
        function GA = addU(obj,U)
            GA = latt.GaussianAtom(obj.a,obj.b,obj.U*eye(3) + U);
        end
    end
    methods(Static)
        function obj = loadobj(s)
           % force construct on load for compatibility with earlier
           % versions
           obj = latt.GaussianAtom(s.a,s.b,s.U);
        end
        
        function [GA,Elements,elindex] = pdbImport(Atoms,AnisoTemps)
            %
            % Interpretation of temperature factors
            %
            % Debye-waller factor:
            %
            %  T = exp(-2*pi^2*s'*U*s)
            %  T = exp(-2*pi^2*s^2*trace(U)/3)
            %  T = exp(-2*pi^2*s^2*B/(8*pi^2)) = exp(-s^2*B/4)
            %
            % Probability density for displacement:
            %
            %  P(u) = exp(-u'*inv(U)*u/2)/sqrt((2*pi)^3*det(U))
            %
            %  U -> Uiso*I
            %
            %  P(u) = exp(-u^2/(2*Uiso))/(2*pi*Uiso)^(3/2)
            %  P(u) = exp(-u^2/(2*B/(8*pi^2)))/(2*pi*B/(8*pi^2))^(3/2)
            %  P(u) = exp(-4*pi^2*u^2/B)*((4*pi)/B)^(3/2)
            
            if nargin < 2 || isempty(AnisoTemps)
                AnisoTemps = table([],[],[],[],[],[],[],...
                    'VariableNames',{'serial','u00','u11','u22','u01','u02','u12'});
            end
            
            AnisoTemps = AnisoTemps(ismember(AnisoTemps.serial,Atoms.serial),:);
            
            % all element types present
            [elname,~,elindex] = unique(Atoms.element);
            
            % look up element information
            Elements = cellfun(@model.atom.Element,elname);
            ScatteringFactors = [Elements.ScatteringFactor];

            % tabulate coeffs
            all_coeffs = reshape([ScatteringFactors.coeffs_coh],[],length(ScatteringFactors))';

            % reshape coeffs
            a = all_coeffs(:,1:2:end);
            b = all_coeffs(:,2:2:end);
            
            [ia,locb] = ismember(AnisoTemps.serial,Atoms.serial);
            assert(all(ia),'some atom records are missing?');

            Uiso = Atoms.tempFactor/(8*pi^2);
            
            u11 = NaN*ones(size(Atoms,1),1);
            u22 = NaN*ones(size(Atoms,1),1);
            u33 = NaN*ones(size(Atoms,1),1);
            u12 = NaN*ones(size(Atoms,1),1);
            u13 = NaN*ones(size(Atoms,1),1);
            u23 = NaN*ones(size(Atoms,1),1);
            
            u11(locb) = AnisoTemps.u00*1E-4;
            
            u22(locb) = AnisoTemps.u11*1E-4;
            u33(locb) = AnisoTemps.u22*1E-4;
            u12(locb) = AnisoTemps.u01*1E-4;
            u13(locb) = AnisoTemps.u02*1E-4;
            u23(locb) = AnisoTemps.u12*1E-4;
            
            GA = latt.GaussianAtom.empty();
            for j=1:size(Atoms,1)
                aj = a(elindex(j),:);
                bj = b(elindex(j),:);
                if isnan(u11(j))
                    Uj = Uiso(j);
                else
                    Uj = [u11(j),u12(j),u13(j);u12(j),u22(j),u23(j);u13(j),u23(j),u33(j)];
                    
                end
                GA(j,1) = latt.GaussianAtom(aj,bj,Uj);
            end
            
        end
    end
end

