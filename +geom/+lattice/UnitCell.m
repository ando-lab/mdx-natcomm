classdef UnitCell 
    % UnitCell
    %
    % note: this is an old piece of code that's used mainly by proc.
    % It will be replaced by latt.Basis / latt.OrientedBasis
    properties
        a
        b
        c
        alpha
        beta
        gamma
    end
    properties(Dependent)
        vCell % unit cell volume
        orthogonalizationMatrix
        a1
        a2
        a3
        a1star
        a2star
        a3star
    end
    properties(Dependent,Hidden=true)
        alphaRadian
        betaRadian
        gammaRadian
    end
    methods
        function obj = UnitCell(a,b,c,alpha,beta,gamma)
            
            if nargin==6
                obj.a = a;
                obj.b = b;
                obj.c = c;
                obj.alpha = alpha;
                obj.beta = beta;
                obj.gamma = gamma;
            elseif nargin==1 && length(a)==6
                obj.a = a(1);
                obj.b = a(2);
                obj.c = a(3);
                obj.alpha = a(4);
                obj.beta = a(5);
                obj.gamma = a(6);
            else
                error('incorrect number of arguments');
            end
        end
        
        function value = get.alphaRadian(obj)
            value = obj.alpha*pi/180;
        end
        
        function value = get.betaRadian(obj)
            value = obj.beta*pi/180;
        end
        
        function value = get.gammaRadian(obj)
            value = obj.gamma*pi/180;
        end
        
        function value = get.vCell(obj)
            % see eqn 5-6 on pg. 233 in Rupp book
            cosAlpha = cos(obj.alphaRadian);
            cosBeta = cos(obj.betaRadian);
            cosGamma = cos(obj.gammaRadian);
            value = obj.a*obj.b*obj.c*...
                sqrt(1 - cosAlpha^2 - cosBeta^2 - cosGamma^2 + ...
                2*cosAlpha*cosBeta*cosGamma);
        end
        
        function value = get.a1(obj)
            value = [ obj.a ,...
                0     ,...
                0     ];
        end
        
        function value = get.a2(obj)
            cosGamma = cos(obj.gammaRadian);
            sinGamma = sin(obj.gammaRadian);
            
            value = [ obj.b*cosGamma ,...
                obj.b*sinGamma ,...
                0 ];
        end
        
        function value = get.a3(obj)
            cosAlpha = cos(obj.alphaRadian);
            cosBeta = cos(obj.betaRadian);
            cosGamma = cos(obj.gammaRadian);
            sinGamma = sin(obj.gammaRadian);
            
            value = [ obj.c*cosBeta ,...
                obj.c*(cosAlpha - cosBeta*cosGamma)/sinGamma ,...
                obj.vCell/(obj.a*obj.b*sinGamma) ];
        end
        
        function value = get.orthogonalizationMatrix(obj)
            % see eqn 5-5 on pg. 233 in Rupp book
            value = [obj.a1',obj.a2',obj.a3'];
        end
        
        function value = get.a1star(obj)
            value = 2*pi*cross(obj.a2,obj.a3)/obj.vCell;
        end
        
        function value = get.a2star(obj)
            value = 2*pi*cross(obj.a3,obj.a1)/obj.vCell;
        end
        
        function value = get.a3star(obj)
            value = 2*pi*cross(obj.a1,obj.a2)/obj.vCell;
        end
    end
end