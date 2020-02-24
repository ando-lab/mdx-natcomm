classdef Basis
    properties
        a
        b
        c
        alpha
        beta
        gamma
    end
    properties(Dependent)
        volume % unit cell volume
        orthogonalizationMatrix
    end
    properties(Dependent,Hidden=true)
        a1
        a2
        a3
        alphaRadian
        betaRadian
        gammaRadian
    end
    methods
        function obj = Basis(a,b,c,alpha,beta,gamma)
            % obj@util.properyValueConstructor(varargin{:});
            assert(nargin==6,'incorrect number of arguments');
            obj.a = a;
            obj.b = b;
            obj.c = c;
            obj.alpha = alpha;
            obj.beta = beta;
            obj.gamma = gamma;
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
        
        function value = get.volume(obj)
            value = det(obj.orthogonalizationMatrix);
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
            
            % volume: see eqn 5-6 on pg. 233 in Rupp book
            vol = obj.a*obj.b*obj.c*...
                sqrt(1 - cosAlpha^2 - cosBeta^2 - cosGamma^2 + ...
                2*cosAlpha*cosBeta*cosGamma);
            
            value = [ obj.c*cosBeta ,...
                obj.c*(cosAlpha - cosBeta*cosGamma)/sinGamma ,...
                vol/(obj.a*obj.b*sinGamma) ];
        end
        
        function value = get.orthogonalizationMatrix(obj)
            % see eqn 5-5 on pg. 233 in Rupp book
            value = [obj.a1',obj.a2',obj.a3'];
        end
        
        function newObj = invert(obj)
            oMinv = inv(obj.orthogonalizationMatrix);
            
            astar = norm(oMinv(1,:));
            bstar = norm(oMinv(2,:));
            cstar = norm(oMinv(3,:));
            
            alphastar = (180/pi)*atan2(...
                norm(cross(oMinv(2,:),oMinv(3,:))),...
                dot(oMinv(2,:),oMinv(3,:)));
            betastar = (180/pi)*atan2(...
                norm(cross(oMinv(1,:),oMinv(3,:))),...
                dot(oMinv(1,:),oMinv(3,:)));
            gammastar = (180/pi)*atan2(...
                norm(cross(oMinv(1,:),oMinv(2,:))),...
                dot(oMinv(1,:),oMinv(2,:)));
            
            newObj = latt.Basis(astar,bstar,cstar,alphastar,betastar,gammastar);
        end
        
        function [X,Y,Z] = frac2lab(obj,x,y,z)
            M = obj.orthogonalizationMatrix;
            [X,Y,Z] = latt.Basis.transformCoordinates(M,x,y,z);
        end
        
        function [x,y,z] = lab2frac(obj,X,Y,Z)
            M = obj.orthogonalizationMatrix;
            [x,y,z] = latt.Basis.transformCoordinates(inv(M),X,Y,Z);
        end
        
        
    end
    
    methods(Static, Access = protected)
        function [X,Y,Z] = transformCoordinates(M,x,y,z)
            X = zeros(size(x));
            Y = zeros(size(y));
            Z = zeros(size(z));
            xyz = [x(:)';y(:)';z(:)'];
            XYZ = M*xyz;
            X(:) = XYZ(1,:);
            Y(:) = XYZ(2,:);
            Z(:) = XYZ(3,:);
        end
    end
end