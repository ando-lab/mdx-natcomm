classdef AffineTransformation
    % AffineTransformation (operator)
    %
    % some info at: http://www.crystallography.fr/mathcryst/pdf/Gargnano/Aroyo_Gargnano_1.pdf
    %
    % I = AffineTransformation(r,t);
    
    properties
        r (3,3) = eye(3); % rotation matrix
        t (3,1) = zeros(3,1); % translation vector
    end
    properties(Dependent = true)
        fourMatrix % four matrix
        isPureRotation
        isPureTranslation
        isIdentity
    end
    
    methods
        function obj = AffineTransformation(r,t)
            %AffineTransformation Construct an instance of this class
            if nargin == 0
                return; % accept defaul values
            end
            
            if all(size(r) == [3,3])
                obj.r = r;
                if nargin==2
                    obj.t = t;
                end
            elseif size(r,2) == 4 % four-matrix
                assert(nargin==1,'second argument not expected');
                if size(r,1) == 4 % check proper four matrix
                    assert(all(r(4,:) == [0,0,0,1]),'not a proper 4-matrix?');
                end
                obj.r = r(1:3,1:3);
                obj.t = r(1:3,4);
            else
                error('did not recognize input');
            end
            
        end
        
        function val = get.isPureRotation(obj)
            val = all(obj.t==0);
        end
        
        function val = get.isPureTranslation(obj)
            eyecmp = obj.r==eye(3);
            val = all(eyecmp(:));
        end
        
        function val = get.isIdentity(obj)
            val = obj.isPureRotation & obj.isPureTranslation;
        end
        
        function val = get.fourMatrix(obj)
            val = [obj.r, obj.t; 0,0,0,1];
        end
        
        function val = mtimes(obj1,obj2)
            if isa(obj1,'symm.AffineTransformation') && isa(obj2,'symm.AffineTransformation') % just in case?
                % compose the transformations
                r1 = obj1.r;
                r2 = obj2.r;
                t1 = obj1.t;
                t2 = obj2.t;
                val = symm.AffineTransformation(r1*r2,r1*t2 + t1);
            elseif isnumeric(obj2) % obj2 is a vector, return a vector
                if obj1.isIdentity % short circuit if identity operator
                    val = obj2;
                    return;
                end
                if size(obj2,1)==3 % 3-vector in. return 3-vector
                    if obj1.isPureRotation % for speed
                        val = obj1.r*obj2;
                    else
                        vecNewFull = obj1.fourMatrix*[obj2; ones(1,size(obj2,2))];
                        val = vecNewFull(1:3,:);
                    end
                else
                    assert(size(obj2,1)==4) % 4-vector in. return a 4-vector
                    val = obj1.fourMatrix*obj2;
                end
            elseif isnumeric(obj1) % obj1 is a covector, return a covector
                if obj2.isIdentity % short circuit if identity operator
                    val = obj1;
                    return;
                end
                if size(obj1,2)==3 % 3-vector in. return 3-vector
                    if obj2.isPureRotation % for speed
                        val = obj1*obj2.r;
                    else
                        vecNewFull = [obj1, zeros(size(obj1,1),1)]*obj2.fourMatrix;
                        val = vecNewFull(:,1:3);
                    end
                else
                    assert(size(obj1,2)==4) % 4-vector in. return 4-vector
                    val = obj1*obj2.fourMatrix;
                end
            else
                error('argument type not recognized');
            end
        end
        
        function ATinv = inv(obj)
            % inverse of the operator
            rinv = inv(obj.r);
            tinv = -rinv*obj.t;
            ATinv = symm.AffineTransformation(rinv,tinv);
        end
        
    end
    
end
