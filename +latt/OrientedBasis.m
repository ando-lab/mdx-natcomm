classdef OrientedBasis < latt.Basis
    properties
        orientationVector = [0,0,0]
    end
    properties(Dependent)
        orientationMatrix
    end
    methods
        function obj = OrientedBasis(a,b,c,alpha,beta,gamma,orientationVector)
            obj@latt.Basis(a,b,c,alpha,beta,gamma);
            if nargin>6
                obj.orientationVector = orientationVector;
            end
        end
        function value = get.orientationMatrix(obj)
            value = rvec2mat(obj.orientationVector);
        end
        function [newObj,err] = invert(obj)
            
            % after inverting, need to realign again.
            theseAxes = obj.orientationMatrix*obj.orthogonalizationMatrix;
            targetAxes = inv(theseAxes)';
            newBasis = invert@latt.Basis(obj);
            [newObj,err] = obj.realign(newBasis,targetAxes(:,1),targetAxes(:,2),targetAxes(:,3));
        end
        function [X,Y,Z] = frac2lab(obj,x,y,z)
            M = obj.orientationMatrix*obj.orthogonalizationMatrix;
            [X,Y,Z] = latt.Basis.transformCoordinates(M,x,y,z);
        end
        function [x,y,z] = lab2frac(obj,X,Y,Z)
            M = obj.orientationMatrix*obj.orthogonalizationMatrix;
            [x,y,z] = latt.Basis.transformCoordinates(inv(M),X,Y,Z);
        end
    end
    methods(Static)
        function [newObj,err] = realign(BasisObj,a1,a2,a3)
            % [newObj,err] = realign(BasisObj,a1,a2,a3)
            %
            % a1, a2, a3 are the target axes
            %
            % one of the axes may be [], signalling that only the other two
            % axes are fitted. It may be necessary to fix a sign ambiguity
            % later so that the rotation matrix has a determinant of 1
            
            % the variable dim signals which dimension was omitted in the
            % fit. dimOmitted = 0 means that no dimensions were omitted
            dimOmitted = 0;
            if isempty(a1) || all(a1==0)
                a1 = zeros(3,1);
                dimOmitted = 1;
            elseif isempty(a2) || all(a2==0)
                a2 = zeros(3,1);
                dimOmitted = 2;
            elseif isempty(a3) || all(a3==0)
                a3 = zeros(3,1);
                dimOmitted = 3;
            end
            
            theseAxes = BasisObj.orthogonalizationMatrix;
            
            targetAxes = [a1(:),a2(:),a3(:)];
            
            [u,~,v] = svd(theseAxes*targetAxes','econ');
            R = v*u';
            
            if dimOmitted && det(R) < 0
                R(:,dimOmitted) = -R(:,dimOmitted);
            end
            
            assert(det(R)>0,'coordinate system is not right-handed?');
            err = R*theseAxes - targetAxes;
            if dimOmitted
                err(:,dimOmitted) = 0;
            end
            
            l = rmat2vec(R);
            
            newObj = latt.OrientedBasis(...
                BasisObj.a,BasisObj.b,BasisObj.c,...
                BasisObj.alpha,BasisObj.beta,BasisObj.gamma,...
                l);

        end
    end
end

function l = rmat2vec(R)
% Convert the rotation matrix into axis-angle form. Note: this is kind of a
% hack, and could probably be simplified. But it works.

tolEye = 1E-12;

if all(abs(diag(R)-1) < tolEye)
    % R is the identity matrix
    l = [0,0,0];
    return;
end

% find eigen-decomposition of the rotation matrix, R
[v,d] = eig(R);

% find the eigenvalue which is real
[~,ix] = min(abs(diag(d) - 1));
% get the corresponding eigenvector, call it lz
lz = v(:,ix);

% calculate an orthonormal basis lx, ly, lz
ly = sum(v(:,(1:3)~=ix),2); ly = ly/sqrt(dot(ly,ly));
lx = cross(ly,lz);

% rotate lx and find the angle in the x-y plane
lxr = R*lx;
ang = (180/pi)*atan2(dot(lxr,ly),dot(lxr,lx));

% store the rotation matrix as a row vector, whose length is
% the angle in degrees, and direction is the rotation axis
l = lz(:)'*ang;

end

function R = rvec2mat(v)
v = v(:)'; % force row vector;
dphi_degrees = norm(v);
if dphi_degrees==0
    v = [1,0,0];
else
    v = v/norm(v);
end
dphi_radians = dphi_degrees*pi/180;

% rotation matrix (rodrigues formula)
R = cos(dphi_radians)*eye(3) + sin(dphi_radians)*...
    [0, -v(3), v(2); v(3),0,-v(1); -v(2), v(1),0] + ...
    (1-cos(dphi_radians))*(v'*v);
end
