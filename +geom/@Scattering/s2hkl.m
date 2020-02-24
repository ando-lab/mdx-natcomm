function [h,k,l] = s2hkl(sx,sy,sz,Crystal,Spindle,phi)
Rc = Crystal.rotationMatrix;
if isempty(Rc), Rc = eye(3); end % if setting unknown

Oc = Crystal.orthogonalizationMatrix;

if nargin>4
    Rs = Spindle.rotationMatrix(phi);
else
    Rs = eye(3);
end

B = Rs*Rc*Oc;

h = sx*B(1,1) + sy*B(2,1) + sz*B(3,1);
k = sx*B(1,2) + sy*B(2,2) + sz*B(3,2);
l = sx*B(1,3) + sy*B(2,3) + sz*B(3,3);
end