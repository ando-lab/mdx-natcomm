function [sx,sy,sz] = hkl2s0(h,k,l,Crystal)
% calculate sx, sy, sz at phi = 0
Rc = Crystal.rotationMatrix;
if isempty(Rc), Rc = eye(3); end % if setting unknown
Oc = Crystal.orthogonalizationMatrix;

B = Rc*Oc;
Binv = inv(B);
sx = h*Binv(1,1) + k*Binv(2,1) + l*Binv(3,1);
sy = h*Binv(1,2) + k*Binv(2,2) + l*Binv(3,2);
sz = h*Binv(1,3) + k*Binv(2,3) + l*Binv(3,3);
end