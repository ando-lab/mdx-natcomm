function [sx,sy,sz] = hkl2s(h,k,l,Crystal,Spindle,phi)
Rc = Crystal.rotationMatrix;
Oc = Crystal.orthogonalizationMatrix;
if length(phi)==1
    Rs = Spindle.rotationMatrix(phi);
    
    B = Rs*Rc*Oc;
    Binv = inv(B);
    sx = h*Binv(1,1) + k*Binv(2,1) + l*Binv(3,1);
    sy = h*Binv(1,2) + k*Binv(2,2) + l*Binv(3,2);
    sz = h*Binv(1,3) + k*Binv(2,3) + l*Binv(3,3);
else % assume there is one phi for each of h,k,l
    sx = zeros(size(h));
    sy = zeros(size(k));
    sz = zeros(size(l));
    for j=1:length(phi)
        Rs = Spindle.rotationMatrix(phi(j));
        B = Rs*Rc*Oc;
        Binv = inv(B);
        sx(j) = h(j)*Binv(1,1) + k(j)*Binv(2,1) + l(j)*Binv(3,1);
        sy(j) = h(j)*Binv(1,2) + k(j)*Binv(2,2) + l(j)*Binv(3,2);
        sz(j) = h(j)*Binv(1,3) + k(j)*Binv(2,3) + l(j)*Binv(3,3);
    end
end

end