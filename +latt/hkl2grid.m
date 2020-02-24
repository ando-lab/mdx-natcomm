function Iout = hkl2grid(hin,kin,lin,Iin,P,SG)

if nargin < 6 || isempty(SG)
    symexp = false;
else
    symexp = true;
end

Iout = NaN*ones(P.N);

if ~symexp
    [n1,n2,n3] = P.frac2ind(hin,kin,lin,false);
    isincl = n1 >= 1 & n2 >= 1 & n3 >= 1 & ...
        n1 <= P.N(1) & n2 <= P.N(2) & n3 <= P.N(3);
    ind = sub2ind(P.N,n1(isincl),n2(isincl),n3(isincl));
    Iout(ind) = Iin(isincl);
else % symmetry-expand
    
    Ops = SG.LaueGroup.generalPositions;
    hkl0 = [hin(:),kin(:),lin(:)];
    
    for n=1:numel(Ops)
        hkl = hkl0*Ops(n).r;
        [n1,n2,n3] = P.frac2ind(hkl(:,1),hkl(:,2),hkl(:,3),false);
        isincl = n1 >= 1 & n2 >= 1 & n3 >= 1 & ...
            n1 <= P.N(1) & n2 <= P.N(2) & n3 <= P.N(3);
        ind = sub2ind(P.N,n1(isincl),n2(isincl),n3(isincl));
        Iout(ind) = Iin(isincl);
    end
    
end

end