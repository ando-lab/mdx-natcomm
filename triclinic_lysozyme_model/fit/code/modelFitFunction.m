function resid = modelFitFunction(I,sigma,Wk,BZgrid,LatticeNetwork,Mass,springConstants)

LatticeNetwork.springConstants = springConstants;
K = LatticeNetwork.Hessian();
L = nm.LatticeModel(Mass,K,BZgrid.N);
modeArray = L.normalModes();

I1 = L.onePhononScattering(modeArray,Wk);
Icalc = cell2mat(shiftdim(I1,-1));

% fit background (constant + linear)
[dh,dk,dl] = BZgrid.grid;
A = ones(prod(BZgrid.N),4);
A(:,2) = dh(:);
A(:,3) = dk(:);
A(:,4) = dl(:);

B = zeros(size(I));

for jj=1:size(Icalc,1)
    x = (A.*repmat(1./sigma(jj,:)',1,size(A,2)))\((I(jj,:) - Icalc(jj,:))./sigma(jj,:))';
    B(jj,:) = A*x;
end

% residual
resid = (I - Icalc - B)./sigma;

end