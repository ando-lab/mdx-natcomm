function resid = adpFitFunction(Umeas,P,M,Kfun,varargin)
% calculate the residual between a set of measured B-factors (or ADPs) and
% those predicted by an elastic model

% Umeas is a cell array of ADPs (or Uiso)
numPts = size(P,1)/3;
assert(size(Umeas,1) == 3*numPts);
assert(size(Umeas,2) == 3);

K = Kfun(varargin{:});
L = nm.LatticeModel(M,K,[1,1,1]);

modeArray = L.normalModes();
covMat = L.covarianceMatrix(modeArray,[0,0,0]);

Ucalc = adpCalc(P,covMat);

resid = Umeas-Ucalc;

end