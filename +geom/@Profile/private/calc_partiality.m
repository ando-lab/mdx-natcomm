function R = calc_partiality(phi,e1,phi1,phi2,sigmaM,m2)
% Feb 16, 2017: modified to use m2 as input rather than xds_parm

%m2 = xds_parm.rotation_axis;
nrefl = size(phi,1);

% xi = m2 dot e1
xi = sum(repmat(m2,nrefl,1).*e1,2);

a = abs(xi)/(sqrt(2)*sigmaM);

R = 0.5 * (erf(a.*(phi2-phi)) - erf(a.*(phi1-phi)));

end