function [c1,c2] = calc_reflection_center(s0, m2, f, ed, p0star)
%
% columns of c1, c2: phi,x,y,p1,p2,p3
%
% where p1,p2,p3 are the reciprocal space coordinates once rotated into the
% diffraction condition
%
% xparm = structure containing xds parameters
%
% p0star = some arbitrary point of interest in reciprocal space of the
%          unrotated crystal (inverse angstrom)
%
% c1 = [phi,x,y]; <- first set of solutions
% c2 = [phi,x,y]; <- second set of solutions
%
% partly vectorized on May 18, 2015
%   (several operations could be more succinctly written as matrix
%    multiplications)
%
% removed x0, and y0 (x and y are now relative to detector origin)
%
% moved to scatterbrain project: Feb 17, 2017
%
% removed xparm. now s0, m2, f, ed are passed as arguments

%s0 = xparm.incident_wavevector;
%m2 = xparm.rotation_axis;
m1 = cross(m2,s0); m1 = m1/sqrt(dot(m1,m1));
m3 = cross(m1,m2);

%f = xparm.detector.f;
%d1 = xparm.detector.ed(:,1)';
%d2 = xparm.detector.ed(:,2)';
%d3 = xparm.detector.ed(:,3)';

d1 = ed(:,1)';
d2 = ed(:,2)';
d3 = ed(:,3)';

p0star_dot_p0star = ...
    p0star(:,1).*p0star(:,1) + ...
    p0star(:,2).*p0star(:,2) + ...
    p0star(:,3).*p0star(:,3);

p0star_dot_m2 = ...
    p0star(:,1)*m2(1) + ...
    p0star(:,2)*m2(2) + ...
    p0star(:,3)*m2(3);

p0star_dot_m1 = ...
    p0star(:,1)*m1(1) + ...
    p0star(:,2)*m1(2) + ...
    p0star(:,3)*m1(3);

p0star_dot_m3 = ...
    p0star(:,1)*m3(1) + ...
    p0star(:,2)*m3(2) + ...
    p0star(:,3)*m3(3);

rho_squared = p0star_dot_p0star - p0star_dot_m2.*p0star_dot_m2;
%rho = sqrt(rho_squared); % distance to rot. axis

pstar_dot_m3 = (-0.5*p0star_dot_p0star - p0star_dot_m2*dot(s0,m2))/dot(s0,m3);
pstar_dot_m2 = p0star_dot_m2;

% check not in blind region

isoutside = rho_squared < (pstar_dot_m3.*pstar_dot_m3) | ...
    p0star_dot_p0star > 4*dot(s0,s0);
pstar_dot_m1 = sqrt(rho_squared - pstar_dot_m3.*pstar_dot_m3);
pstar_dot_m1(isoutside) = NaN;
cos_phi = (pstar_dot_m1.*p0star_dot_m1 + ...
                   pstar_dot_m3.*p0star_dot_m3)./rho_squared;
sin_phi = (pstar_dot_m1.*p0star_dot_m3 - ...
                   pstar_dot_m3.*p0star_dot_m1)./rho_squared;
phi = atan2(sin_phi,cos_phi)*180/pi;
pstar = [m1(1)*pstar_dot_m1 + m2(1)*pstar_dot_m2 + m3(1)*pstar_dot_m3,...
    m1(2)*pstar_dot_m1 + m2(2)*pstar_dot_m2 + m3(2)*pstar_dot_m3,...
    m1(3)*pstar_dot_m1 + m2(3)*pstar_dot_m2 + m3(3)*pstar_dot_m3];
S_dot_d1 = (pstar(:,1) + s0(1))*d1(1) + ...
           (pstar(:,2) + s0(2))*d1(2) + ...
           (pstar(:,3) + s0(3))*d1(3);
S_dot_d2 = (pstar(:,1) + s0(1))*d2(1) + ...
           (pstar(:,2) + s0(2))*d2(2) + ...
           (pstar(:,3) + s0(3))*d2(3);
S_dot_d3 = (pstar(:,1) + s0(1))*d3(1) + ...
           (pstar(:,2) + s0(2))*d3(2) + ...
           (pstar(:,3) + s0(3))*d3(3);       

x = f*S_dot_d1./S_dot_d3;
y = f*S_dot_d2./S_dot_d3;
c1 = [phi,x,y,pstar];

pstar_dot_m1 = -1*sqrt(rho_squared - pstar_dot_m3.*pstar_dot_m3);
pstar_dot_m1(isoutside) = NaN;
cos_phi = (pstar_dot_m1.*p0star_dot_m1 + ...
                   pstar_dot_m3.*p0star_dot_m3)./rho_squared;
sin_phi = (pstar_dot_m1.*p0star_dot_m3 - ...
                   pstar_dot_m3.*p0star_dot_m1)./rho_squared;
phi = atan2(sin_phi,cos_phi)*180/pi;
pstar = [m1(1)*pstar_dot_m1 + m2(1)*pstar_dot_m2 + m3(1)*pstar_dot_m3,...
    m1(2)*pstar_dot_m1 + m2(2)*pstar_dot_m2 + m3(2)*pstar_dot_m3,...
    m1(3)*pstar_dot_m1 + m2(3)*pstar_dot_m2 + m3(3)*pstar_dot_m3];

S_dot_d1 = (pstar(:,1) + s0(1))*d1(1) + ...
           (pstar(:,2) + s0(2))*d1(2) + ...
           (pstar(:,3) + s0(3))*d1(3);
S_dot_d2 = (pstar(:,1) + s0(1))*d2(1) + ...
           (pstar(:,2) + s0(2))*d2(2) + ...
           (pstar(:,3) + s0(3))*d2(3);
S_dot_d3 = (pstar(:,1) + s0(1))*d3(1) + ...
           (pstar(:,2) + s0(2))*d3(2) + ...
           (pstar(:,3) + s0(3))*d3(3);       

x = f*S_dot_d1./S_dot_d3;
y = f*S_dot_d2./S_dot_d3;
c2 = [phi,x,y,pstar];

end
