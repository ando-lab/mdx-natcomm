function [h,k,l] = ewald2hkl(Source,Crystal,Spindle,phi,smax)
% list of miller indices for all reciprocal lattice cells
% with s(h,k,l) < smax that have one or more edges intersecting the ewald
% sphere.

Rc = Crystal.rotationMatrix;
Oc = Crystal.orthogonalizationMatrix;
Rs = Spindle.rotationMatrix(phi);
B = Rs*Rc*Oc;
Binv = inv(B);

% coordinates of ewald sphere
sphereRadius = 1/Source.wavelength;
sphereCenter = -Source.wavevector;

% define range of h, l, and l to calculate

% planes bounding the sphere:
hlim = (sort(sphereCenter*B(:,1) + ...
    [1,-1]*sphereRadius/sqrt(dot(Binv(1,:),Binv(1,:)))));
klim = (sort(sphereCenter*B(:,2) + ...
    [1,-1]*sphereRadius/sqrt(dot(Binv(2,:),Binv(2,:)))));
llim = (sort(sphereCenter*B(:,3) + ...
    [1,-1]*sphereRadius/sqrt(dot(Binv(3,:),Binv(3,:)))));

% apply resolution limits to bounding planes
hmax = smax/sqrt(Binv(1,:)*Binv(1,:)');
kmax = smax/sqrt(Binv(2,:)*Binv(2,:)');
lmax = smax/sqrt(Binv(3,:)*Binv(3,:)');
hlim(abs(hlim) > hmax) = sign(hlim(abs(hlim) > hmax))*hmax;
klim(abs(klim) > kmax) = sign(klim(abs(klim) > kmax))*kmax;
llim(abs(llim) > lmax) = sign(llim(abs(llim) > lmax))*lmax;

% define reciprocal unit cell edges
hedges = (ceil(hlim(1) - 0.5) + 0.5):(floor(hlim(2) + 0.5) - 0.5);
kedges = (ceil(klim(1) - 0.5) + 0.5):(floor(klim(2) + 0.5) - 0.5);
ledges = (ceil(llim(1) - 0.5) + 0.5):(floor(llim(2) + 0.5) - 0.5);

% find intersections for edges in the l-direction
[hgrid,kgrid,lgrid] = ndgrid(hedges,kedges,0);

lineOrigin = [hgrid(:),kgrid(:),lgrid(:)]*Binv;
lineDirection = [0,0,1]*Binv;

lsph = lineSphereIntersection(lineOrigin,lineDirection,...
    sphereCenter,sphereRadius,lmax);

hasInt1 = ~isnan(lsph(:,1));
hasInt2 = ~isnan(lsph(:,2));

% l index of adjacent cells
lsph = round(lsph);

% miller indices of cells adjacent to intersecting edges
indl = [hgrid(hasInt1)+.5,kgrid(hasInt1)+.5,lsph(hasInt1,1);
        hgrid(hasInt1)-.5,kgrid(hasInt1)+.5,lsph(hasInt1,1);
        hgrid(hasInt1)+.5,kgrid(hasInt1)-.5,lsph(hasInt1,1);
        hgrid(hasInt1)-.5,kgrid(hasInt1)-.5,lsph(hasInt1,1);
        hgrid(hasInt2)+.5,kgrid(hasInt2)+.5,lsph(hasInt2,2);
        hgrid(hasInt2)-.5,kgrid(hasInt2)+.5,lsph(hasInt2,2);
        hgrid(hasInt2)+.5,kgrid(hasInt2)-.5,lsph(hasInt2,2);
        hgrid(hasInt2)-.5,kgrid(hasInt2)-.5,lsph(hasInt2,2)];

% apply symmetry limits
indl = indl(sum((indl*Binv).^2,2) <= smax^2,:);     

% find intersections for edges in the k-direction
[hgrid,kgrid,lgrid] = ndgrid(hedges,0,ledges);

lineOrigin = [hgrid(:),kgrid(:),lgrid(:)]*Binv;
lineDirection = [0,1,0]*Binv;

ksph = lineSphereIntersection(lineOrigin,lineDirection,...
    sphereCenter,sphereRadius,kmax);

hasInt1 = ~isnan(ksph(:,1));
hasInt2 = ~isnan(ksph(:,2));

% k index of adjacent cells
ksph = round(ksph);

% miller indices of cells adjacent to intersecting edges
indk = [hgrid(hasInt1)+.5,ksph(hasInt1,1),lgrid(hasInt1)+.5;
        hgrid(hasInt1)+.5,ksph(hasInt1,1),lgrid(hasInt1)-.5;
        hgrid(hasInt1)-.5,ksph(hasInt1,1),lgrid(hasInt1)+.5;
        hgrid(hasInt1)-.5,ksph(hasInt1,1),lgrid(hasInt1)-.5;
        hgrid(hasInt2)+.5,ksph(hasInt2,2),lgrid(hasInt2)+.5;
        hgrid(hasInt2)+.5,ksph(hasInt2,2),lgrid(hasInt2)-.5;
        hgrid(hasInt2)-.5,ksph(hasInt2,2),lgrid(hasInt2)+.5;
        hgrid(hasInt2)-.5,ksph(hasInt2,2),lgrid(hasInt2)-.5];

% apply symmetry limits
indk = indk(sum((indk*Binv).^2,2) <= smax^2,:);    

% find intersections for edges in the h-direction
[hgrid,kgrid,lgrid] = ndgrid(0,kedges,ledges);

lineOrigin = [hgrid(:),kgrid(:),lgrid(:)]*Binv;
lineDirection = [1,0,0]*Binv;

hsph = lineSphereIntersection(lineOrigin,lineDirection,...
    sphereCenter,sphereRadius,hmax);

hasInt1 = ~isnan(hsph(:,1));
hasInt2 = ~isnan(hsph(:,2));

% h index of adjacent cells
hsph = round(hsph);

% miller indices of cells adjacent to intersecting edges
indh = [hsph(hasInt1,1),kgrid(hasInt1)+.5,lgrid(hasInt1)+.5;
        hsph(hasInt1,1),kgrid(hasInt1)+.5,lgrid(hasInt1)-.5;
        hsph(hasInt1,1),kgrid(hasInt1)-.5,lgrid(hasInt1)+.5;
        hsph(hasInt1,1),kgrid(hasInt1)-.5,lgrid(hasInt1)-.5;
        hsph(hasInt2,2),kgrid(hasInt2)+.5,lgrid(hasInt2)+.5;
        hsph(hasInt2,2),kgrid(hasInt2)+.5,lgrid(hasInt2)-.5;
        hsph(hasInt2,2),kgrid(hasInt2)-.5,lgrid(hasInt2)+.5;
        hsph(hasInt2,2),kgrid(hasInt2)-.5,lgrid(hasInt2)-.5];

% apply symmetry limits
indh = indh(sum((indh*Binv).^2,2) <= smax^2,:);

% unique list of miller indices with edges that intersect the ewald sphere
hkl = unique([indh;indk;indl],'rows');

% assign output variables
h = hkl(:,1);
k = hkl(:,2);
l = hkl(:,3);
end