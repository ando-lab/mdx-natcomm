function t = lineSphereIntersection(lineOrigin,lineDirection,sphereCenter,sphereRadius,tmax)

if nargin==4
    tmax = Inf;
end

t = NaN*ones(size(lineOrigin,1),2);
lineOrigin = lineOrigin - repmat(sphereCenter,size(lineOrigin,1),1);

a = dot(lineDirection,lineDirection);
b = 2*(lineDirection(1)*lineOrigin(:,1) + ...
       lineDirection(2)*lineOrigin(:,2) + ...
       lineDirection(3)*lineOrigin(:,3));
c = sum(lineOrigin.^2,2) - sphereRadius*sphereRadius;
disc = b.*b - 4*a*c;
ainv = 1/a;

for j=1:size(lineOrigin,1)
    if disc(j) >= 0
        t(j,:) = -0.5*ainv*([b(j),b(j)] + [1,-1]*sqrt(disc(j)));
    end % disc==0, one solution
end
t(abs(t)>tmax) = NaN;

end