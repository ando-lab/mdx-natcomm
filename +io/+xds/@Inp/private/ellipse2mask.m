function m = ellipse2mask(untrusted_ellipse,nx,ny)

% default values for converting vertices to polygons
th = linspace(0,2*pi,97);
x_circ = (cos(th) + 1)/2; y_circ = (sin(th) + 1)/2;

m = true(nx,ny);
for j=1:size(untrusted_ellipse,1)
    x1 = untrusted_ellipse(j,1);
    x2 = untrusted_ellipse(j,2);
    y1 = untrusted_ellipse(j,3);
    y2 = untrusted_ellipse(j,4);
    
    % the xds convention here is unknown, but the following works most of
    % the time
    x = (x_circ*(x2-x1) + x1);
    y = (y_circ*(y2-y1) + y1);
    
    m = m & not(poly2mask(y,x,nx,ny)); 
end

end
