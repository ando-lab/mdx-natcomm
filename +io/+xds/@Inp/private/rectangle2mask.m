function m = rectangle2mask(untrusted_rectangle,nx,ny)
% default values for converting vertices to polygons
x_sq = [0,1,1,0,0]; y_sq = [1,1,0,0,1];
m = true(nx,ny);
for j=1:size(untrusted_rectangle,1)
    x1 = untrusted_rectangle(j,1);
    x2 = untrusted_rectangle(j,2);
    y1 = untrusted_rectangle(j,3);
    y2 = untrusted_rectangle(j,4);

    % the xds convention is that x1 < x < x2, y1 < y < y2
    % while the convention for convention for poly2mask is 
    % x1 < x <= x2, y1 < y <= y2
    x = (x_sq*(x2-x1 - 1) + x1 + 0.5); % was (x_sq*(x2-x1) + x1)
    y = (y_sq*(y2-y1 - 1) + y1 + 0.5); % was (y_sq*(y2-y1) + y1)
    
    m = m & not(poly2mask(y,x,nx,ny)); 
end

end
