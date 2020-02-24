function draw_detector(xds_inp)
%%
%figure(1);clf;
%%
x_sq = [0,1,1,0,0]; y_sq = [1,1,0,0,1];
th = linspace(0,2*pi,100);
x_circ = (cos(th) + 1)/2; y_circ = (sin(th) + 1)/2;

x0 = xds_inp.orgx;
y0 = xds_inp.orgy;
nx = xds_inp.nx;
ny = xds_inp.ny;
qx = xds_inp.qx;
qy = xds_inp.qy;
f = xds_inp.detector_distance;
x000 = f*xds_inp.incident_beam_direction(1)/xds_inp.incident_beam_direction(3);
y000 = f*xds_inp.incident_beam_direction(2)/xds_inp.incident_beam_direction(3);

% draw detector face
%patch((x_sq*nx - x0)*qx,(y_sq*ny - y0)*qy,'w');hold all;
patch((x_sq*nx),(y_sq*ny),'w','facecolor','none','edgecolor','r');hold all;

set(gca,'YDir','reverse');
%
for j=1:size(xds_inp.untrusted_rectangle,1)
    x1 = xds_inp.untrusted_rectangle(j,1);
    x2 = xds_inp.untrusted_rectangle(j,2);
    y1 = xds_inp.untrusted_rectangle(j,3);
    y2 = xds_inp.untrusted_rectangle(j,4);
    patch((x_sq*(x2-x1) + x1),(y_sq*(y2-y1) + y1),...
        [1,1,1]*.7,'facecolor','none','edgecolor','g');
end

for j=1:size(xds_inp.untrusted_ellipse,1)
    x1 = xds_inp.untrusted_ellipse(j,1);
    x2 = xds_inp.untrusted_ellipse(j,2);
    y1 = xds_inp.untrusted_ellipse(j,3);
    y2 = xds_inp.untrusted_ellipse(j,4);
    patch((x_circ*(x2-x1) + x1),(y_circ*(y2-y1) + y1),...
        [1,1,1]*.7,'facecolor','none','edgecolor','g');
end
plot(x0,y0,'b+'); hold all;
plot(x000/qx + x0,y000/qx + y0,'ro'); %<-- beam center
axis image;
end