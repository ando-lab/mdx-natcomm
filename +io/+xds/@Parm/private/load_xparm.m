function p = load_xparm(fn)
%  file definition from: 
%
%  http://xds.mpimf-heidelberg.mpg.de/html_doc/xds_files.html#XPARM.XDS
%
%  The file comprises 11 lines of data values that code for
% 
%     1.  The first line only contains one keyword, XPARM.XDS,
%     distinguishing the new file format from the older versions (released
%     before March 30, 2013) which do not have this line.
%
%     2.  Starting image number (STARTING_FRAME=), spindle angle at start
%     (STARTING_ANGLE=), oscillation range, and laboratory coordinates of
%     the rotation axis.
%
%     3.  Wavelength (Å) and laboratory coordinates of the incident beam
%     wavevector.
%
%     4.  Space group number and unit cell parameters (Å and degrees).
%
%     5.  Laboratory coordinates of the unit cell a-axis of the unrotated
%     crystal (at spindle dial 0 degrees).
%
%     6.  Laboratory coordinates of the unit cell b-axis of the unrotated
%     crystal (at spindle dial 0 degrees).
%
%     7.  Laboratory coordinates of the unit cell c-axis of the unrotated
%     crystal (at spindle dial 0 degrees).
%
%     8.  Number of detector segments, number of "fast" pixels (NX=),
%     number of "slow" pixels (NY=) in a data image, and pixel sizes (mm)
%     (QX=, QY=) along the "fast" and "slow" directions, respectively.
%
%     9.  Specifies the origin of the detector coordinate system with
%     respect to the laboratory system by 3 numbers, ORGX, ORGY, F. ORGX
%     and ORGY are in pixel units, F is in mm.
%
%     10. Laboratory coordinates of the unit vector along the detector
%     X-axis.
%
%     11. Laboratory coordinates of the unit vector along the detector
%     Y-axis.
%
%     12. Laboratory coordinates of the unit vector along the detector
%     normal.
% 
%   Now, for each detector segment the following two lines of information
%   are provided.
%
%     The 5 numbers of this line, iseg x1 x2 y1 y2, define the pixel
%     numbers IX,IY belonging to segment #iseg as x1<=IX<=x2, y1<=IY<=y2.
%
%     The 9 numbers of this line, ORGXS ORGYS FS EDS(:,1) EDS(:,2),
%     describe origin and orientation of segment #iseg with respect to the
%     detector coordinate system.
if nargin==0 || isempty(fn)
    fn = 'XPARM.XDS';
end

fid = fopen(fn);
header_text = fgetl(fid);
cols = textscan(fid,'%f %f %f %f %f %f %f %f %f',...
    'Delimiter',' ','MultipleDelimsAsOne',1);
fclose(fid);
h.data = [cols{:}];
             
p.header_text = header_text;

p.starting_frame = h.data(1,1);
p.starting_angle = h.data(1,2);
p.oscillation_range = h.data(1,3);
p.rotation_axis = h.data(1,4:6);

p.lambda = h.data(2,1);
p.incident_wavevector = h.data(2,2:4);

p.space_group_number = h.data(3,1);
p.a = h.data(3,2);
p.b = h.data(3,3);
p.c = h.data(3,4);
p.alpha = h.data(3,5);
p.beta = h.data(3,6);
p.gamma = h.data(3,7);

p.a_axis = h.data(4,1:3);

p.b_axis = h.data(5,1:3);

p.c_axis = h.data(6,1:3);

p.detector.nseg = h.data(7,1);
p.detector.nx = h.data(7,2);
p.detector.ny = h.data(7,3);
p.detector.qx = h.data(7,4);
p.detector.qy = h.data(7,5);

p.detector.orgx = h.data(8,1);
p.detector.orgy = h.data(8,2);
p.detector.f = h.data(8,3);

p.detector.ed(:,1) = h.data(9,1:3);

p.detector.ed(:,2) = h.data(10,1:3);

p.detector.ed(:,3) = h.data(11,1:3);

% need to loop over number of segments

for j=1:p.detector.nseg
    
  % iseg x1 x2 y1 y2
    p.detector.segments(j).iseg     = h.data(12 + 2*(j-1),1);
    p.detector.segments(j).x1       = h.data(12 + 2*(j-1),2);
    p.detector.segments(j).x2       = h.data(12 + 2*(j-1),3);
    p.detector.segments(j).y1       = h.data(12 + 2*(j-1),4);
    p.detector.segments(j).y2       = h.data(12 + 2*(j-1),5);
    
  % ORGXS ORGYS FS EDS(:,1) EDS(:,2)
    p.detector.segments(j).orgxs    = h.data(13 + 2*(j-1),1);
    p.detector.segments(j).orgys    = h.data(13 + 2*(j-1),2);
    p.detector.segments(j).fs       = h.data(13 + 2*(j-1),3);
    p.detector.segments(j).eds(:,1) = h.data(13 + 2*(j-1),4:6);
    p.detector.segments(j).eds(:,2) = h.data(13 + 2*(j-1),7:9);

end

end