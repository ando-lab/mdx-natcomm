function xparm = header2xparm(parm_geom,parm_detector,parm_segment)

% now, convert into xparm format
xparm.header_text = 'INTEGRATE';

xparm.starting_frame = parm_geom.starting_angle_and_frame(2);
xparm.starting_angle = parm_geom.starting_angle_and_frame(1);
xparm.oscillation_range = parm_geom.oscillation_range;
xparm.rotation_axis = parm_geom.rotation_axis;

xparm.lambda = parm_geom.x_ray_wavelength;

d = parm_geom.incident_beam_direction;
d = d/sqrt(dot(d,d));

xparm.incident_wavevector = d/xparm.lambda;

xparm.space_group_number = parm_geom.space_group_number;
xparm.a = parm_geom.unit_cell_constants(1);
xparm.b = parm_geom.unit_cell_constants(2);
xparm.c = parm_geom.unit_cell_constants(3);
xparm.alpha = parm_geom.unit_cell_constants(4);
xparm.beta = parm_geom.unit_cell_constants(5);
xparm.gamma = parm_geom.unit_cell_constants(6);

xparm.a_axis = []; % not present in INTEGRATE.wedge header lines
xparm.b_axis = [];
xparm.c_axis = [];

xparm.detector.nseg = parm_detector.number_of_detector_segments;
xparm.detector.nx = parm_detector.nx_ny_qx_qy(1);
xparm.detector.ny = parm_detector.nx_ny_qx_qy(2);
xparm.detector.qx = parm_detector.nx_ny_qx_qy(3);
xparm.detector.qy = parm_detector.nx_ny_qx_qy(4);

xparm.detector.orgx = parm_detector.orgx_orgy(1);
xparm.detector.orgy = parm_detector.orgx_orgy(2);
xparm.detector.f = parm_detector.detector_distance;

xparm.detector.ed(:,1) = parm_detector.direction_of_detector_x_axis;
xparm.detector.ed(:,2) = parm_detector.direction_of_detector_y_axis;
xparm.detector.ed(:,3) = cross(xparm.detector.ed(:,1),xparm.detector.ed(:,2));

% need to loop over number of segments
for j=1:xparm.detector.nseg
    xparm.detector.segments(j).iseg     = j;
    xparm.detector.segments(j).x1       = parm_segment(j).seg_x1_x2_y1_y2(1);
    xparm.detector.segments(j).x2       = parm_segment(j).seg_x1_x2_y1_y2(2);
    xparm.detector.segments(j).y1       = parm_segment(j).seg_x1_x2_y1_y2(3);
    xparm.detector.segments(j).y2       = parm_segment(j).seg_x1_x2_y1_y2(4);
    xparm.detector.segments(j).orgxs    = parm_segment(j).seg_orgx_orgy(1);
    xparm.detector.segments(j).orgys    = parm_segment(j).seg_orgx_orgy(2);
    xparm.detector.segments(j).fs       = parm_segment(j).seg_distance;
    xparm.detector.segments(j).eds(:,1) = parm_segment(j).direction_of_segment_x_axis;
    xparm.detector.segments(j).eds(:,2) = parm_segment(j).direction_of_segment_y_axis;
end
end