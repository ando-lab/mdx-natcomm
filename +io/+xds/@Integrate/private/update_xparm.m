function xparm = update_xparm(parm,startingFrame,xparm_in)

xparm = xparm_in;

xparm.starting_frame = startingFrame;
xparm.starting_angle = (startingFrame - xparm_in.starting_frame)*...
    xparm_in.oscillation_range + xparm_in.starting_angle;
xparm.rotation_axis = parm.rotation_axis;
xparm.incident_wavevector = parm.direct_beam_coords;
xparm.space_group_number = parm.space_group_number;
xparm.a = parm.unit_cell(1);
xparm.b = parm.unit_cell(2);
xparm.c = parm.unit_cell(3);
xparm.alpha = parm.unit_cell(4);
xparm.beta = parm.unit_cell(5);
xparm.gamma = parm.unit_cell(6);
xparm.a_axis = parm.a_axis;
xparm.b_axis = parm.b_axis;
xparm.c_axis = parm.c_axis;
xparm.detector.orgx = parm.detector_origin_pixels(1);
xparm.detector.orgy = parm.detector_origin_pixels(2);
xparm.detector.f = parm.f;

end