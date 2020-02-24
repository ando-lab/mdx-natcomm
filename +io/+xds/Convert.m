classdef Convert
    %CONVERT class for converting between xds input and parm and geom objects
    % (geom.Crystal, geom.Detector, geom.Source, geom.Spindle, io.ImageSeries, etc)
    methods(Static)
        function IS = integrate2ImageSeries(xdsWedges)
            IS = io.ImageSeries.empty();
            for j=1:length(xdsWedges)
                frames = [xdsWedges(j).table.image];
                frameRange = frames([1,end]);
                excludedFrames = setdiff(frameRange(1):frameRange(2),frames);
                IS(j) = io.ImageSeries('frameRange',frameRange,...
                    'excludedFrames',excludedFrames);
            end
        end
        function IS = inp2ImageSeries(xdsInp)
            NameTemplate = io.xds.FileNameTemplate(...
                xdsInp.name_template_of_data_frames);
            excludes = [];
            for j=1:size(xdsInp.exclude_data_range,1)
                thisRange = xdsInp.exclude_data_range(j,:);
                excludes = [excludes,thisRange(1):thisRange(2)];
            end
            excludes = unique(excludes);
            IS = io.ImageSeries('fileNameTemplate',NameTemplate.format,...
                'frameRange',xdsInp.data_range,'excludedFrames',excludes);
            
        end
        function DE = integrate2DiffractionExperiment(xdsWedges)
            DE = geom.DiffractionExperiment.empty();
            for j=1:length(xdsWedges)
                DE(j) = io.xds.Convert.parm2DiffractionExperiment(...
                    xdsWedges(j).xparm);
                DE(j).Spindle = geom.RotationSeries(DE(j).Spindle);
                DE(j).Spindle.seriesFrameRange = xdsWedges(j).dataRange;
            end
        end
        function DE = inp2DiffractionExperiment(xdsInp)
            DE = geom.DiffractionExperiment();
            DE.Crystal = io.xds.Convert.inp2Crystal(xdsInp);
            DE.Source = io.xds.Convert.inp2Source(xdsInp);
            DE.Detector = io.xds.Convert.inp2Detector(xdsInp);
            DE.Spindle = io.xds.Convert.inp2Spindle(xdsInp);
            % apply parallax corrections to xShift and yShift
            [dx,dy] = geom.Corrections.parallax(DE.Source,DE.Detector);
            % rounding dx and dy to two decimal places reduces file size
            dx = 0.01*round(dx*100);
            dy = 0.01*round(dy*100);
            DE.Detector.xShift = DE.Detector.xShift + dx;
            DE.Detector.yShift = DE.Detector.yShift + dy;
        end
        function DE = parm2DiffractionExperiment(xdsParm)
            DE = geom.DiffractionExperiment();
            DE.Crystal = io.xds.Convert.parm2Crystal(xdsParm);
            DE.Source = io.xds.Convert.parm2Source(xdsParm);
            DE.Detector = io.xds.Convert.parm2Detector(xdsParm);
            DE.Spindle = io.xds.Convert.parm2Spindle(xdsParm);
        end
        function Crystal = inp2Crystal(xdsInp)
            Crystal = geom.Crystal();
            if ~isempty(xdsInp.space_group_number)
                Crystal.spaceGroupNumber = xdsInp.space_group_number;
            end
            if ~isempty(xdsInp.unit_cell_constants)
                Crystal.a = xdsInp.unit_cell_constants(1);
                Crystal.b = xdsInp.unit_cell_constants(2);
                Crystal.c = xdsInp.unit_cell_constants(3);
                Crystal.alpha = xdsInp.unit_cell_constants(4);
                Crystal.beta = xdsInp.unit_cell_constants(5);
                Crystal.gamma = xdsInp.unit_cell_constants(6);
            end
            Crystal.a_axis = xdsInp.unit_cell_a_axis;
            Crystal.b_axis = xdsInp.unit_cell_b_axis;
            Crystal.c_axis = xdsInp.unit_cell_c_axis;
        end
        function Crystal = parm2Crystal(xdsParm)
            Crystal = geom.Crystal();
            Crystal.spaceGroupNumber = xdsParm.space_group_number;
            Crystal.a = xdsParm.a;
            Crystal.b = xdsParm.b;
            Crystal.c = xdsParm.c;
            Crystal.alpha = xdsParm.alpha;
            Crystal.beta = xdsParm.beta;
            Crystal.gamma = xdsParm.gamma;
            Crystal.a_axis = xdsParm.a_axis;
            Crystal.b_axis = xdsParm.b_axis;
            Crystal.c_axis = xdsParm.c_axis;
        end
        function Detector = inp2Detector(xdsInp)
            if strcmpi(xdsInp.detector,'pilatus') && ...
                    xdsInp.nx == 2463 && xdsInp.ny == 2527 && ...
                    xdsInp.qx == 0.172 && xdsInp.qy == 0.172
                Detector = geom.Pilatus6m();
            else
                warning('detector type %s in xdsInp not recognized',...
                    xdsInp.detector);
                Detector = geom.Detector('name',xdsInp.detector);
                Detector.nx = xdsInp.nx;
                Detector.ny = xdsInp.ny;
                Detector.qx = xdsInp.qx;
                Detector.qy = xdsInp.qy;
            end
            
            if ~isempty(xdsInp.x_geo_corr)
                try
                    xCorr = io.DectrisCBF(xdsInp.x_geo_corr);
                    Detector.xShift = 0.01*double(xCorr.image);
                catch
                    warning('could not load correction file\n:\t%s',...
                        xdsInp.x_geo_corr);
                    Detector.xShift = 0;
                end
            else
                Detector.xShift = 0;
            end
            
            if ~isempty(xdsInp.y_geo_corr)
                try
                    yCorr = io.DectrisCBF(xdsInp.y_geo_corr);
                    Detector.yShift = 0.01*double(yCorr.image);
                catch
                    warning('could not load correction file\n:\t%s',...
                        xdsInp.y_geo_corr);
                    Detector.yShift = 0;
                end
            else
                Detector.yShift = 0;
            end
            
            Detector.f = xdsInp.detector_distance;
            Detector.orgx = xdsInp.orgx;
            Detector.orgy = xdsInp.orgy;
            
            Detector.sensorThickness = xdsInp.sensor_thickness;
            Detector.mask = io.xds.Inp.getMask(xdsInp);
        end
        function Detector = parm2Detector(xdsParm)
            Detector = geom.Detector();
            Detector.nx = xdsParm.detector.nx;
            Detector.ny = xdsParm.detector.ny;
            Detector.qx = xdsParm.detector.qx;
            Detector.qy = xdsParm.detector.qy;
            Detector.ed = xdsParm.detector.ed;
            Detector.f = xdsParm.detector.f;
            Detector.orgx = xdsParm.detector.orgx;
            Detector.orgy = xdsParm.detector.orgy;
        end
        function Source = inp2Source(xdsInp)
            Source = geom.Source();
            Source.fractionOfPolarization = ...
                xdsInp.fraction_of_polarization;
            Source.polarizationPlaneNormal = ...
                xdsInp.polarization_plane_normal;
            Source.wavelength = xdsInp.x_ray_wavelength;
            Source.direction = xdsInp.incident_beam_direction;
        end
        function Source = parm2Source(xdsParm)
            Source = geom.Source();
            Source.wavelength = xdsParm.lambda;
            Source.direction = xdsParm.incident_wavevector*xdsParm.lambda;
        end
        function Spindle = inp2Spindle(xdsInp)
            Spindle = geom.Spindle();
            Spindle.rotationAxis=xdsInp.rotation_axis;
            Spindle.startingFrame=xdsInp.starting_frame;
            Spindle.startingAngle=xdsInp.starting_angle;
            Spindle.oscillationRange=xdsInp.oscillation_range;
        end
        function Spindle = parm2Spindle(xdsParm)
            Spindle = geom.Spindle();
            Spindle.rotationAxis=xdsParm.rotation_axis;
            Spindle.startingFrame=xdsParm.starting_frame;
            Spindle.startingAngle=xdsParm.starting_angle;
            Spindle.oscillationRange=xdsParm.oscillation_range;
        end
    end
end