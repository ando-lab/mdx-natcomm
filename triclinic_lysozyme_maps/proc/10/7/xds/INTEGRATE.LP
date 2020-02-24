
 ***** INTEGRATE ***** (VERSION Mar 15, 2019  BUILT=20190806)   6-Jan-2020

 INPUT PARAMETER VALUES
 ----------------------
 MAXIMUM_NUMBER_OF_PROCESSORS=   0  (0: automatic choice of OpenMP threads)
 number of OpenMP threads used  20
 MAXIMUM_NUMBER_OF_JOBS=   0  (0:automatic choice of forked main tasks)
 number of forked INTEGRATE tasks   1
 CLUSTER_NODES=
 SECONDS=     0
 NAME_TEMPLATE_OF_DATA_FRAMES=~/data/lysozyme/nitrate/ubatch4/lys_nitr_10_7_????.cbf CBF
 DATA_RANGE=       1     500
 NUMBER_OF_IMAGES_IN_CACHE=   21
 MINIMUM_VALID_PIXEL_VALUE=     0    OVERLOAD=   1048500
 BACKGROUND_PIXEL=    6.00      SIGNAL_PIXEL=    3.00
 MAXIMUM_ERROR_OF_SPOT_POSITION=   3.0
 MAXIMUM_ERROR_OF_SPINDLE_POSITION=   2.0
 MINPK=  75.00000
 CUT=  2.00
 PROFILE_FITTING= TRUE
 SPACE_GROUP_NUMBER=    1
 UNIT_CELL_CONSTANTS=    27.42    32.11    34.50  88.648 108.486 111.840
 ROTATION_AXIS= -0.99986  0.01639 -0.00268
 OSCILLATION_RANGE=  0.100000 DEGREES
 STARTING_ANGLE= 115.000 STARTING_FRAME=       1
 X-RAY_WAVELENGTH=   0.97700 ANGSTROM
 INCIDENT_BEAM_DIRECTION=  0.000629  0.005163  1.023528
 NUMBER OF DETECTOR SEGMENTS   1
 NX=  2463 NY=  2527 QX= 0.17200 QY= 0.17200
 ORGX=   1265.55  ORGY=   1237.34
 DETECTOR_DISTANCE=   231.324
 DIRECTION_OF_DETECTOR_X-AXIS=   1.00000   0.00000   0.00000
 DIRECTION_OF_DETECTOR_Y-AXIS=   0.00000   1.00000   0.00000
 SEGMENT=    1 2463    1 2527  ! x1,x2,y1,y2 (pixel units)
    SEGMENT_ORGX=      0.00   SEGMENT_ORGY=      0.00
    SEGMENT_DISTANCE=     0.00
    DIRECTION_OF_SEGMENT_X-AXIS=  1.000000  0.000000  0.000000
    DIRECTION_OF_SEGMENT_Y-AXIS=  0.000000  1.000000  0.000000
 RESOLUTION RANGE RECORDED BY DETECTOR (ANGSTROM)    44.486     1.212
 NUMBER OF TRUSTED DETECTOR PIXELS    5197940
 MEAN CONTENTS OF TRUSTED PIXELS IN BACKGROUND TABLE    136.346
 MEAN VALUE OF NON-XRAY BACKGROUND (OFFSET)=    0.00
 NUMBER OF X-RAY COUNTS EQUIVALENT TO PIXEL CONTENTS IN A DATA IMAGE
 (PIXEL VALUE - OFFSET)/GAIN ,    WITH GAIN=    1.00



 ******************************************************************************
                     PROCESSING OF IMAGES        1 ...      20
 ******************************************************************************

 USING  20 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ   SIGMAB   SIGMAR
     1   0  1.067  5187034    0    182     109     1  0.04199  0.13152
     2   0  1.041  5186827    0    168     115     1  0.04495  0.15134
     3   0  1.047  5187786    0    169     108     0  0.04504  0.15974
     4   0  1.045  5187980    0    156     106     1  0.04144  0.14623
     5   0  1.049  5187828    0    155     107     2  0.04292  0.13403
     6   0  1.048  5188191    0    159     102     1  0.04208  0.12811
     7   0  1.045  5188339    0    160      99     1  0.04264  0.14986
     8   0  1.052  5188386    0    151     103     0  0.03987  0.14762
     9   0  1.057  5188574    0    154      99     0  0.04235  0.14833
    10   0  1.055  5188193    0    152     102     0  0.04286  0.15484
    11   0  1.051  5188547    0    158     102     0  0.04519  0.15173
    12   0  1.056  5188552    0    159     100     4  0.04319  0.12664
    13   0  1.057  5188073    0    160      99     1  0.04598  0.14350
    14   0  1.062  5187919    0    162     105     1  0.04187  0.14415
    15   0  1.059  5186885    0    164     112     2  0.04276  0.12842
    16   0  1.059  5187119    0    169     114     1  0.04308  0.16078
    17   0  1.062  5187907    0    173      98     2  0.04043  0.12739
    18   0  1.061  5188013    0    160      97     0  0.04067  0.14546
    19   0  1.063  5189335    0    154      88     0  0.04364  0.15782
    20   0  1.064  5189922    0    146      83     0  0.04131  0.16486

   516 OUT OF    519 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.55
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.411    32.086    34.502  88.666 108.485 111.850
 REC. CELL PARAMETERS   0.041675  0.033775  0.030739  83.822  70.589  67.372
 COORDINATES OF UNIT CELL A-AXIS   -15.212   -14.287    17.771
 COORDINATES OF UNIT CELL B-AXIS    12.090    26.579    13.299
 COORDINATES OF UNIT CELL C-AXIS   -20.971    19.920   -18.810
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.068    -0.082     0.068
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.045
 LAB COORDINATES OF ROTATION AXIS -0.999862  0.016385 -0.002675
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000692  0.005073  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.32   1244.14
 DETECTOR ORIGIN (PIXELS) AT                     1265.41   1237.47
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.28
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     516 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   54     36     65     71     56     31     54     75     62
 SIGMAB (degree)       0.035  0.037  0.039  0.038  0.035  0.035  0.034  0.033  0.033
 SIGMAR (degree)       0.037  0.035  0.037  0.036  0.035  0.035  0.035  0.036  0.033


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    516 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  2  3  4  3  2  0  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  1  1  2  3  2  1  0  0  0  0     0  0  0  1  4 11 16 11  6  2  1  0  0
     0  0  0  0  1  0  1  0  0  0  0  0  0     0  0  0  0  1  3  4  2  1  1  0  0  0     0  0  0  1  7 18 24 17  9  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  1  4 13 18 13  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  1  4  7  6  4  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  2  4  4  3  1  0  0  0  0     0  0  0  0  2  4  4  2  1  0  0  0  0
     0  0  0  1  5 11 14 10  5  2  0  0  0     0  0  0  2  9 19 22 14  7  3  1  0  0     0  0  0  2  8 16 16  9  4  2  0  0  0
     0  0  0  3 16 39 52 35 17  7  2  1  0     0  0  0  5 26 60 72 45 21  9  3  1  0     0  0  0  4 18 38 39 23 10  4  1  0  0
     0  0  0  5 25 62 76 53 25 11  4  1  0     0  0  1  7 40 89100 64 29 13  4  1  0     0  0  0  5 24 48 49 30 13  6  2  0  0
     0  0  0  2 14 43 56 39 20  9  3  1  0     0  0  0  4 24 60 69 44 21  9  3  1  0     0  0  0  3 16 32 34 20  9  3  1  0  0
     0  0  0  1  5 16 21 18  9  4  1  0  0     0  0  0  1  9 22 25 17  8  3  1  0  0     0  0  0  1  5 12 12  7  3  1  0  0  0
     0  0  0  0  1  4  6  5  2  1  0  0  0     0  0  0  0  2  4  6  5  2  1  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  3  2  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6  8  7  4  1  1  0  0  0     0  0  0  1  2  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  8 13 11  5  2  1  0  0  0     0  0  0  1  2  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  7 11 10  6  2  1  0  0  0     0  0  0  1  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  6  3  1  0  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    686 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES       21 ...      40
 ******************************************************************************

 USING  20 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ   SIGMAB   SIGMAR
    21   0  1.067  5189628    0    140      90     7  0.04119  0.09288
    22   0  1.068  5189117    0    143      94     5  0.04332  0.10264
    23   0  1.066  5188766    0    156      97     3  0.03969  0.11034
    24   0  1.069  5188264    0    157     101     4  0.04144  0.09801
    25   0  1.071  5187710    0    168     104     3  0.04201  0.09509
    26   0  1.070  5188154    0    174     103     4  0.04010  0.08056
    27   0  1.074  5187921    0    172     102     4  0.03967  0.08288
    28   0  1.073  5188113    0    167     104     4  0.04037  0.09169
    29   0  1.075  5187927    0    170     106     1  0.04158  0.11698
    30   0  1.077  5188128    0    158     109     1  0.04401  0.13256
    31   0  1.077  5188547    0    154     106     1  0.04347  0.15902
    32   0  1.080  5189204    0    147      95     0  0.04454  0.15179
    33   0  1.076  5189227    0    143      95     1  0.04353  0.12199
    34   0  1.081  5190158    0    139      84     1  0.04286  0.11804
    35   0  1.082  5189909    0    136      91     0  0.04476  0.12815
    36   0  1.082  5189878    0    132      79     1  0.04103  0.08706
    37   0  1.084  5189525    0    140      83     0  0.04234  0.09526
    38   0  1.086  5189466    0    153      89     1  0.04196  0.10029
    39   0  1.085  5189559    0    143      87     1  0.04022  0.08974
    40   0  1.088  5189792    0    145      88     3  0.03982  0.09432

   479 OUT OF    481 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.56
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.13
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.404    32.088    34.502  88.638 108.480 111.837
 REC. CELL PARAMETERS   0.041676  0.033767  0.030736  83.862  70.607  67.396
 COORDINATES OF UNIT CELL A-AXIS   -15.216   -14.273    17.769
 COORDINATES OF UNIT CELL B-AXIS    12.086    26.584    13.297
 COORDINATES OF UNIT CELL C-AXIS   -20.963    19.928   -18.809
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.068    -0.064     0.057
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.045
 LAB COORDINATES OF ROTATION AXIS -0.999862  0.016385 -0.002675
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000495  0.005138  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.25   1244.11
 DETECTOR ORIGIN (PIXELS) AT                     1265.60   1237.36
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.24
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     479 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   48     33     60     68     54     31     55     65     53
 SIGMAB (degree)       0.034  0.036  0.039  0.037  0.035  0.036  0.034  0.032  0.032
 SIGMAR (degree)       0.036  0.036  0.038  0.035  0.033  0.035  0.035  0.034  0.032


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    479 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     1  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0 -1     0  0  0  1  3  8 10  7  3  1  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  1  5 12 17 12  6  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  1  4  9 13 10  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  4  6  4  3  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  2  4  4  3  1  0  0 -1  0     0  0  0  0  2  4  4  3  1  0  0  0  0
     0  0  0  1  4  9 10  7  3  1  0  0  0     0  0  0  1  9 19 21 14  6  3  1  0  0     0  0  0  2  8 16 16 10  5  2  0  0  0
     0  0  0  2 12 32 39 27 13  5  2  0  0     0  0  0  5 25 61 68 44 20  8  3  0  0     0  0  0  4 20 45 46 27 12  5  2  0  0
     0  0  0  3 19 51 64 45 21  9  3  1  0     0  0  0  6 36 87100 64 30 12  4  1  0     0  0  0  5 27 58 62 36 17  7  2  0  0
     0  0  0  2 14 36 49 36 18  7  3  1  0     0  0  0  4 24 59 71 47 22  9  3  1  0     0  0  0  3 18 40 43 26 11  4  1  0  0
     0  0  0  0  4 14 21 16  9  4  1  0  0     0  0  0  1  8 22 28 18 10  4  1  0  0     0  0  0  1  6 14 15  8  4  1  0  0  0
     0  0  0  0  1  3  6  5  3  1  0  0  0     0  0  0  0  2  4  7  5  2  1  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0 -1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  8  6  4  1  0  0  0  0     0  0  0  1  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  9 14 13  7  3  1  0  0  0     0  0  0  1  2  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3 10 16 15  8  4  1  0  0  0     0  0  0  1  1  2  1  1  0  0  0  0 -1     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  6 11 10  6  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  3  3  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    637 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES       41 ...      60
 ******************************************************************************

 USING  20 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ   SIGMAB   SIGMAR
    41   0  1.091  5190024    0    140      89     1  0.04278  0.11690
    42   0  1.091  5189812    0    147      87     1  0.04434  0.09956
    43   0  1.096  5189249    0    151      96     2  0.04004  0.10290
    44   0  1.096  5187963    0    155     104     3  0.04092  0.09163
    45   0  1.097  5188048    0    155     100     4  0.04040  0.08670
    46   0  1.097  5188774    0    160      94     2  0.04419  0.09669
    47   0  1.099  5188413    0    154     102     1  0.04312  0.11056
    48   0  1.103  5188687    0    152     102     1  0.04191  0.12805
    49   0  1.105  5189606    0    141      91     1  0.04354  0.13010
    50   0  1.106  5189612    0    151      87     0  0.04546  0.14316
    51   0  1.108  5189237    0    151      92     5  0.04455  0.09516
    52   0  1.105  5188602    0    147      99     1  0.04114  0.10164
    53   0  1.108  5189091    0    152      95     2  0.03906  0.08701
    54   0  1.112  5189361    0    141      90     4  0.04126  0.07975
    55   0  1.111  5189915    0    133      86     6  0.04445  0.08660
    56   0  1.113  5189866    0    126      87     0  0.04243  0.12534
    57   0  1.114  5190374    0    124      77     0  0.04453  0.12121
    58   0  1.118  5190401    0    140      78     1  0.04347  0.10519
    59   0  1.114  5189483    0    149      92     1  0.04198  0.10216
    60   0  1.113  5189203    0    151      96     1  0.04016  0.09322

   467 OUT OF    468 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.55
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.404    32.095    34.498  88.640 108.487 111.846
 REC. CELL PARAMETERS   0.041682  0.033762  0.030741  83.853  70.598  67.386
 COORDINATES OF UNIT CELL A-AXIS   -15.219   -14.276    17.763
 COORDINATES OF UNIT CELL B-AXIS    12.088    26.590    13.302
 COORDINATES OF UNIT CELL C-AXIS   -20.954    19.927   -18.814
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.068    -0.055     0.062
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.045
 LAB COORDINATES OF ROTATION AXIS -0.999862  0.016385 -0.002675
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000615  0.005149  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.27   1244.07
 DETECTOR ORIGIN (PIXELS) AT                     1265.47   1237.31
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.23
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     467 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   48     29     56     67     55     30     54     60     57
 SIGMAB (degree)       0.034  0.036  0.038  0.038  0.036  0.035  0.034  0.032  0.033
 SIGMAR (degree)       0.036  0.034  0.036  0.035  0.033  0.036  0.036  0.035  0.033


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    467 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  1  3  8 10  6  3  1  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  2  2  2  1  0  0  0  0     0  0  0  1  4 11 15 10  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  0  3  7 10  8  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  3  4  3  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  3  4  2  1  0  0  0  0     0  0  0  0  2  4  4  3  1  0  0  0  0
     0  0  0  1  4 10 12  7  3  1  0  0  0     0  0  0  2  8 20 22 13  6  2  1  0  0     0  0  0  2  9 19 19 11  5  2  1  0  0
     0  0  0  2 13 32 40 25 12  5  2  0  0     0  0  1  5 25 61 71 42 21  8  3  1  0     0  0  1  5 21 47 48 27 13  5  2  0  0
     0  0  0  4 18 46 57 39 19  8  3  1  0     0  0  1  8 38 87100 66 30 12  4  1  0     0  0  1  7 32 63 66 41 18  7  2  0  0
     0  0  0  2 11 31 41 29 15  6  2  0  0     0  0  0  4 25 60 73 48 23  9  3  1  0     0  0  0  4 22 44 48 29 13  5  2  0  0
     0  0  0  0  3 12 15 12  6  3  1  0  0     0  0  0  1  7 21 26 20  9  4  1  0  0     0  0  0  1  6 15 16 11  5  2  0  0  0
     0  0  0  0  1  2  4  3  2  1  0  0  0     0  0  0  0  2  4  5  4  2  1  0  0  0     0  0  0  0  1  3  3  2  1  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0 -1  0     0  0  0  0  0  1  1  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  3  3  2  1  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  8  8  4  2  1  0  0  0     0  0  0  1  2  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  9 15 15  8  4  1  0  0  0     0  0  0  1  2  3  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3 11 18 16 10  4  2  1  0  0     0  0  0  1  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 11 11  6  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  3  3  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    640 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES       61 ...      80
 ******************************************************************************

 USING  20 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ   SIGMAB   SIGMAR
    61   0  1.118  5187818    0    156     112     4  0.04116  0.09288
    62   0  1.120  5188349    0    149     106     4  0.04189  0.10773
    63   0  1.119  5189150    0    153      90     7  0.04344  0.08110
    64   0  1.120  5188358    0    152     100     8  0.04331  0.07845
    65   0  1.122  5188291    0    155     104    10  0.04099  0.08129
    66   0  1.119  5188217    0    170     108     7  0.04254  0.08660
    67   0  1.126  5188444    0    167     111     5  0.03984  0.08901
    68   0  1.125  5188449    0    159     109     2  0.04129  0.11004
    69   0  1.126  5190025    0    156      91     0  0.04326  0.12270
    70   0  1.127  5190394    0    150      87     1  0.04076  0.12196
    71   0  1.128  5189784    0    137      87     1  0.04296  0.10823
    72   0  1.125  5189710    0    136      93     1  0.04575  0.11213
    73   0  1.126  5189795    0    141      90     1  0.04230  0.11431
    74   0  1.125  5190151    0    149      84     0  0.04355  0.12914
    75   0  1.133  5190063    0    142      87     0  0.04290  0.12959
    76   0  1.131  5189293    0    148      91     1  0.04157  0.11126
    77   0  1.134  5189090    0    140     100     1  0.04129  0.11340
    78   0  1.135  5189852    0    140      86     0  0.04087  0.12265
    79   0  1.134  5190077    0    129      87     0  0.04298  0.12425
    80   0  1.130  5190279    0    135      85     0  0.04086  0.13095

   489 OUT OF    490 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.56
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.13
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.407    32.088    34.493  88.674 108.493 111.869
 REC. CELL PARAMETERS   0.041690  0.033778  0.030750  83.803  70.574  67.348
 COORDINATES OF UNIT CELL A-AXIS   -15.225   -14.279    17.761
 COORDINATES OF UNIT CELL B-AXIS    12.092    26.584    13.294
 COORDINATES OF UNIT CELL C-AXIS   -20.946    19.909   -18.832
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.055    -0.042     0.065
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.045
 LAB COORDINATES OF ROTATION AXIS -0.999862  0.016385 -0.002675
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.001079  0.005148  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.32   1244.08
 DETECTOR ORIGIN (PIXELS) AT                     1264.90   1237.32
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.22
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     489 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   51     30     61     66     52     31     60     70     55
 SIGMAB (degree)       0.034  0.036  0.038  0.037  0.035  0.035  0.033  0.032  0.032
 SIGMAR (degree)       0.038  0.034  0.037  0.037  0.034  0.035  0.035  0.035  0.033


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    489 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0 -1  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  3  5  3  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     1  0  0  1  4 12 15 11  5  2  1  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  3  4  3  1  0  0  0  0     0  0  0  1  6 16 22 16  7  3  1  0  0
     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  0  4 12 16 11  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  4  6  4  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  1  3  3  2  0  0  0  0  0     0  0  0  0  1  3  3  2  1  0  0  0  0
     0  0  0  1  4 11 14  8  4  1  0  0  0     0  0  0  1  7 18 20 11  5  2  0  0  0     0  0  0  1  6 14 14  8  4  1  0  0  0
     0  0  0  2 14 37 46 31 15  6  2  0 -1     0  0  0  4 24 58 67 42 20  8  3  0  0     0  0  0  4 19 38 38 22 11  4  1  0  0
     0  0  0  4 21 54 69 49 23  9  3  1  0     0  0  0  6 36 87100 68 30 13  4  1  0     0  0  0  5 25 53 53 34 15  6  2  0  0
     0  0  0  2 15 40 51 35 16  7  2  0  0     0  0  0  3 24 62 71 47 22  9  3  1  0     0  0  0  3 16 36 36 22 10  4  1  0  0
     0  0  0  1  5 13 18 14  7  3  1  0  0     0  0  0  1  8 19 25 18  9  4  1  0  0     0  0  0  1  5 10 12  8  4  2  0  0  0
     0  0  0  0  1  3  4  2  1  0  0  0  0     0  0  0  0  2  4  5  4  2  1  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  6  4  2  0  0  0  0     0  0  0  1  2  3  2  1  1  0  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0
     0  0  0  2  8 13 11  6  3  1  0  0  0     0  0  0  1  3  3  3  2  1  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  2  7 14 11  7  3  1  0  0  0     0  0  0  1  2  3  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  1  4  8  6  4  2  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    653 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES       81 ...     100
 ******************************************************************************

 USING  20 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ   SIGMAB   SIGMAR
    81   0  1.131  5189351    0    138      92     0  0.04086  0.12300
    82   0  1.139  5188649    0    144      99     5  0.04270  0.10477
    83   0  1.134  5188793    0    156      98     5  0.04114  0.09425
    84   0  1.140  5188422    0    159     104     6  0.04066  0.09121
    85   0  1.145  5187818    0    159     109     6  0.03928  0.09072
    86   0  1.146  5188594    0    157     101     3  0.04087  0.10443
    87   0  1.144  5189639    0    154      89     3  0.04396  0.10172
    88   0  1.146  5189814    0    152      86     2  0.04057  0.09308
    89   0  1.146  5189023    0    153      92     1  0.03960  0.10459
    90   0  1.150  5187920    0    153     105     4  0.04001  0.09277
    91   0  1.152  5188540    0    160      99     1  0.04220  0.09695
    92   0  1.152  5189349    0    153      89     3  0.04218  0.11288
    93   0  1.154  5189475    0    140      87     4  0.03996  0.09702
    94   0  1.156  5189803    0    138      83     5  0.04336  0.09408
    95   0  1.155  5189907    0    144      83     2  0.04091  0.10436
    96   0  1.154  5189386    0    142      91     1  0.04250  0.10227
    97   0  1.159  5189152    0    140      99     6  0.03762  0.09035
    98   0  1.155  5189432    0    140      94     3  0.03935  0.09069
    99   0  1.157  5190003    0    148      86     3  0.03968  0.09368
   100   0  1.162  5189476    0    138      90     4  0.04313  0.10014

   511 OUT OF    513 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.58
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.15
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.407    32.094    34.495  88.671 108.489 111.872
 REC. CELL PARAMETERS   0.041689  0.033772  0.030747  83.807  70.579  67.347
 COORDINATES OF UNIT CELL A-AXIS   -15.230   -14.277    17.759
 COORDINATES OF UNIT CELL B-AXIS    12.092    26.590    13.295
 COORDINATES OF UNIT CELL C-AXIS   -20.944    19.909   -18.838
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.049    -0.031     0.065
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.045
 LAB COORDINATES OF ROTATION AXIS -0.999862  0.016385 -0.002675
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000991  0.005135  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.38   1244.05
 DETECTOR ORIGIN (PIXELS) AT                     1265.08   1237.31
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.23
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     511 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   53     33     59     69     53     30     57     79     64
 SIGMAB (degree)       0.034  0.036  0.037  0.037  0.035  0.035  0.033  0.032  0.033
 SIGMAR (degree)       0.037  0.037  0.036  0.035  0.034  0.034  0.036  0.036  0.034


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    511 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  3  4  3  1  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  2  1  1  1  0  0  0     0  0  0  1  4 11 14 10  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  3  1  1  0  0  0     0  0  0  1  5 16 21 14  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  1  3 10 14 11  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  4  5  4  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  1  3  4  2  1  0  0  0  0     0  0  0  0  2  3  3  2  1  0  0  0  0
     0  0  0  1  4 10 12  9  4  1  0  0  0     0  0  0  1  7 18 20 13  5  2  1  0  0     0  0  0  1  6 14 14  9  3  1  0  0  0
     0  0  0  2 14 37 45 31 15  6  2  0  0     0  0  0  4 26 62 66 44 20  8  3  1  0     0  0  1  4 19 42 40 25 11  4  1  0  0
     0  0  0  4 19 55 67 46 22  9  3  1  0     0  0  0  7 34 90100 65 30 12  4  1  0     0  0  0  6 24 57 57 34 16  6  2  0  0
     0  0  0  2 11 35 46 33 17  7  2  1  0     0  0  0  4 21 58 69 46 22  9  3  1  0     0  0  0  3 15 36 39 24 11  5  1  0  0
     0  0  0  1  4 12 15 12  7  3  1  0  0     0  0  0  1  7 19 23 16  9  3  1  0  0     0  0  0  1  5 12 13  8  4  1  0  0  0
     0  0  0  0  1  3  3  3  2  1  0  0  0     0  0  0  0  1  3  5  4  2  1  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  6  6  3  1  0  0  0  0     0  0  0  1  2  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 13 11  6  3  1  0  0  0     0  0  0  1  2  3  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 14 13  7  3  1  0  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  8  8  5  2  1  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  3  3  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    644 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      101 ...     120
 ******************************************************************************

 USING  20 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ   SIGMAB   SIGMAR
   101   0  1.162  5188997    0    145      93     0  0.04219  0.10988
   102   0  1.160  5189330    0    155      90     1  0.04513  0.11794
   103   0  1.165  5189414    0    159      88     1  0.04317  0.11115
   104   0  1.167  5189129    0    144      93     3  0.03967  0.08722
   105   0  1.168  5190114    0    145      81     3  0.03921  0.08254
   106   0  1.170  5190987    0    149      74     4  0.04032  0.08397
   107   0  1.170  5190430    0    141      79     2  0.04194  0.09880
   108   0  1.170  5190598    0    142      82     3  0.04116  0.09862
   109   0  1.171  5191070    0    137      76     3  0.04339  0.11567
   110   0  1.174  5191498    0    136      72     0  0.04258  0.14541
   111   0  1.170  5191143    0    140      79     3  0.04306  0.12421
   112   0  1.169  5189841    0    147      92     5  0.04058  0.10206
   113   0  1.171  5189237    0    156     103     3  0.03974  0.09462
   114   0  1.172  5189175    0    162      98     5  0.03881  0.07478
   115   0  1.174  5189586    0    158      99     6  0.03728  0.08806
   116   0  1.172  5190659    0    142      84     1  0.03714  0.10933
   117   0  1.174  5191428    0    137      76     3  0.04145  0.11549
   118   0  1.175  5191299    0    139      78     0  0.03870  0.11567
   119   0  1.177  5190460    0    151      82     0  0.03785  0.09672
   120   0  1.179  5189297    0    161      98     3  0.03974  0.09582

   468 OUT OF    472 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.56
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.398    32.095    34.488  88.629 108.484 111.833
 REC. CELL PARAMETERS   0.041684  0.033758  0.030749  83.872  70.607  67.403
 COORDINATES OF UNIT CELL A-AXIS   -15.229   -14.254    17.764
 COORDINATES OF UNIT CELL B-AXIS    12.087    26.595    13.293
 COORDINATES OF UNIT CELL C-AXIS   -20.933    19.919   -18.827
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.056    -0.017     0.044
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.045
 LAB COORDINATES OF ROTATION AXIS -0.999862  0.016385 -0.002675
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000397  0.005196  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.40   1244.07
 DETECTOR ORIGIN (PIXELS) AT                     1265.87   1237.25
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.17
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     468 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   46     33     63     67     52     26     56     64     51
 SIGMAB (degree)       0.034  0.036  0.037  0.036  0.035  0.035  0.034  0.033  0.032
 SIGMAR (degree)       0.038  0.034  0.035  0.034  0.034  0.037  0.038  0.036  0.033


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    468 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  1  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  1  0  0  0     0  0  0  1  2  6  7  5  2  1  0  0  0
     0  0  0  0  0  0  0  0  1  0  0  0 -1     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  0  3  9 11  9  5  2  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  1  0  0  0  0     0  0  0  0  2  6  8  6  3  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  3  3  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  1  1  1  0  0  0 -1  0     0  0  0  0  1  3  3  2  1  0  0 -1  0     0  0  0  0  2  4  4  2  1  0  0  0  0
     0  0  0  0  3  6  8  5  2  1  0  0  0     0  0  0  2  8 19 20 13  5  2  0  0  0     0  0  0  2 10 21 20 12  5  2  0  0  0
     0  0  0  2  9 26 32 22 10  4  1  0  0     0  0  0  5 24 61 68 44 20  8  3  0  0     0  0  1  6 26 57 57 34 15  7  2  0  0
     0  0  0  2 15 42 50 36 17  7  2  0  0     0  0  1  5 37 89100 67 31 13  4  1  0     0  0  1  6 34 77 82 51 23  9  3  0  0
     0  0  0  1 11 31 38 26 14  5  2  0  0     0  0  0  3 26 66 74 48 24  9  3  1  0     0  0  0  4 24 56 58 34 17  6  2  0  0
     0  0  0  0  3  9 13 11  6  3  0  0  0     0  0  0  1  8 19 25 18  9  4  1  0  0     0  0  0  1  7 16 18 12  5  2  0  0  0
     0  0  0  0  1  2  3  3  1  0  0  0  0     0  0  0  0  2  4  5  4  2  0  0  0  0     0  0  0  0  1  3  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0 -1 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  3  3  1  0  0  0  0  0    -1  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 11  9  5  2  1  0  0  0     0  0  0  1  2  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3 13 25 22 12  5  2  0  0  0     0  0  0  1  3  4  3  1  1  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  3 15 31 30 18  8  3  1  0  0     0  0  0  1  2  4  4  2  1  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  2  9 22 21 12  5  2  0  0  0     0  0  0  0  1  3  2  1  1  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  3  6  6  4  1  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    651 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      121 ...     140
 ******************************************************************************

 USING  20 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ   SIGMAB   SIGMAR
   121   0  1.176  5188930    0    164     105     2  0.04054  0.10842
   122   0  1.174  5189158    0    169     105     3  0.03909  0.09417
   123   0  1.178  5189554    0    172      96     4  0.04113  0.09073
   124   0  1.179  5188665    0    167     101     5  0.04105  0.08266
   125   0  1.183  5188020    0    169     109     6  0.04086  0.07540
   126   0  1.182  5188444    0    154     109     2  0.03996  0.09989
   127   0  1.182  5189866    0    154      95     1  0.04143  0.12694
   128   0  1.182  5189651    0    142      90     2  0.04349  0.12291
   129   0  1.182  5188960    0    134      92     3  0.04167  0.10066
   130   0  1.187  5189560    0    155      88     3  0.04105  0.10030
   131   0  1.184  5189155    0    149      96     4  0.04005  0.09538
   132   0  1.184  5188817    0    138      98     3  0.04196  0.09022
   133   0  1.184  5190258    0    135      80     1  0.04515  0.12052
   134   0  1.187  5190586    0    130      76     0  0.04233  0.11182
   135   0  1.191  5190611    0    137      78     0  0.04260  0.11147
   136   0  1.188  5189713    0    150      90     0  0.04104  0.10096
   137   0  1.191  5188939    0    144     101     0  0.04276  0.12326
   138   0  1.187  5190412    0    154      82     0  0.04216  0.13572
   139   0  1.188  5190039    0    146      90     2  0.04680  0.12752
   140   0  1.190  5189475    0    140      90     3  0.04001  0.11912

   473 OUT OF    477 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.56
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.404    32.096    34.494  88.630 108.490 111.839
 REC. CELL PARAMETERS   0.041678  0.033759  0.030744  83.866  70.600  67.395
 COORDINATES OF UNIT CELL A-AXIS   -15.235   -14.252    17.770
 COORDINATES OF UNIT CELL B-AXIS    12.089    26.599    13.285
 COORDINATES OF UNIT CELL C-AXIS   -20.932    19.919   -18.841
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.047    -0.006     0.036
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.045
 LAB COORDINATES OF ROTATION AXIS -0.999862  0.016385 -0.002675
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000644  0.005096  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.29   1244.09
 DETECTOR ORIGIN (PIXELS) AT                     1265.45   1237.39
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.23
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     473 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   46     32     63     67     49     30     54     66     56
 SIGMAB (degree)       0.034  0.036  0.037  0.036  0.035  0.036  0.034  0.033  0.032
 SIGMAR (degree)       0.036  0.032  0.035  0.034  0.036  0.039  0.037  0.035  0.032


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    473 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  2  7  9  6  3  1  1  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  1  4 12 13 10  5  2  1  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  1  1  2  1  0  0  0  0     0  0  0  0  3  8 11  9  4  2  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  3  5  3  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  2  1  0  0  0 -1  0     0  0  0  0  2  3  4  2  1  0  0 -1  0     0  0  0  0  2  4  3  2  1  0  0  0  0
     0  0  0  1  3  9 10  7  3  1  0  0  0     0  0  0  2  8 18 20 12  5  2  0  0  0     0  0  0  2  9 17 18 10  4  2  0  0  0
     0  0  0  2 10 31 36 25 11  4  2  0  0     0  0  1  5 26 62 67 42 18  8  3  0  0     0  0  0  5 24 51 51 29 13  5  2  0  0
     0  0  0  3 18 47 57 40 20  8  3  1  0     0  0  1  7 39 89100 63 30 13  4  1  0     0  0  1  6 32 69 71 41 19  8  2  0  0
     0  0  0  2 12 33 44 31 16  6  2  0  0     0  0  0  4 25 60 73 47 23  9  3  1  0     0  0  0  4 21 44 49 28 14  5  2  0  0
     0  0  0  1  4 11 16 11  6  2  1  0  0     0  0  0  1  7 18 23 16  8  3  1  0  0     0  0  0  1  6 13 15  9  4  1  0  0  0
     0  0  0  0  1  2  4  2  1  0  0  0  0     0  0  0  0  1  3  5  3  1  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  6  3  1  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  9 17 16  9  4  1  0  0  0     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3 13 23 21 12  5  2  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  8 14 14  7  3  1  0  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  4  4  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    657 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      141 ...     160
 ******************************************************************************

 USING  20 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ   SIGMAB   SIGMAR
   141   0  1.190  5189019    0    138      98     4  0.04360  0.10317
   142   0  1.189  5188836    0    165      95     1  0.04369  0.09911
   143   0  1.189  5188394    0    158     101     2  0.04114  0.09748
   144   0  1.190  5187493    0    164     109     0  0.04088  0.10399
   145   0  1.195  5187967    0    167     106     2  0.04318  0.10944
   146   0  1.195  5188768    0    157      97     3  0.04218  0.10620
   147   0  1.199  5189154    0    148      94     4  0.04337  0.11186
   148   0  1.200  5188935    0    141      99     0  0.04319  0.13551
   149   0  1.200  5189096    0    149      99     0  0.04491  0.13561
   150   0  1.203  5188577    0    149     103     1  0.04382  0.11875
   151   0  1.207  5189006    0    148      94     1  0.04216  0.11614
   152   0  1.203  5189121    0    141      97     2  0.04161  0.09886
   153   0  1.208  5189640    0    150      95     0  0.04124  0.11455
   154   0  1.208  5189471    0    165      92     0  0.04421  0.11270
   155   0  1.203  5188380    0    171     104     0  0.04334  0.11003
   156   0  1.205  5188075    0    172     112     1  0.04086  0.10134
   157   0  1.203  5188541    0    185     105     2  0.03965  0.08945
   158   0  1.207  5187962    0    174     110     2  0.03781  0.08710
   159   0  1.209  5187929    0    171     110     4  0.03755  0.07947
   160   0  1.211  5188474    0    175     109     2  0.03965  0.10354

   508 OUT OF    513 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.52
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.08
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.409    32.079    34.498  88.630 108.492 111.853
 REC. CELL PARAMETERS   0.041677  0.033781  0.030742  83.860  70.596  67.381
 COORDINATES OF UNIT CELL A-AXIS   -15.240   -14.249    17.774
 COORDINATES OF UNIT CELL B-AXIS    12.096    26.583    13.270
 COORDINATES OF UNIT CELL C-AXIS   -20.929    19.925   -18.844
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.056     0.002     0.030
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.045
 LAB COORDINATES OF ROTATION AXIS -0.999862  0.016385 -0.002675
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000540  0.005096  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.25   1244.11
 DETECTOR ORIGIN (PIXELS) AT                     1265.54   1237.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.24
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     508 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   57     31     58     72     59     34     57     70     60
 SIGMAB (degree)       0.034  0.035  0.037  0.036  0.035  0.035  0.033  0.032  0.032
 SIGMAR (degree)       0.036  0.032  0.035  0.035  0.035  0.036  0.036  0.035  0.032


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    508 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  3  4  2  1  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  0  4 10 13  9  4  2  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  3  3  3  1  0  0  0  1     0  0  0  1  6 16 21 16  7  3  1  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  1  4 12 16 11  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  4  5  4  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  2  3  3  2  1  0  0  0  0     0  0  0  0  2  3  2  1  0  0  0  0  0
     0  0  0  1  4 11 13  8  4  1  0  0  0     0  0  0  2  8 17 19 12  5  2  0  0  0     0  0  0  2  7 12 12  6  3  1  0  0  0
     0  0  0  2 15 38 47 31 14  6  2  0  0     0  0  0  5 27 61 66 40 18  7  2  0  0     0  0  0  4 19 38 37 21  9  3  1  0  0
     0  0  0  3 21 58 73 51 23  9  4  1  0     0  0  1  6 38 90100 64 28 11  5  1  0     0  0  1  5 27 55 54 32 14  6  2  0  0
     0  0  0  2 14 41 52 36 16  7  2  0  0     0  0  0  4 24 61 71 46 20  9  3  0  0     0  0  0  3 16 36 39 23 10  4  1  0  0
     0  0  0  1  4 13 17 12  6  3  1  0  0     0  0  0  1  7 19 24 15  7  3  1  0  0     0  0  0  1  5 12 13  8  3  1  0  0  0
     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  0  1  3  4  3  1  0  0  0  0     0  0  0  0  1  1  2  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  4  3  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  6  9  8  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 13 11  6  2  1  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  7  4  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    681 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      161 ...     180
 ******************************************************************************

 USING  20 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ   SIGMAB   SIGMAR
   161   0  1.211  5188710    0    184     110     2  0.04223  0.12320
   162   0  1.212  5188374    0    178     115     2  0.04145  0.10173
   163   0  1.215  5188214    0    160     111     2  0.04075  0.09426
   164   0  1.214  5189066    0    157      99     2  0.04163  0.10309
   165   0  1.214  5189296    0    151      98     1  0.03877  0.10069
   166   0  1.219  5190201    0    145      84     4  0.04378  0.08640
   167   0  1.219  5189923    0    151      90     1  0.04567  0.12854
   168   0  1.219  5189054    0    155     104     2  0.04330  0.12054
   169   0  1.220  5189320    0    156      96     3  0.04275  0.09301
   170   0  1.211  5189994    0    151      90     2  0.04028  0.11200
   171   0  1.220  5190801    0    146      83     2  0.03992  0.10726
   172   0  1.220  5190510    0    144      84     2  0.04156  0.11040
   173   0  1.216  5189575    0    145      94     2  0.04132  0.11640
   174   0  1.219  5189534    0    148      95     6  0.04133  0.08312
   175   0  1.219  5189650    0    152      88     6  0.04027  0.07360
   176   0  1.221  5190050    0    149      84     8  0.03933  0.06615
   177   0  1.220  5189720    0    147      95     3  0.04056  0.10281
   178   0  1.227  5190249    0    153      80     1  0.04259  0.10151
   179   0  1.226  5189923    0    149      84     2  0.04137  0.09654
   180   0  1.227  5189927    0    153      90     2  0.04232  0.10585

   475 OUT OF    481 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.51
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.15
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.406    32.085    34.496  88.631 108.488 111.859
 REC. CELL PARAMETERS   0.041682  0.033775  0.030743  83.859  70.599  67.375
 COORDINATES OF UNIT CELL A-AXIS   -15.241   -14.250    17.768
 COORDINATES OF UNIT CELL B-AXIS    12.093    26.591    13.273
 COORDINATES OF UNIT CELL C-AXIS   -20.928    19.920   -18.846
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.048     0.006     0.034
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.045
 LAB COORDINATES OF ROTATION AXIS -0.999862  0.016385 -0.002675
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000657  0.005081  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.28   1244.07
 DETECTOR ORIGIN (PIXELS) AT                     1265.42   1237.40
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.22
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     475 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   46     35     61     69     60     33     55     56     50
 SIGMAB (degree)       0.035  0.036  0.036  0.036  0.035  0.035  0.033  0.032  0.033
 SIGMAR (degree)       0.036  0.031  0.034  0.035  0.035  0.037  0.037  0.035  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    475 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  1  3  9 12  8  4  1  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  1  5 14 18 13  6  2  1  0  0
     0  0  0  0  0  1  0  1  0  1  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  1  3  9 12  9  4  2  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  1  3  5  3  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  0  0  0  0  0 -1  0  0  0  0     0 -1  0  0  0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0 -1  0  0  0  0  0  0     0  0  0  0  0 -1  0  0  0  0  0  0  0     0  0  0  0  0 -1  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  1  2  2  0  0  0 -1  0     0  0  0  0  2  3  4  3  1  0  0 -1  0     0  0  0  0  2  3  3  2  0  0  0 -1  0
     0  0  0  1  3 10 12  7  4  1  0  0  0     0  0  0  1  8 20 21 12  6  2  1  0  0     0  0  0  1  7 16 15  8  4  1  0  0  0
     0  0  0  2 14 36 43 30 13  5  2  0  0     0  0  0  5 27 63 69 43 20  8  3  0  0     0  0  0  4 21 44 43 24 11  4  1  0  0
     0  0  0  3 20 54 66 45 21  9  3  1  0     0  0  0  6 38 89100 65 29 12  4  1  0     0  0  0  5 27 57 58 36 15  6  2  0  0
     0  0  0  1 11 34 46 31 15  6  2  1  0     0  0  0  4 21 59 70 46 21  8  3  1  0     0  0  0  3 16 37 39 25 10  4  1  0  0
     0  0  0  0  4 12 16 11  6  2  1  0  0     0  0  0  1  7 19 24 16  9  3  1  0  0     0  0  0  1  5 11 13  8  4  1  0  0  0
     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  0  1  3  4  3  1  1  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0 -1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  5  4  2  1  0  0  0  0     0  0  0  1  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 12 10  5  2  1  0  0  0     0  0  0  1  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  8 15 13  8  3  1  0  0  0     0  0  0  1  1  1  1  0  0  0  0  0  0     1  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  5  9  9  5  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  2  3  3  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    657 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      181 ...     200
 ******************************************************************************

 USING  20 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ   SIGMAB   SIGMAR
   181   0  1.226  5189590    0    156      94     2  0.04071  0.09884
   182   0  1.229  5189504    0    158      91     3  0.04153  0.08124
   183   0  1.233  5188935    0    159     101     5  0.04116  0.08242
   184   0  1.231  5188783    0    152     100     5  0.04077  0.07684
   185   0  1.230  5189944    0    152      91     0  0.04081  0.11392
   186   0  1.229  5190276    0    150      88     2  0.04353  0.09930
   187   0  1.232  5188673    0    141     108     2  0.04045  0.09339
   188   0  1.228  5189489    0    153      90     4  0.03922  0.07271
   189   0  1.229  5189999    0    161      90     5  0.04329  0.09446
   190   0  1.230  5189407    0    156      97     3  0.04016  0.09116
   191   0  1.230  5188688    0    159      97     4  0.04056  0.08508
   192   0  1.234  5189224    0    157      96     1  0.04390  0.11092
   193   0  1.233  5189666    0    152      92     0  0.04108  0.10755
   194   0  1.231  5190690    0    138      85     3  0.04036  0.09675
   195   0  1.231  5190784    0    140      82     0  0.03954  0.10222
   196   0  1.230  5190908    0    128      75     1  0.03972  0.08797
   197   0  1.238  5190368    0    126      80     1  0.03943  0.08615
   198   0  1.238  5190496    0    128      82     8  0.04033  0.07339
   199   0  1.236  5191247    0    135      80     3  0.04079  0.10153
   200   0  1.241  5190579    0    141      79     2  0.04166  0.09137

   482 OUT OF    484 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.53
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.11
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.413    32.090    34.502  88.621 108.492 111.835
 REC. CELL PARAMETERS   0.041664  0.033763  0.030737  83.877  70.602  67.403
 COORDINATES OF UNIT CELL A-AXIS   -15.244   -14.235    17.788
 COORDINATES OF UNIT CELL B-AXIS    12.096    26.601    13.261
 COORDINATES OF UNIT CELL C-AXIS   -20.928    19.919   -18.858
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.036     0.015     0.005
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.045
 LAB COORDINATES OF ROTATION AXIS -0.999862  0.016385 -0.002675
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000262  0.005203  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.28   1244.11
 DETECTOR ORIGIN (PIXELS) AT                     1265.93   1237.27
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.31
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     482 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   50     31     61     71     52     25     55     70     58
 SIGMAB (degree)       0.033  0.035  0.036  0.035  0.035  0.035  0.034  0.032  0.032
 SIGMAR (degree)       0.037  0.033  0.035  0.036  0.036  0.039  0.037  0.035  0.032


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    482 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  3  3  2  1  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  1  4 11 14 10  5  2  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  1  6 18 23 17  8  3  1  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  2  2  2  1  0  0  0  0     0  0  0  0  3 12 16 12  6  2  1  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  0  2  4  5  4  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0 -1  0     0  0  0  0  2  4  4  2  1  0  0 -1  0     0  0  0  0  2  3  3  2  1  0  0  0  0
     0  0  0  1  4 12 12  8  4  1  0  0  0     0  0  0  2  9 21 21 12  5  2  0  0  0     0  0  0  2  7 15 14  7  3  1  0  0  0
     0  0  0  2 16 42 47 33 15  6  2  0  0     0  0  0  5 28 65 67 43 19  8  3  1  0     0  0  0  4 20 42 39 24 10  4  1  0  0
     0  0  0  3 21 60 75 52 25 10  3  1  0     0  0  0  6 36 90100 65 30 12  4  1  0     0  0  0  5 25 57 56 34 16  6  2  0  0
     0  0  0  2 13 41 52 37 18  7  2  1  0     0  0  0  4 23 60 69 43 20  8  3  1  0     0  0  0  3 16 38 40 22 10  4  1  0  0
     0  0  0  0  4 13 15 12  5  3  1  0  0     0  0  0  1  7 18 20 14  6  3  1  0  0     0  0  0  1  5 11 12  7  3  1  0  0  0
     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  0  1  3  3  2  1  0  0  0  0     0  0  0  0  1  1  2  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0 -1  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  5  4  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 12 10  6  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  9 17 15  8  4  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  6 11 10  5  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  3  3  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    651 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      201 ...     220
 ******************************************************************************

 USING  20 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ   SIGMAB   SIGMAR
   201   0  1.239  5190370    0    149      81     5  0.03701  0.08134
   202   0  1.240  5189271    0    161      97     6  0.04024  0.07901
   203   0  1.241  5188685    0    160     101     4  0.03872  0.08531
   204   0  1.241  5189523    0    161      96     4  0.04014  0.09771
   205   0  1.241  5189542    0    156      94     5  0.04181  0.09352
   206   0  1.238  5189201    0    152      88     6  0.04227  0.08435
   207   0  1.236  5188838    0    147      99     6  0.04233  0.08164
   208   0  1.237  5188286    0    144     109     3  0.04384  0.11979
   209   0  1.238  5188799    0    157     100     2  0.04462  0.11416
   210   0  1.238  5189484    0    152      96     3  0.04224  0.11624
   211   0  1.239  5188921    0    154      95     4  0.04573  0.11712
   212   0  1.241  5188187    0    155     102     4  0.04607  0.11007
   213   0  1.238  5188388    0    157      98     8  0.04453  0.09099
   214   0  1.235  5188523    0    146      97     0  0.04254  0.13059
   215   0  1.236  5189385    0    149      92     6  0.04221  0.10678
   216   0  1.238  5189293    0    144      97     3  0.04426  0.12685
   217   0  1.234  5189060    0    135      93     4  0.04223  0.11606
   218   0  1.241  5188530    0    152     100     2  0.04330  0.10323
   219   0  1.238  5188440    0    162     101     1  0.04320  0.11329
   220   0  1.240  5188361    0    166     104     1  0.03992  0.12891

   484 OUT OF    487 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.63
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.11
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.401    32.097    34.488  88.608 108.486 111.827
 REC. CELL PARAMETERS   0.041676  0.033753  0.030747  83.898  70.614  67.417
 COORDINATES OF UNIT CELL A-AXIS   -15.241   -14.230    17.777
 COORDINATES OF UNIT CELL B-AXIS    12.085    26.610    13.268
 COORDINATES OF UNIT CELL C-AXIS   -20.919    19.911   -18.852
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.025     0.021     0.010
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.045
 LAB COORDINATES OF ROTATION AXIS -0.999862  0.016385 -0.002675
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000275  0.004982  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.31   1244.10
 DETECTOR ORIGIN (PIXELS) AT                     1265.95   1237.56
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.19
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     484 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   50     34     65     67     59     31     54     62     51
 SIGMAB (degree)       0.034  0.035  0.036  0.035  0.035  0.035  0.033  0.032  0.032
 SIGMAR (degree)       0.040  0.033  0.034  0.036  0.037  0.042  0.040  0.037  0.032


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    484 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  3  7  9  6  3  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0 -1     0  0  0  1  3 11 14 11  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  2  8 11  9  4  2  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  2  4  4  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0 -1  0  0  0  0  0  0  0  0     0  1  0  0  0  0  0  0  0  0  0  0  0     0  0  0 -1  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0 -1 -1 -1  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  3  4  2  1  0  0 -1  0     0  0  0  0  2  4  4  2  1  0  0  0  0
     0  0  0  1  2  6  8  6  3  1  0  0  0     0  0  0  2  8 17 18 12  6  2  0  0  0     0  0  0  2  9 18 17 11  5  2  0  0  0
     0  0  0  2 11 28 36 23 11  4  1  0  0     0  0  0  5 25 56 68 41 19  8  3  1  0     0  0  0  5 24 51 56 32 14  6  2  0  0
     0  0  0  2 15 46 60 43 20  8  3  0  0     0  0  0  6 33 85100 66 30 13  4  1  0     0  0  0  6 33 74 80 49 21  9  3  1  0
     0  0  0  1 10 31 42 33 15  6  2  0  0     0  0  0  3 19 55 66 47 22  9  3  1  0     0  0  0  4 19 47 51 33 15  6  2  0  0
     0  0  0  1  4 10 15 12  6  3  1  0  0     0  0  0  1  7 16 22 16  9  3  1  0  0     0  0  0  1  6 15 16 10  5  2  0  0  0
     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  0  2  3  4  3  2  0  0  0  0     0  0  0  0  1  3  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0 -1  0 -1  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  8  7  4  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2 11 22 24 14  6  2  1  0  0     0  0  0  1  2  5  4  3  1  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  3 16 36 36 21  9  4  1  0  0     0  0  0  1  4  8  7  4  2  1  0  0  0     0  0  0  1  1  1  1  0  0  0  0  0 -1
     0  0  0  2 10 23 23 14  6  2  1  0  0     0  0  0  1  2  6  5  3  1  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  3  8  7  4  2  1  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    654 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      221 ...     240
 ******************************************************************************

 USING  20 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ   SIGMAB   SIGMAR
   221   0  1.242  5187976    0    160     113     6  0.04345  0.09929
   222   0  1.242  5187876    0    157     109     9  0.04109  0.08499
   223   0  1.248  5188268    0    154     108     6  0.04348  0.09713
   224   0  1.250  5188742    0    154     107     7  0.04167  0.09231
   225   0  1.246  5189228    0    167      93     2  0.04173  0.09086
   226   0  1.247  5188439    0    159     105     1  0.04288  0.11336
   227   0  1.248  5187658    0    163     108     1  0.03947  0.09585
   228   0  1.251  5188475    0    160     102     2  0.04177  0.09405
   229   0  1.257  5188533    0    151     103     0  0.04377  0.10984
   230   0  1.255  5188949    0    156     101     0  0.04338  0.13116
   231   0  1.256  5188741    0    155     106     0  0.04280  0.12575
   232   0  1.255  5188448    0    155     101     2  0.03955  0.08927
   233   0  1.255  5188947    0    155      99     5  0.04171  0.08543
   234   0  1.255  5190371    0    150      80     1  0.04101  0.10202
   235   0  1.263  5190616    0    148      84     0  0.04229  0.11006
   236   0  1.261  5189667    0    153      95     2  0.04034  0.10934
   237   0  1.264  5188549    0    154     105     0  0.04087  0.11253
   238   0  1.259  5189616    0    160      91     1  0.04180  0.09398
   239   0  1.260  5190092    0    158      88     3  0.04195  0.10770
   240   0  1.259  5189221    0    157      96     1  0.04476  0.11315

   508 OUT OF    510 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.54
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.11
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.411    32.077    34.499  88.624 108.494 111.851
 REC. CELL PARAMETERS   0.041672  0.033782  0.030741  83.866  70.596  67.385
 COORDINATES OF UNIT CELL A-AXIS   -15.246   -14.233    17.786
 COORDINATES OF UNIT CELL B-AXIS    12.093    26.592    13.248
 COORDINATES OF UNIT CELL C-AXIS   -20.923    19.913   -18.865
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.027     0.024     0.005
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.045
 LAB COORDINATES OF ROTATION AXIS -0.999862  0.016385 -0.002675
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000406  0.005267  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.27   1244.10
 DETECTOR ORIGIN (PIXELS) AT                     1265.74   1237.18
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.28
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     508 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   53     31     53     68     58     36     61     72     64
 SIGMAB (degree)       0.034  0.035  0.036  0.035  0.035  0.035  0.033  0.033  0.033
 SIGMAR (degree)       0.038  0.032  0.034  0.035  0.036  0.040  0.039  0.035  0.032


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    508 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  3  4  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  1  5 13 17 12  5  2  1  0 -1
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  4  5  4  2  1  0  0  0     0  0  0  1  7 22 28 21 10  4  1  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  2  3  4  3  1  0  0  0  0     0  0  0  1  5 16 19 13  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  2  6  7  4  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0 -1 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  1  2  2  1  0  0  0 -1  0     0  0  0  0  2  4  3  2  1  0 -1 -1  0     0  0  0  0  1  3  2  1  0  0  0  0  0
     0  0  0  1  5 11 13  8  3  1  0  0  0     0  0  0  2  9 19 20 12  5  2  0  0  0     0  0  0  2  8 14 12  7  3  1  0  0  0
     0  0  0  3 17 43 49 34 15  6  2  0  0     0  0  0  5 28 63 66 42 18  8  3  0  0     0  0  0  4 20 40 39 23  9  4  1  0  0
     0  0  0  4 24 68 81 56 27 11  4  1  0     0  0  0  7 38 91100 64 30 12  4  1  0     0  0  0  6 26 54 57 35 16  6  2  0  0
     0  0  0  2 17 45 56 38 17  7  3  0  0     0  0  0  4 26 61 71 46 21  9  3  0  0     0  0  0  3 17 37 41 25 11  5  1  0  0
    -1  0  0  1  5 15 17 11  6  2  0  0  0     0  0  0  1  7 19 21 15  8  3  1  0  0     0  0  0  1  5 12 11  7  3  1  0  0  0
     0  0  0  0  1  3  3  2  1  0  0 -1  0     0  0  0  0  1  3  3  2  1  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0 -1  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
    -1  0  0  1  3  4  3  2  1  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  1  0  0  0  0
     0  0  0  2  6 10  8  5  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  1  0  0  0  0
     0  0  0  2  7 14 12  8  4  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  5  9 10  6  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  2  3  3  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    675 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      241 ...     260
 ******************************************************************************

 USING  20 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ   SIGMAB   SIGMAR
   241   0  1.260  5188928    0    158      97     1  0.03979  0.10023
   242   0  1.263  5189470    0    163      94     1  0.04440  0.12031
   243   0  1.269  5189134    0    150      99     1  0.04348  0.10901
   244   0  1.269  5190420    0    152      82     1  0.04248  0.12326
   245   0  1.271  5189502    0    148      95     3  0.04139  0.10000
   246   0  1.269  5189503    0    159      93     4  0.04219  0.08592
   247   0  1.268  5189278    0    160      94     2  0.04171  0.10425
   248   0  1.268  5189595    0    144      98     6  0.04092  0.08831
   249   0  1.264  5189861    0    151      92     1  0.04067  0.11510
   250   0  1.263  5189687    0    149      91     2  0.04264  0.10202
   251   0  1.263  5189085    0    146      96     3  0.04226  0.09469
   252   0  1.262  5190119    0    146      83     4  0.03812  0.08730
   253   0  1.265  5190089    0    147      82     3  0.03996  0.08767
   254   0  1.257  5189384    0    132      92     6  0.04207  0.08611
   255   0  1.257  5189580    0    136      84     2  0.04257  0.10807
   256   0  1.255  5190114    0    135      86     2  0.04242  0.11679
   257   0  1.260  5190293    0    148      81     4  0.04484  0.11702
   258   0  1.258  5189205    0    154      93     2  0.04326  0.11388
   259   0  1.260  5188390    0    156     103     2  0.04008  0.10726
   260   0  1.258  5189002    0    157      94     0  0.04257  0.12011

   468 OUT OF    475 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.56
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.13
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.408    32.093    34.497  88.637 108.495 111.878
 REC. CELL PARAMETERS   0.041687  0.033773  0.030744  83.841  70.587  67.351
 COORDINATES OF UNIT CELL A-AXIS   -15.247   -14.240    17.776
 COORDINATES OF UNIT CELL B-AXIS    12.094    26.610    13.251
 COORDINATES OF UNIT CELL C-AXIS   -20.922    19.899   -18.877
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.007     0.029     0.009
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.045
 LAB COORDINATES OF ROTATION AXIS -0.999862  0.016385 -0.002675
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000853  0.005134  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.35   1244.06
 DETECTOR ORIGIN (PIXELS) AT                     1265.23   1237.32
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.27
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     468 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     36     61     64     43     27     58     67     56
 SIGMAB (degree)       0.033  0.035  0.036  0.035  0.034  0.035  0.034  0.032  0.032
 SIGMAR (degree)       0.037  0.031  0.033  0.035  0.037  0.039  0.037  0.036  0.033


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    468 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  3  4  3  1  0  0  0  0
     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  1  4 11 15 10  5  2  1  0  0
     0  0  0  0  1  1  2  1  1  1  0  0  0     0  0  0  0  1  3  4  3  2  1  0  0  0     0  0  0  1  5 16 22 16  7  3  1  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  3  3  2  1  0  0  0  0     0  0  0  0  3 11 16 11  5  2  1  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  2  4  6  4  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0 -1  0  0  0  0  0  0  0  0     0  0  0 -1 -1  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0 -1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  2  3  3  2  1  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0
     0  0  0  1  4 11 12  8  4  1  0  0  0     0  0  0  2  9 19 18 11  5  2  0  0  0     0  0  0  2  7 13 11  7  3  1  0  0  0
     0  0  0  2 14 40 49 33 15  6  2  0  0     0  0  0  5 26 66 70 43 19  8  3  0  0     0  0  1  4 20 43 40 23 10  4  1  0  0
     0  0  0  3 20 57 73 52 23 10  4  1 -1     0  0  0  7 37 90100 66 29 12  4  1 -1     0  0  0  6 26 56 56 35 15  6  2  0 -1
     0  0  0  2 13 39 50 34 16  7  2  0  0     0  0  0  4 23 59 68 43 20  8  3  0  0     0  0  0  3 16 36 38 22 10  4  1  0  0
     0  0  0  1  4 11 14 11  5  2  0  0  0     0  0  0  1  7 16 19 14  6  3  1  0  0     0  0  0  1  4 10 10  7  3  1  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  1  3  3  2  1  0  0  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  3  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  6 10  8  5  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 14 12  7  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  8  8  5  2  1  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0 -1  0  0  0  1  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    633 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      261 ...     280
 ******************************************************************************

 USING  20 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ   SIGMAB   SIGMAR
   261   0  1.260  5189227    0    154      91     1  0.04638  0.11436
   262   0  1.260  5189897    0    152      82     3  0.04330  0.08982
   263   0  1.260  5189719    0    156      88     2  0.04057  0.11110
   264   0  1.258  5189436    0    158      90     5  0.04052  0.09664
   265   0  1.262  5190045    0    156      85     4  0.03959  0.08908
   266   0  1.261  5189984    0    150      86     4  0.04102  0.08872
   267   0  1.263  5189950    0    156      87     4  0.04004  0.08706
   268   0  1.261  5189460    0    157      89     0  0.04137  0.08525
   269   0  1.259  5189804    0    164      85     2  0.04015  0.07556
   270   0  1.259  5189378    0    167      93     1  0.03883  0.07594
   271   0  1.264  5189536    0    168      87     3  0.03935  0.07047
   272   0  1.268  5189860    0    159      94     2  0.04428  0.09575
   273   0  1.265  5189526    0    162      94     3  0.04034  0.08518
   274   0  1.266  5189271    0    164      93     4  0.04036  0.07481
   275   0  1.268  5189429    0    163      89     1  0.04131  0.09853
   276   0  1.267  5189197    0    180      93     1  0.04102  0.08651
   277   0  1.272  5188607    0    171     101     0  0.04125  0.11141
   278   0  1.274  5187604    0    177     112     2  0.03968  0.09243
   279   0  1.271  5188592    0    168     106     0  0.04121  0.10449
   280   0  1.266  5190321    0    153      78     0  0.03849  0.09086

   502 OUT OF    506 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.50
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.411    32.087    34.497  88.631 108.496 111.858
 REC. CELL PARAMETERS   0.041675  0.033773  0.030744  83.855  70.590  67.375
 COORDINATES OF UNIT CELL A-AXIS   -15.247   -14.229    17.790
 COORDINATES OF UNIT CELL B-AXIS    12.096    26.608    13.241
 COORDINATES OF UNIT CELL C-AXIS   -20.921    19.899   -18.878
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.007     0.030    -0.011
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.045
 LAB COORDINATES OF ROTATION AXIS -0.999862  0.016385 -0.002675
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000644  0.005269  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.29   1244.11
 DETECTOR ORIGIN (PIXELS) AT                     1265.44   1237.18
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.27
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     502 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   50     33     62     74     65     33     56     62     57
 SIGMAB (degree)       0.034  0.035  0.035  0.035  0.035  0.035  0.034  0.033  0.033
 SIGMAR (degree)       0.035  0.030  0.034  0.035  0.036  0.038  0.038  0.036  0.031


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    502 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  4 11 13 10  5  2  1  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  3  3  2  1  0  0  0  0     0  0  0  1  5 16 22 14  7  3  1  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  0  4 11 15 10  5  2  1  0 -1
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  2  5  5  3  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0
     0  0 -1  0  0  0  1  0  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  2  1  0  0  0  0  0     0  0  0  0  2  3  3  2  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0
     0  0  0  1  3 10 12  7  3  1  0  0  0     0  0  0  2  7 19 19 11  5  2  0  0  0     0  0  0  2  6 12 12  6  3  1  0  0  0
     0  0  0  2 14 37 44 31 15  6  2  0  0     0  0  0  6 27 62 66 40 18  8  3  0  0     0  0  0  5 21 43 41 22 10  4  1  0  0
     0  0  0  3 20 56 71 47 23 10  4  1  0     0  0  0  7 39 90100 63 30 12  4  1  0     0  0  0  6 28 58 59 35 16  7  2  0  0
     0  0  0  2 12 37 48 32 16  6  2  0 -1     0  0  0  4 21 57 69 45 21  8  3  1  0     0  0  0  3 16 37 41 25 12  5  1  0  0
     0  0  0  1  5 12 14 10  5  2  0  0  0     0  0  0  1  8 18 20 15  7  3  1  0  0     0  0  0  1  6 12 13  9  4  1  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  1  3  3  2  1  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  3  3  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 10  9  5  2  1  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 14 12  7  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  9  8  5  2  1  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0
     0  0  0  0  1  3  3  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    669 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      281 ...     300
 ******************************************************************************

 USING  20 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ   SIGMAB   SIGMAR
   281   0  1.270  5190961    0    130      74     2  0.03968  0.08503
   282   0  1.274  5191885    0    123      64     2  0.04412  0.10524
   283   0  1.282  5191384    0    117      76     0  0.04660  0.12539
   284   0  1.283  5191062    0    124      78     0  0.04638  0.13595
   285   0  1.285  5190089    0    138      87     1  0.04355  0.12322
   286   0  1.281  5189469    0    145      88     9  0.03905  0.07264
   287   0  1.281  5190317    0    156      85     1  0.03999  0.10631
   288   0  1.277  5190070    0    153      83     1  0.04582  0.10516
   289   0  1.273  5189063    0    146      96     1  0.03987  0.11337
   290   0  1.273  5189457    0    153      91     7  0.03790  0.07790
   291   0  1.273  5189757    0    150      88     2  0.04187  0.10917
   292   0  1.269  5189791    0    158      90     2  0.04377  0.11649
   293   0  1.268  5189755    0    155      96     1  0.04521  0.12862
   294   0  1.271  5190314    0    143      85     6  0.04001  0.09407
   295   0  1.266  5190650    0    150      79     4  0.03702  0.07723
   296   0  1.269  5189866    0    143      88     4  0.03940  0.07402
   297   0  1.275  5189846    0    143      95     5  0.04040  0.09121
   298   0  1.277  5190788    0    143      85     1  0.04351  0.11998
   299   0  1.273  5190062    0    136      92     3  0.04032  0.10673
   300   0  1.278  5190623    0    138      82     0  0.03835  0.10946

   453 OUT OF    457 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.54
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.14
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.413    32.088    34.502  88.626 108.496 111.890
 REC. CELL PARAMETERS   0.041683  0.033780  0.030739  83.848  70.588  67.342
 COORDINATES OF UNIT CELL A-AXIS   -15.247   -14.232    17.789
 COORDINATES OF UNIT CELL B-AXIS    12.100    26.613    13.228
 COORDINATES OF UNIT CELL C-AXIS   -20.926    19.897   -18.884
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.002     0.030    -0.019
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.045
 LAB COORDINATES OF ROTATION AXIS -0.999862  0.016385 -0.002675
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000755  0.005175  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.26   1244.09
 DETECTOR ORIGIN (PIXELS) AT                     1265.27   1237.29
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.31
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     453 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     28     52     61     44     28     57     71     55
 SIGMAB (degree)       0.033  0.035  0.036  0.035  0.034  0.035  0.033  0.032  0.032
 SIGMAR (degree)       0.036  0.031  0.034  0.035  0.036  0.038  0.037  0.036  0.033


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    453 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  2  5  6  3  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  4  3  1  0  0  0  0     0  0  0  0  5 15 21 15  7  3  1  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  2  4  6  4  2  1  0  0  0     0  0  0  1  8 24 31 23 10  4  1  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  2  4  5  3  1  0  0  0  0     0  0  0  1  6 17 21 14  7  3  1  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  2  6  7  4  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0


     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  3  2  1  0  0  0  0  0     0  0  0  0  2  3  2  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  1  6 14 16  8  4  1  0  0  0     0  0  0  2  9 18 17  9  4  1  0  0  0     0  0  0  1  5 10  8  4  2  0  0  0  0
     0  0  0  3 17 43 54 38 17  7  2  0  0     0  0  0  5 26 59 66 43 18  7  2  0  0     0  0  0  4 16 31 30 18  8  3  1  0  0
     0  0  0  3 24 66 83 57 26 11  4  1  0     0  0  0  6 39 87100 65 30 12  4  1  0     0  0  0  4 23 44 46 27 13  5  2  0  0
     0  0  0  2 17 45 54 35 17  7  2  0  0     0  0  0  5 26 61 67 41 19  8  3  0  0     0  0  0  3 14 30 31 18  8  3  1  0  0
     0  0  0  1  5 14 15 10  5  2  0  0  0     0  0  0  1  8 17 19 13  6  2  1  0  0     0  0  0  1  5  9  8  5  2  1  0  0  0
     0  0  0  0  2  2  2  1  1  0  0  0  0     0  0  0  0  1  3  3  2  1  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  1  4  6  5  3  1  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  1  0  1  0  0  0  0  0  0
     0  0  0  1  5  8  7  4  2  1  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  1  1  0  1  0  0  0  0  0
     0  0  0  0  3  5  5  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    620 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      301 ...     320
 ******************************************************************************

 USING  20 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ   SIGMAB   SIGMAR
   301   0  1.273  5190650    0    142      81     2  0.04215  0.09728
   302   0  1.280  5189907    0    152      88     1  0.04285  0.09572
   303   0  1.273  5188586    0    167     101     4  0.04296  0.10070
   304   0  1.274  5187725    0    174     112     2  0.04048  0.11234
   305   0  1.278  5188283    0    171     102     2  0.04300  0.10991
   306   0  1.279  5187809    0    180     115     2  0.04083  0.09761
   307   0  1.281  5188335    0    171     107     5  0.04003  0.10818
   308   0  1.281  5188494    0    168     104     7  0.03895  0.08012
   309   0  1.278  5188453    0    162     108     4  0.04128  0.11040
   310   0  1.284  5188062    0    153     107     4  0.04336  0.10744
   311   0  1.279  5188256    0    155     111     7  0.04443  0.10075
   312   0  1.282  5187910    0    175     110     8  0.04139  0.08490
   313   0  1.284  5188522    0    159     101     7  0.04079  0.08179
   314   0  1.275  5189766    0    170      86     3  0.04035  0.09696
   315   0  1.281  5189045    0    161      92     3  0.04035  0.09461
   316   0  1.278  5188155    0    148     101     2  0.03959  0.11235
   317   0  1.279  5188902    0    155      93     4  0.04412  0.09944
   318   0  1.280  5189658    0    162      91     3  0.04268  0.11319
   319   0  1.272  5188854    0    166      96     3  0.04144  0.11700
   320   0  1.275  5188295    0    163     101     2  0.04174  0.12839

   514 OUT OF    518 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.54
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.411    32.086    34.497  88.619 108.495 111.880
 REC. CELL PARAMETERS   0.041681  0.033779  0.030744  83.860  70.593  67.355
 COORDINATES OF UNIT CELL A-AXIS   -15.245   -14.230    17.790
 COORDINATES OF UNIT CELL B-AXIS    12.094    26.614    13.228
 COORDINATES OF UNIT CELL C-AXIS   -20.924    19.894   -18.880
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.007     0.029    -0.017
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.045
 LAB COORDINATES OF ROTATION AXIS -0.999862  0.016385 -0.002675
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000576  0.005176  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.34   1244.09
 DETECTOR ORIGIN (PIXELS) AT                     1265.58   1237.29
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.28
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     514 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   54     37     61     71     64     35     57     65     58
 SIGMAB (degree)       0.033  0.035  0.035  0.035  0.035  0.035  0.034  0.033  0.032
 SIGMAR (degree)       0.036  0.031  0.032  0.034  0.037  0.039  0.038  0.035  0.031


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    514 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  3  4  2  1  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  0  4 12 15 12  5  2  1  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  3  4  3  2  1  1  0  0     0  0  0  1  6 18 24 18  9  4  1  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  3  4  2  1  0  0  0  0     0  0  0  1  5 13 17 11  5  2  1  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  1  2  5  6  4  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  2  3  3  1  1  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0
     0  0  0  1  4 12 13  7  4  1  0  0  0     0  0  0  2  8 20 19 10  5  2  0  0  0     0  0  0  2  6 12 11  5  2  1  0  0  0
     0  0  0  2 14 40 49 34 16  6  2  0  0     0  0  0  5 27 63 67 43 20  8  3  1  0     0  0  1  4 18 39 37 22 10  4  1  0  0
     0  0  0  3 23 62 76 53 25 11  4  1  0     0  0  0  7 38 89100 66 30 13  4  1  0     0  0  1  5 25 53 55 33 15  6  2  0  0
     0  0  0  2 16 41 52 34 16  6  2  0  0     0  0  0  4 24 58 70 45 20  8  3  1  0     0  0  0  3 14 33 37 23 10  4  1  0  0
     0  0  0  1  6 14 16 11  5  2  1  0  0     0  0  0  2  8 19 20 15  7  3  1  0  0     0  0  0  1  4 11 11  8  4  1  0  0  0
     0  0  0  0  2  3  3  2  1  0  0  0  0     0  0  0  0  2  3  3  2  1  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  8  7  4  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  6 11 10  6  3  1  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  6  7  4  1  1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  0  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0 -1 -1 -1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    680 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      321 ...     340
 ******************************************************************************

 USING  20 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ   SIGMAB   SIGMAR
   321   0  1.282  5188006    0    166     111     3  0.04084  0.12231
   322   0  1.282  5188430    0    160     102     6  0.04271  0.09739
   323   0  1.284  5189320    0    154      94     0  0.04401  0.13853
   324   0  1.278  5189471    0    149      96     1  0.03901  0.11450
   325   0  1.277  5189307    0    140      98     3  0.04421  0.12564
   326   0  1.273  5189983    0    138      86     7  0.04178  0.09092
   327   0  1.271  5190068    0    145      88     7  0.04389  0.09883
   328   0  1.273  5189549    0    140      96     5  0.04237  0.11120
   329   0  1.279  5189198    0    144     101     5  0.04477  0.11079
   330   0  1.283  5189545    0    148      94     4  0.04111  0.11947
   331   0  1.277  5190132    0    143      83     4  0.04293  0.10685
   332   0  1.278  5190200    0    133      85     2  0.04280  0.12449
   333   0  1.276  5190494    0    122      77     0  0.04208  0.11696
   334   0  1.277  5190802    0    130      76     0  0.04092  0.12277
   335   0  1.278  5190598    0    139      83     0  0.04115  0.12639
   336   0  1.283  5190051    0    144      86     2  0.04178  0.11777
   337   0  1.280  5189242    0    154      96     6  0.03879  0.09599
   338   0  1.276  5188410    0    159     104     2  0.03976  0.09964
   339   0  1.277  5188177    0    164     111     5  0.03961  0.08733
   340   0  1.279  5188549    0    160     107     3  0.03969  0.10062

   474 OUT OF    480 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.68
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.16
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.401    32.102    34.472  88.619 108.478 111.778
 REC. CELL PARAMETERS   0.041660  0.033736  0.030759  83.907  70.624  67.466
 COORDINATES OF UNIT CELL A-AXIS   -15.238   -14.205    17.799
 COORDINATES OF UNIT CELL B-AXIS    12.080    26.625    13.256
 COORDINATES OF UNIT CELL C-AXIS   -20.915    19.881   -18.858
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.008     0.023    -0.027
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.045
 LAB COORDINATES OF ROTATION AXIS -0.999862  0.016385 -0.002675
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000112  0.005038  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.28   1244.07
 DETECTOR ORIGIN (PIXELS) AT                     1266.14   1237.46
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.15
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     474 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     30     59     65     47     32     56     67     57
 SIGMAB (degree)       0.034  0.035  0.035  0.035  0.033  0.033  0.034  0.033  0.033
 SIGMAR (degree)       0.041  0.031  0.034  0.036  0.040  0.046  0.040  0.037  0.033


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    474 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  0  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  4  5  4  2  1  0  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  1  0  1  1  0  0  0  0     0  0  0  0  2  5  7  6  3  1  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  1  4  5  5  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  1  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0 -1 -1 -1 -1  0  0  0  0     0  0  0  0 -1 -1 -1 -1 -1 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0
     0  0  0  0  0  0  1  1  0  0  0 -1  0     0  0  0  0  1  2  2  1  1  0  0 -1  0     0  0  0  0  2  3  4  2  1  0 -1  0  0
     0  0  0  0  2  6  8  5  2  1  0  0  0     0  0  0  1  7 15 18 10  5  2  0  0  0     0  0  0  2 10 20 22 13  5  2  0  0  0
     0  0  0  1  9 24 27 18  8  3  1  0  0     0  0  0  4 25 61 66 42 18  7  3  0  0     0  0  0  6 29 66 68 43 18  7  2  0  0
     0  0  0  2 12 34 43 30 14  5  2  0 -1     0  0  0  6 33 84100 68 31 13  5  1  0     0  0  1  8 39 87 96 64 28 12  4  1  0
     0  0  0  1  7 22 27 19  9  4  1  0  0     0  0  0  4 23 55 65 45 21  9  3  0  0     0  0  0  5 27 59 64 41 19  8  3  1  0
     0  0  0  1  3  7  9  6  3  1  0  0  0     0  0  0  2  8 18 21 14  7  3  1  0  0     0  0  0  2 10 19 22 14  6  2  1  0  0
     0  0  0  0  1  2  1  1  0  0  0 -1  0     0  0  0  1  2  4  3  2  1  0  0 -1  0     0  0  0  1  2  4  3  2  1  0  0 -1  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0
     0  0  0 -1  0  0  0  0  0  0  0  0  0     0  0  0 -1  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0 -1 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -2  0
     0  0  0  1  2  3  3  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  2  7 14 14  8  3  1  0  0  0     0  0  0  0  2  4  5  2  1  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  3 18 38 40 25 10  4  1  0  0     0  0  0  1  6 13 14 10  4  1  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0
     0  0  0  5 24 51 58 36 16  6  2  0  0     0  0  0  2  9 20 23 13  6  2  1  0  0     0  0  0  1  2  4  4  2  1  0  0  0  0
     0  0  0  3 16 36 39 23 10  4  1  0  0     0  0  0  1  6 14 15  9  4  1  1  0  0     0  0  0  0  2  3  3  1  1  0  0  0  0
     0  0  0  1  5 11 12  7  3  1  0  0  0     0  0  0  0  2  4  5  2  1  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  1  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    637 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      341 ...     360
 ******************************************************************************

 USING  20 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ   SIGMAB   SIGMAR
   341   0  1.276  5189137    0    161      97     1  0.04115  0.10253
   342   0  1.278  5189997    0    159      88     0  0.04011  0.10766
   343   0  1.280  5189235    0    152      94     2  0.03904  0.08804
   344   0  1.274  5188521    0    153     102     3  0.03983  0.08619
   345   0  1.275  5188922    0    155      94     1  0.04038  0.08346
   346   0  1.273  5188787    0    152      96     2  0.04043  0.09368
   347   0  1.274  5188618    0    155     103     2  0.04104  0.09920
   348   0  1.270  5188872    0    161      99     3  0.04026  0.09434
   349   0  1.273  5189410    0    158      95     1  0.04252  0.11413
   350   0  1.273  5190199    0    155      85     0  0.04354  0.12284
   351   0  1.273  5189869    0    149      91     1  0.04284  0.12051
   352   0  1.270  5189581    0    147      93     2  0.04035  0.10617
   353   0  1.273  5189649    0    152      93     0  0.04145  0.12531
   354   0  1.276  5189474    0    152      96     2  0.04064  0.10759
   355   0  1.272  5189469    0    155      93     3  0.03977  0.09962
   356   0  1.268  5189282    0    164      95     2  0.03823  0.10330
   357   0  1.272  5188739    0    164     102     2  0.04257  0.11064
   358   0  1.273  5188825    0    153     100     1  0.03981  0.12882
   359   0  1.270  5189937    0    150      88     1  0.04103  0.13869
   360   0  1.272  5190296    0    136      89     4  0.04234  0.11570

   498 OUT OF    502 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.59
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.406    32.101    34.492  88.647 108.498 111.854
 REC. CELL PARAMETERS   0.041684  0.033759  0.030750  83.837  70.581  67.373
 COORDINATES OF UNIT CELL A-AXIS   -15.239   -14.219    17.796
 COORDINATES OF UNIT CELL B-AXIS    12.089    26.633    13.229
 COORDINATES OF UNIT CELL C-AXIS   -20.925    19.867   -18.898
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.041     0.031    -0.033
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.045
 LAB COORDINATES OF ROTATION AXIS -0.999862  0.016385 -0.002675
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000906  0.005209  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.35   1244.11
 DETECTOR ORIGIN (PIXELS) AT                     1265.16   1237.27
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.25
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     498 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   48     30     58     71     62     32     59     66     56
 SIGMAB (degree)       0.034  0.034  0.035  0.035  0.035  0.036  0.035  0.033  0.032
 SIGMAR (degree)       0.037  0.031  0.034  0.035  0.038  0.042  0.039  0.035  0.032


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    498 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  0  3  9 11  8  4  1  0  0  0
     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  1  4 13 17 12  6  2  1  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  0  3  9 12  8  4  2  0  0 -1
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  1  2  1  0  0  0  0  0     0  0  0  0  2  4  5  3  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0 -1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  2  3  3  2  1  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0
     0  0  0  0  3  9 12  7  3  1  0  0  0     0  0  0  2  8 17 21 11  5  2  0  0  0     0  0  0  1  7 14 15  8  4  1  0  0  0
     0  0  0  2 11 32 41 29 13  5  2  0 -1     0  0  0  4 24 61 69 45 19  8  3  0 -1     0  0  1  5 21 48 51 31 13  6  2  0 -1
     0  0  0  2 14 47 59 42 20  8  3  1  0     0  0  1  6 33 87100 66 30 13  4  1  0     0  0  1  6 29 66 71 45 21  8  3  1  0
     0  0  0  1 12 31 38 27 13  5  2  0  0     0  0  0  3 24 56 66 43 21  8  3  1  0     0  0  0  3 19 42 47 29 13  5  2  0  0
     0  0  0  1  5 11 13  8  4  2  1  0  0     0  0  0  2  8 18 22 14  7  3  1  0  0     0  0  0  1  6 13 15 10  5  2  1  0  0
     0  0  0  0  2  3  2  2  0  0  0  0  0     0  0  0  0  2  4  4  2  1  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0 -1 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  4  4  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  8 16 14  9  4  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3 10 22 22 12  6  2  1  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 13 13  8  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  4  4  3  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    665 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      361 ...     380
 ******************************************************************************

 USING  20 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ   SIGMAB   SIGMAR
   361   0  1.272  5190075    0    125      90     1  0.04023  0.11312
   362   0  1.273  5190885    0    128      77     4  0.03995  0.09824
   363   0  1.269  5190231    0    135      85     4  0.04325  0.11672
   364   0  1.272  5189864    0    128      89     3  0.04408  0.10994
   365   0  1.272  5190045    0    137      84     6  0.04032  0.09293
   366   0  1.277  5188916    0    145      93     2  0.04049  0.11881
   367   0  1.274  5188409    0    151     107     1  0.04594  0.11678
   368   0  1.277  5188748    0    153     104     2  0.04474  0.11745
   369   0  1.270  5189475    0    156      95     0  0.04409  0.14143
   370   0  1.271  5189491    0    154      90     0  0.04360  0.13450
   371   0  1.275  5189322    0    153      94     4  0.04366  0.11753
   372   0  1.271  5189319    0    156      92     1  0.04286  0.11024
   373   0  1.267  5189372    0    156      96     3  0.04276  0.11300
   374   0  1.268  5189101    0    158      97     1  0.04291  0.11052
   375   0  1.266  5188633    0    167     106     0  0.04793  0.12919
   376   0  1.268  5188246    0    170     109     0  0.04133  0.12361
   377   0  1.268  5188412    0    159     105     1  0.04286  0.11502
   378   0  1.269  5189794    0    159      92     1  0.04217  0.11424
   379   0  1.268  5190296    0    153      85     2  0.04172  0.11565
   380   0  1.271  5189621    0    148      92     6  0.04283  0.09163

   480 OUT OF    484 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.60
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.14
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.403    32.096    34.499  88.641 108.499 111.916
 REC. CELL PARAMETERS   0.041708  0.033780  0.030745  83.820  70.576  67.309
 COORDINATES OF UNIT CELL A-AXIS   -15.237   -14.234    17.781
 COORDINATES OF UNIT CELL B-AXIS    12.094    26.630    13.219
 COORDINATES OF UNIT CELL C-AXIS   -20.931    19.873   -18.897
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.038     0.026    -0.024
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.045
 LAB COORDINATES OF ROTATION AXIS -0.999862  0.016385 -0.002675
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.001206  0.005135  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.40   1244.06
 DETECTOR ORIGIN (PIXELS) AT                     1264.82   1237.31
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.27
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     480 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   51     33     57     65     54     29     56     66     57
 SIGMAB (degree)       0.033  0.034  0.035  0.035  0.034  0.035  0.034  0.033  0.032
 SIGMAR (degree)       0.038  0.033  0.034  0.035  0.034  0.039  0.038  0.037  0.034


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    480 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  0  2  4  5  4  2  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  4  4  3  1  1  0  0  0     0  0  0  1  5 15 19 14  6  3  1  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  2  5  6  4  2  1  0  0  0     0  0  0  1  7 22 29 20  9  4  1  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  2  4  5  3  1  1  0  0  0     0  0  0  1  6 16 19 13  6  2  1  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  2  3  2  1  0  0  0  0  0     0  0  0  1  3  7  7  4  2  0  0  0  0
     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  1  1  1  1  0  0  0  0  0  0     0  0  0  1  2  2  2  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0


     0  0  0  0  0  0 -1  0  0  0  0  0  0     0  0  0 -1  0 -1 -1 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0
     0  0  0  0  1  2  2  1  0  0 -1 -1  0     0  0  0  0  1  3  2  1  0  0  0 -1  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  1  5 12 14  9  4  2  0  0  0     0  0  0  1  8 17 17 10  4  1  0  0  0     0  0  0  1  6 10  8  5  2  0  0  0  0
     0  0  0  3 14 42 53 36 15  7  2  0  0     0  0  0  5 25 62 68 40 17  7  2  0  0     0  0  0  4 18 37 34 19  7  3  1  0  0
     0  0  0  3 21 62 78 54 25 11  4  1  0     0  0  0  6 37 89100 63 29 12  4  1  0     0  0  0  5 24 50 51 30 13  5  2  0  0
     0  0  0  3 16 42 50 34 16  6  2  0  0     0  0  0  4 24 57 64 42 20  7  2  0  0     0  0  0  3 15 29 32 20  9  4  1  0  0
     0  0  0  2  6 13 14 10  4  2  0  0  0     0  0  0  1  7 16 19 13  6  2  0  0  0     0  0  0  1  4  9 10  7  3  1  0  0  0
     0  0  0  1  2  2  2  1  0  0  0 -1  0     0  0  0  0  1  2  3  2  1  0  0 -1  0     0  0  0  0  1  2  2  1  1  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
    -1  0  0  0  2  2  2  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  1  0  1  0  0  0  0  0
     0  0  0  1  5  8  7  3  1  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  6 10  9  5  3  1  0  0 -1     0  0  0  0  1  1  1  1  0  0  0  0 -1     0  0  0  0  1  0  0  0  0  0  0  0  0
     0  0  0  1  3  6  7  4  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  1  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    639 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      381 ...     400
 ******************************************************************************

 USING  20 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ   SIGMAB   SIGMAR
   381   0  1.266  5189403    0    152      93     2  0.04305  0.12029
   382   0  1.274  5189517    0    144      94     3  0.04038  0.11394
   383   0  1.274  5189889    0    148      86     1  0.04136  0.13137
   384   0  1.271  5189693    0    150      94     4  0.04359  0.12464
   385   0  1.271  5189822    0    147      91     0  0.04404  0.14667
   386   0  1.270  5188877    0    159      97     1  0.04465  0.13843
   387   0  1.266  5189027    0    144      98     1  0.04432  0.12575
   388   0  1.266  5189743    0    143      90     4  0.04177  0.10462
   389   0  1.267  5189895    0    141      89     0  0.04207  0.12770
   390   0  1.272  5189852    0    147      95     3  0.03741  0.11123
   391   0  1.267  5189673    0    153      90     4  0.04126  0.09616
   392   0  1.263  5189472    0    156      89     1  0.03786  0.10932
   393   0  1.265  5188677    0    153      98     1  0.04271  0.11060
   394   0  1.268  5189197    0    159      96     3  0.04006  0.11108
   395   0  1.263  5189305    0    145      98     1  0.04037  0.11909
   396   0  1.259  5189451    0    132      96     9  0.03748  0.07074
   397   0  1.259  5190695    0    136      82     1  0.03980  0.11142
   398   0  1.259  5190796    0    137      79     1  0.03987  0.09948
   399   0  1.267  5189997    0    126      86     5  0.04146  0.09461
   400   0  1.267  5189935    0    144      93     0  0.04367  0.14915

   450 OUT OF    456 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.54
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.15
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.413    32.098    34.494  88.620 108.495 111.858
 REC. CELL PARAMETERS   0.041671  0.033761  0.030746  83.867  70.595  67.378
 COORDINATES OF UNIT CELL A-AXIS   -15.235   -14.219    17.810
 COORDINATES OF UNIT CELL B-AXIS    12.093    26.636    13.215
 COORDINATES OF UNIT CELL C-AXIS   -20.935    19.874   -18.883
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.037     0.012    -0.052
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.045
 LAB COORDINATES OF ROTATION AXIS -0.999862  0.016385 -0.002675
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000578  0.005122  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.30   1244.16
 DETECTOR ORIGIN (PIXELS) AT                     1265.54   1237.43
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.29
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     450 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   48     29     55     60     43     29     55     68     55
 SIGMAB (degree)       0.034  0.035  0.035  0.035  0.035  0.036  0.035  0.033  0.033
 SIGMAR (degree)       0.037  0.032  0.034  0.034  0.038  0.042  0.039  0.036  0.033


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    450 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  3  4  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  0  4 10 14 11  5  2  1  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  3  3  1  0  0  0  0     0  0  0  0  4 15 21 16  8  3  1  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  0  4 10 14 10  5  2  1  0  0
     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  1  2  5  5  3  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  1  1  2  2  1  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  1  1  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0


     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0 -1  0     0  0  0  0  2  3  3  2  1  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0
     0  0  0  1  4 11 12  8  4  2  0  0  0     0  0  0  2  8 19 20 12  6  2  0  0  0     0  0  0  2  7 13 14  9  4  1  0  0  0
     0  0  0  2 14 37 48 34 15  6  2  0  0     0  0  0  5 27 61 70 46 21  8  3  0  0     0  0  0  5 21 43 46 29 12  5  1  0  0
     0  0  0  3 19 53 69 50 24 10  3  1  0     0  0  0  6 35 84100 67 32 13  5  1  0     0  0  1  5 26 59 63 41 19  7  3  0  0
     0  0  0  2 13 35 43 29 15  6  2  0  0     0  0  0  4 24 56 64 41 20  8  3  0  0     0  0  0  3 17 36 40 25 11  4  2  0  0
     0  0  0  2  5 12 13  8  3  1  0  0  0     0  0  0  2  8 18 19 11  5  2  0  0  0     0  0  0  1  5 12 13  8  3  1  0  0  0
     0  0  0  1  2  3  2  1  0  0  0 -1  0     0  0  0  1  2  3  3  2  1  0  0 -1  0     0  0  0  0  1  2  2  1  1  0  0  0  0
     0  0  0  0  1  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  4  3  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 13 12  8  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  9 18 17 10  4  2  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  6 11 12  6  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  4  4  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    633 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      401 ...     420
 ******************************************************************************

 USING  20 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ   SIGMAB   SIGMAR
   401   0  1.265  5188843    0    157     102     1  0.04489  0.12551
   402   0  1.263  5187961    0    161     107     1  0.04031  0.10129
   403   0  1.259  5188093    0    174     103     0  0.04257  0.10741
   404   0  1.262  5187857    0    169     109     0  0.04128  0.11287
   405   0  1.259  5188189    0    173     100     0  0.04016  0.10710
   406   0  1.264  5187525    0    166     111     0  0.04186  0.10120
   407   0  1.266  5188081    0    160     108     2  0.03884  0.09809
   408   0  1.258  5188773    0    151     100     1  0.04099  0.11359
   409   0  1.257  5189884    0    150      84     1  0.04537  0.11511
   410   0  1.258  5189259    0    143      90     1  0.04179  0.10189
   411   0  1.260  5188930    0    141      95     5  0.04004  0.09228
   412   0  1.255  5189153    0    148      91     4  0.03961  0.07909
   413   0  1.255  5189218    0    147      92     4  0.04010  0.09225
   414   0  1.257  5189209    0    157      91     2  0.04200  0.09450
   415   0  1.260  5188853    0    158      99     0  0.04364  0.11716
   416   0  1.256  5188146    0    154     104     3  0.03854  0.08590
   417   0  1.252  5188853    0    146     101     1  0.03919  0.09131
   418   0  1.248  5190703    0    145      81     1  0.04189  0.09011
   419   0  1.250  5189883    0    138      95     2  0.04244  0.10878
   420   0  1.247  5190384    0    141      81     4  0.04003  0.08080

   500 OUT OF    504 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.57
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.408    32.099    34.494  88.622 108.492 111.886
 REC. CELL PARAMETERS   0.041687  0.033767  0.030745  83.856  70.595  67.349
 COORDINATES OF UNIT CELL A-AXIS   -15.232   -14.223    17.802
 COORDINATES OF UNIT CELL B-AXIS    12.093    26.639    13.208
 COORDINATES OF UNIT CELL C-AXIS   -20.938    19.868   -18.886
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.046     0.010    -0.051
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.045
 LAB COORDINATES OF ROTATION AXIS -0.999862  0.016385 -0.002675
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000817  0.005182  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.30   1244.03
 DETECTOR ORIGIN (PIXELS) AT                     1265.23   1237.22
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.28
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     500 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   52     33     61     73     63     31     56     65     52
 SIGMAB (degree)       0.034  0.035  0.035  0.035  0.034  0.035  0.035  0.033  0.032
 SIGMAR (degree)       0.037  0.032  0.034  0.035  0.037  0.040  0.038  0.037  0.035


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    500 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  4  4  3  1  0  0  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  0  4 12 17 12  6  2  1  0  0
     0  0  0  0  0  0  1  1  1  0  0  0  0     0  0  0  0  1  3  4  3  2  1  0  0  0     0  0  0  1  6 19 25 18  9  4  1  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  3  4  2  1  0  0  0  0     0  0  0  1  4 14 17 12  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  1  3  5  5  4  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  1  2  2  2  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0 -1  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0 -1 -1  0     0  0  0  0  2  3  3  1  0  0 -1 -1  0     0  0  0  0  1  2  2  1  0  0  0  0  0
     0  0  0  1  5 13 14  8  4  1  0  0  0     0  0  0  2  8 21 19 11  5  2  0  0  0     0  0  0  2  6 13 11  6  2  1  0  0  0
     0  0  0  2 14 40 52 33 16  7  2  0  0     0  0  1  5 26 61 70 41 20  8  3  0  0     0  0  0  4 19 37 39 23 10  4  1  0  0
     0  0  0  3 22 59 73 51 25 10  4  1  0     0  0  1  6 38 88100 65 31 13  4  1  0     0  0  1  5 26 54 57 35 16  7  2  0  0
     0  0  0  2 14 41 51 33 15  6  2  0  0     0  0  0  4 24 64 72 44 19  8  3  1  0     0  0  0  3 16 40 43 25 11  5  2  0  0
     0  0  0  2  6 13 14  9  4  2  0  0  0     0  0  0  2  8 18 20 13  6  3  1  0  0     0  0  0  1  5 12 13  8  4  2  0  0  0
     0  0  0  1  2  3  3  2  0  0  0 -1  0     0  0  0  1  2  3  3  2  1  0  0 -1  0     0  0  0  0  1  2  2  1  1  0  0  0  0
     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  3  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  9  8  5  2  1  0  0  0     0  0  0  0  0  1  0  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  2  7 13 12  8  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  4  9 10  6  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  3  3  3  1  1  0  0  0     0  0  0  0  1  0  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    662 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      421 ...     440
 ******************************************************************************

 USING  20 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ   SIGMAB   SIGMAR
   421   0  1.251  5190650    0    141      76     2  0.04151  0.09364
   422   0  1.245  5190724    0    132      74     0  0.04121  0.10918
   423   0  1.247  5190718    0    141      76     0  0.04261  0.10584
   424   0  1.244  5189636    0    133      91     0  0.04110  0.10112
   425   0  1.245  5189792    0    134      88     1  0.04178  0.09792
   426   0  1.242  5190725    0    134      78     1  0.04228  0.09892
   427   0  1.238  5189875    0    130      86     0  0.04004  0.11408
   428   0  1.238  5190384    0    123      81     0  0.04163  0.12003
   429   0  1.237  5191074    0    139      74     1  0.04271  0.10930
   430   0  1.237  5190070    0    148      87     1  0.04056  0.11106
   431   0  1.239  5189273    0    153      94     3  0.04158  0.09133
   432   0  1.237  5188724    0    165     101     2  0.04022  0.09318
   433   0  1.235  5189527    0    164      89     3  0.03991  0.08926
   434   0  1.233  5189526    0    166      93     7  0.04155  0.08707
   435   0  1.234  5188618    0    170     104     3  0.04047  0.10104
   436   0  1.232  5188741    0    166     101     3  0.04088  0.09288
   437   0  1.232  5189336    0    170      95     1  0.04069  0.09291
   438   0  1.233  5189317    0    155      90     3  0.03952  0.08781
   439   0  1.230  5189247    0    146      98     4  0.03815  0.09695
   440   0  1.227  5189718    0    154      93     3  0.04204  0.09862

   477 OUT OF    480 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.50
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.14
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.405    32.096    34.488  88.624 108.491 111.872
 REC. CELL PARAMETERS   0.041688  0.033767  0.030751  83.859  70.596  67.363
 COORDINATES OF UNIT CELL A-AXIS   -15.223   -14.216    17.810
 COORDINATES OF UNIT CELL B-AXIS    12.089    26.641    13.200
 COORDINATES OF UNIT CELL C-AXIS   -20.942    19.856   -18.883
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.058    -0.002    -0.065
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.045
 LAB COORDINATES OF ROTATION AXIS -0.999862  0.016385 -0.002675
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000709  0.005082  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.31   1244.12
 DETECTOR ORIGIN (PIXELS) AT                     1265.38   1237.44
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.22
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     477 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   52     30     63     67     56     30     56     63     52
 SIGMAB (degree)       0.034  0.034  0.035  0.035  0.035  0.036  0.035  0.034  0.032
 SIGMAR (degree)       0.036  0.032  0.034  0.035  0.037  0.041  0.038  0.036  0.033


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    477 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  3  4  3  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  0  3 11 15 10  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  3  1  0  0  0  0     0  0  0  1  5 14 19 14  7  3  1  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  0  3 10 13  9  4  2  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  2  4  4  2  1  0  0  0  0
     0  0  0  0  1  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  1  2  2  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  2  0  0  0  0  0     0  0  0  0  2  4  4  3  1  0  0  0  0     0  0  0  0  1  3  2  1  0  0  0  0  0
     0  0  0  1  4 11 15  9  4  1  0  0  0     0  0  0  2 10 20 23 13  5  2  1  0  0     0  0  0  2  7 15 16  9  3  1  0  0  0
     0  0  0  2 12 38 50 33 16  6  2  0  0     0  0  0  5 25 66 74 46 22  8  3  0  0     0  0  0  5 20 48 48 30 13  5  1  0  0
     0  0  0  3 19 51 65 45 21  9  3  1  0     0  0  0  5 36 85100 66 30 13  4  1  0     0  0  0  5 27 60 65 41 18  8  3  0  0
     0  0  0  2 11 34 42 29 13  5  2  0  0     0  0  0  4 22 56 64 43 18  8  3  0  0     0  0  0  4 17 38 42 27 12  5  2  0  0
     0  0  0  1  5 11 13  7  3  1  0  0  0     0  0  0  2  8 17 20 11  6  2  1  0  0     0  0  0  1  6 12 14  8  4  1  0  0  0
     0  0  0  1  3  3  3  1  0  0  0  0  0     0  0  0  1  2  3  3  2  1  0  0 -1  0     0  0  0  0  1  2  2  1  1  0  0  0  0
     0  0  0  1  1  1  0  0  0 -1 -1 -1  0     0  0  0  0  0  0  0  0  0 -1 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0 -1  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  4  4  2  1  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 14 13  7  3  1  0  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  8 17 16 10  5  2  1  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5 10 10  6  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  3  4  2  1  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    652 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      441 ...     460
 ******************************************************************************

 USING  20 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ   SIGMAB   SIGMAR
   441   0  1.231  5190631    0    151      79     0  0.04445  0.10508
   442   0  1.226  5190121    0    140      87     2  0.04027  0.09402
   443   0  1.223  5189713    0    139      92     1  0.03874  0.10312
   444   0  1.229  5190399    0    145      82     0  0.03942  0.09363
   445   0  1.227  5190147    0    154      86     0  0.03885  0.09344
   446   0  1.226  5189400    0    156      91     2  0.04034  0.08921
   447   0  1.226  5188458    0    153     104     3  0.03957  0.08910
   448   0  1.224  5189346    0    161      94     1  0.04202  0.11389
   449   0  1.222  5189663    0    149      93     3  0.04305  0.10479
   450   0  1.220  5190401    0    154      88     1  0.03870  0.11524
   451   0  1.216  5189564    0    155      97     4  0.03976  0.09870
   452   0  1.220  5189219    0    151      96     5  0.04164  0.09319
   453   0  1.214  5188861    0    155     100     3  0.04033  0.10580
   454   0  1.217  5188902    0    154     106     1  0.04026  0.11846
   455   0  1.219  5189087    0    153      99     1  0.04292  0.12160
   456   0  1.214  5189565    0    157      91     4  0.04138  0.08730
   457   0  1.211  5189241    0    159      94     1  0.04185  0.11585
   458   0  1.212  5188961    0    156      97     1  0.04025  0.10349
   459   0  1.209  5189286    0    160      95     2  0.04056  0.11052
   460   0  1.208  5189256    0    152      95     2  0.04106  0.12184

   484 OUT OF    488 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.55
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.11
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.412    32.099    34.484  88.622 108.491 111.839
 REC. CELL PARAMETERS   0.041666  0.033754  0.030753  83.875  70.602  67.399
 COORDINATES OF UNIT CELL A-AXIS   -15.223   -14.209    17.827
 COORDINATES OF UNIT CELL B-AXIS    12.084    26.649    13.198
 COORDINATES OF UNIT CELL C-AXIS   -20.944    19.847   -18.882
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.071    -0.009    -0.077
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.045
 LAB COORDINATES OF ROTATION AXIS -0.999862  0.016385 -0.002675
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000347  0.005207  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.32   1244.07
 DETECTOR ORIGIN (PIXELS) AT                     1265.86   1237.23
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.24
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     484 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   53     36     61     65     52     29     52     67     61
 SIGMAB (degree)       0.034  0.035  0.035  0.035  0.034  0.035  0.035  0.034  0.033
 SIGMAR (degree)       0.039  0.034  0.034  0.036  0.040  0.043  0.040  0.037  0.033


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    484 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  0  3  8 12  9  4  2  0  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0    -1  0  0  0  1  2  2  2  1  1  0  0  0     0  0  0  0  4 12 18 13  6  3  1  0  0
    -1  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  0  3  7 11  8  4  2  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  3  4  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0 -1  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  1  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0
     0  0 -1  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  3  2  1  0  0  0 -1  0     0  0  0  0  2  3  3  1  0  0  0  0  0
     0  0  0  1  3  9 11  7  3  1  0  0  0     0  0  0  2  8 18 19 11  5  2  0  0  0     0  0  0  2  7 15 16  9  4  1  0  0  0
     0  0  0  2 12 33 41 29 14  5  2  0  0     0  0  0  4 25 58 66 44 20  8  3  0  0     0  0  0  4 21 46 48 31 14  5  2  0  0
     0  0  0  2 17 48 62 44 21  8  3  1  0     0  0  0  6 34 84100 66 32 13  5  1  0     0  0  1  6 28 64 72 45 22  8  3  0  0
     0  0  0  2 10 29 38 26 12  5  2  0  0     0  0  0  4 20 52 61 41 20  9  3  0  0     0  0  0  4 17 43 46 29 14  6  2  0  0
     0  0  0  1  5 10 12  8  3  1  0  0  0     0  0  0  2  8 18 21 14  6  2  1  0  0     0  0  1  2  7 14 16 10  5  2  0  0  0
     0  0  0  1  2  3  3  1  0  0  0  0  0     0  0  0  1  3  5  4  2  1  0  0 -1  0     0  0  0  0  2  3  3  2  1  0  0  0  0
     0  0  0  1  1  1  0  0  0  0  0 -1  0     0  0  0  0  1  1  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0 -1  0 -1  0  0  0     0  0  0  0  0  0  0 -1  0 -1  0  0  0     0  0  0  0  0 -1  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  1  3  6  6  3  2  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  2  9 20 20 12  5  2  0  0  0     0  0  0  0  1  3  3  2  1  0  0  0  0     0  0  0  0  0  1  0  0  1  0  0  0  0
     0  0  0  3 13 29 29 18  8  3  1  0 -1     0  0  0  0  2  5  4  2  1  0  0  0  0     0  0  0  0  0  1  0  1  0  0  0  0  0
     0  0  0  2  8 19 20 11  5  2  0  0  0     0  0  0  0  2  3  3  1  1  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  3  6  6  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    653 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      461 ...     480
 ******************************************************************************

 USING  20 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ   SIGMAB   SIGMAR
   461   0  1.213  5190444    0    145      83     0  0.04231  0.14218
   462   0  1.215  5190449    0    149      83     1  0.04009  0.11229
   463   0  1.210  5190337    0    150      85     1  0.04406  0.11063
   464   0  1.211  5189359    0    139      96     3  0.04059  0.08373
   465   0  1.208  5190061    0    146      84     2  0.03972  0.09303
   466   0  1.207  5189853    0    150      88     3  0.04099  0.09271
   467   0  1.208  5189809    0    147      91     5  0.04074  0.08937
   468   0  1.206  5189145    0    153      98     2  0.04202  0.12441
   469   0  1.206  5188897    0    153     108     6  0.03771  0.09807
   470   0  1.204  5189201    0    162      98     2  0.04171  0.10963
   471   0  1.205  5189490    0    150      91     3  0.04345  0.09846
   472   0  1.207  5188895    0    160      94     2  0.04077  0.09781
   473   0  1.205  5188615    0    165     100     7  0.04043  0.09646
   474   0  1.200  5188608    0    160     106     3  0.04109  0.11505
   475   0  1.196  5188182    0    159     115     4  0.04167  0.11093
   476   0  1.191  5188508    0    159     102     6  0.04110  0.08926
   477   0  1.190  5187708    0    160     108     1  0.04170  0.09608
   478   0  1.190  5188500    0    151      98     2  0.04212  0.08631
   479   0  1.192  5189795    0    154      82     4  0.04168  0.07963
   480   0  1.192  5190124    0    150      80     2  0.04200  0.09973

   487 OUT OF    490 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.58
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.08
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.417    32.108    34.496  88.611 108.484 111.868
 REC. CELL PARAMETERS   0.041665  0.033752  0.030741  83.879  70.610  67.372
 COORDINATES OF UNIT CELL A-AXIS   -15.221   -14.221    17.826
 COORDINATES OF UNIT CELL B-AXIS    12.084    26.662    13.194
 COORDINATES OF UNIT CELL C-AXIS   -20.961    19.851   -18.883
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.081    -0.022    -0.074
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.045
 LAB COORDINATES OF ROTATION AXIS -0.999862  0.016385 -0.002675
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000525  0.005134  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.27   1244.03
 DETECTOR ORIGIN (PIXELS) AT                     1265.58   1237.28
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.36
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     487 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   51     31     63     72     57     32     53     66     52
 SIGMAB (degree)       0.034  0.035  0.036  0.035  0.035  0.036  0.036  0.034  0.033
 SIGMAR (degree)       0.038  0.032  0.034  0.035  0.038  0.042  0.039  0.037  0.034


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    487 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  3  4  3  2  0  0  0  0
    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  1  5 13 18 13  6  2  1  0 -1
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  3  4  3  2  1  0  0  0     0  0  0  1  6 18 26 19  9  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  1  4 12 15 12  6  2  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  1  2  5  5  4  1  1  0  0  0
     0  0  0  0  1  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  1  2  2  2  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  1  1  1  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0 -1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0 -1 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0 -1  0     0  0  0  0  2  4  4  2  1  0  0 -1  0     0  0  0  0  1  3  2  1  0  0  0 -1  0
     0  0  0  0  5 13 15  9  4  1  0  0  0     0  0  0  1  9 20 21 12  5  2  0  0  0     0  0  0  1  7 15 14  8  3  1  0 -1  0
     0  0  0  2 17 44 54 37 18  7  2  0  0     0  0  0  4 27 62 69 45 21  8  3  0  0     0  0  0  4 19 40 41 27 12  5  2  0  0
     0  0  0  4 24 60 79 55 25 11  4  1  0     0  0  0  7 38 83100 66 29 13  4  1  0     0  0  0  6 27 53 58 37 16  7  2  0  0
     0  0  0  2 15 39 49 34 16  7  2  0  0     0  0  0  5 23 56 65 43 20  8  3  0  0     0  0  0  4 16 36 40 25 12  4  1  0  0
     0  0  0  2  7 15 15 10  5  2  0  0  0     0  0  0  3 10 20 22 14  7  3  0  0  0     0  0  0  2  6 13 14  9  4  1  0  0  0
     0  0  1  1  3  4  3  2  1  0  0 -1  0     0  0  0  1  3  5  4  3  1  0  0  0  0     0  0  0  1  1  3  3  2  1  0  0  0  0
     0  0  0  1  1  1  0  0  0  0  0  0  0     0  0  0  1  1  1  0  0  0 -1 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0 -1  0  0  0 -1  0  0  0     0  0  0  0  0 -1  0 -1  0 -1  0  0  0     0  0  0  0  0 -1  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  5  4  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  6 12 11  7  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  8 14 14  8  4  2  0  0  0     0  0  0  1  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5 10 10  6  3  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  1     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  4  4  3  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    649 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      481 ...     500
 ******************************************************************************

 USING  20 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ   SIGMAB   SIGMAR
   481   0  1.188  5189903    0    145      87     4  0.03921  0.08330
   482   0  1.186  5189380    0    150      93     3  0.04169  0.08380
   483   0  1.185  5189195    0    150      99     3  0.04209  0.09318
   484   0  1.185  5189242    0    158      96     5  0.03859  0.08870
   485   0  1.184  5189002    0    158      97     3  0.04215  0.10535
   486   0  1.185  5188916    0    160     100     4  0.04377  0.09227
   487   0  1.182  5188893    0    161      95     2  0.03941  0.08703
   488   0  1.180  5189705    0    153      91     5  0.03998  0.07889
   489   0  1.176  5189561    0    165      89     1  0.04198  0.09031
   490   0  1.181  5188684    0    152     101     3  0.04059  0.09784
   491   0  1.180  5188991    0    146      90     4  0.03932  0.06864
   492   0  1.180  5189265    0    148      95     2  0.04264  0.11136
   493   0  1.177  5189263    0    148      94     2  0.04219  0.11236
   494   0  1.174  5189741    0    151      88     3  0.04302  0.10878
   495   0  1.175  5189904    0    158      89     5  0.03988  0.08467
   496   0  1.177  5188996    0    144      97     1  0.04020  0.10811
   497   0  1.172  5189125    0    143      96     4  0.04498  0.10122
   498   0  1.170  5189649    0    134      90     2  0.04153  0.10900
   499   0  1.170  5190162    0    137      84     1  0.04100  0.12311
   500   0  1.172  5189930    0    138      85     0  0.04386  0.12249

   487 OUT OF    492 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.55
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.15
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.407    32.096    34.493  88.608 108.486 111.928
 REC. CELL PARAMETERS   0.041700  0.033780  0.030745  83.859  70.603  67.310
 COORDINATES OF UNIT CELL A-AXIS   -15.214   -14.226    17.812
 COORDINATES OF UNIT CELL B-AXIS    12.080    26.660    13.169
 COORDINATES OF UNIT CELL C-AXIS   -20.963    19.838   -18.889
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.103    -0.024    -0.076
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.045
 LAB COORDINATES OF ROTATION AXIS -0.999862  0.016385 -0.002675
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000966  0.005230  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.27   1244.09
 DETECTOR ORIGIN (PIXELS) AT                     1265.00   1237.22
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.30
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     487 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   49     35     61     66     54     30     55     66     58
 SIGMAB (degree)       0.034  0.035  0.035  0.035  0.035  0.036  0.035  0.034  0.033
 SIGMAR (degree)       0.038  0.033  0.035  0.036  0.037  0.039  0.037  0.037  0.035


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    487 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  3  7  8  6  3  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  2  5  7  5  2  1  0  0  0     0  0  0  1  8 21 31 24 10  4  1  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  2  6  9  7  3  1  0  0  0     0  0  0  1 10 28 39 29 14  6  2  0  0
     0  0  0  0  0  2  3  1  0  0  0  0  0     0  0  0  0  2  6  6  4  2  1  0  0  0     0  0  0  1  7 21 24 16  8  3  1  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  2  3  3  2  1  0  0  0  0     0  0  0  1  4  8  8  5  2  1  0  0  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  1  2  1  1  0  0  0  0  0  0     0  0  0  2  3  3  2  1  0  0  0  0  0
     0  0  0  0  1  0  0  0  0  0  0  0  0     0  0  0  1  1  1  0  0  0  0  0 -1  0     0  0  1  1  2  1  0  0  0  0 -1 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0


     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  3  3  2  1  0  0 -1  0     0  0  0  0  2  3  3  1  1  0  0 -1  0     0  0  0  0  1  2  2  1  0  0  0  0  0
     0  0  0  1  8 17 18 12  6  2  0  0  0     0  0  0  2 11 20 19 12  5  2  0  0  0     0  0  0  2  7 11  9  5  2  1  0  0  0
     0  0  0  3 21 53 71 51 22  9  3  1  0     0  0  0  6 31 66 75 50 21  8  3  0  0     0  0  0  5 19 34 35 20  8  3  1  0  0
     0  0  0  4 27 74 92 64 31 13  5  1  0     0  0  1  8 41 93100 64 30 12  4  1  0     0  0  0  5 24 45 45 27 12  5  1  0  0
     0  0  0  3 19 49 57 37 17  7  2  0 -1     0  0  1  6 27 62 65 39 17  7  2  0 -1     0  0  1  4 15 32 31 18  8  3  1  0  0
     0  0  0  2  8 16 17 11  5  2  0  0  0     0  0  1  2 10 19 20 13  6  2  0  0  0     0  0  0  1  5 11 12  8  4  2  0  0  0
     0  0  0  2  3  3  3  2  1  0  0  0  0     0  0  0  1  2  3  4  2  1  0  0 -1  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  1  1  1  0  0  0 -1 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0 -1 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0 -1  0  0  0  0  0  0     0  0  0  0  0  0 -1 -1  0  0  0  0  0     0  0  0 -1 -1  0 -1  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  7  6  3  1  0  0  0  0     0  0  0  0  1  1  0  1  1  0  0  0  0     0  0  0  0  1  0  0  1  0  0  0  0  0
     0  0  0  2  5  9  8  4  2  1  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0
     0  0  0  1  3  6  6  4  2  1  0  0 -1     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  1  3  3  2  1  1  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  0  0  1  0  1  0  0  0  0
     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  0  0  0  0  1  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  1  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0 -1  0  0  0  0  0  0  0  0  0     0  0  0 -1  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    634 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
  REFERENCE PROFILES AS FUNCTION OF DETECTOR POSITION
 ******************************************************************************

 PROFILES ARE DETERMINED AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0




 ***** RUN-AVERAGE OF PROFILE #  1 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  3  4  3  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  0  5 12 15 10  5  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  2  3  5  4  3  0  0  0  0     0  0  0  2  7 17 22 16  7  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  4  3  0  0  0  0  0     0  0  0  0  4 12 15 11  5  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  4  5  4  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  0  0  2  2  0  0  0  0  0  0
     0  0  0  0  5 11 12  7  3  0  0  0  0     0  0  0  2  9 18 18 10  4  0  0  0  0     0  0  0  3  7 13 13  7  3  0  0  0  0
     0  0  0  3 18 40 47 31 13  5  0  0  0     0  0  0  6 31 64 68 42 17  6  0  0  0     0  0  0  5 23 44 43 25 10  4  0  0  0
     0  0  0  4 24 59 70 49 21  8  0  0  0     0  0  0  8 42 91100 66 28 10  0  0  0     0  0  0  7 31 61 62 39 16  6  0  0  0
     0  0  0  3 16 40 48 33 14  5  0  0  0     0  0  0  5 27 61 69 45 19  7  0  0  0     0  0  0  4 20 41 43 27 11  4  0  0  0
     0  0  0  0  4 11 15 10  5  0  0  0  0     0  0  0  0  7 17 21 13  6  2  0  0  0     0  0  0  0  6 12 13  8  3  0  0  0  0
     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  0  0  2  2  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  4  5  4  3  0  0  0  0  0     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  8 14 12  7  3  3  0  0  0     0  0  0  0  3  5  4  7  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3 10 19 18 11  4  3  0  0  0     0  0  0  0  4  5  5  5  4  0  0  0  0     0  0  0  0  0  3  3  0  0  0  0  0  0
     0  0  0  3  7 13 13  7  4  2  0  0  0     0  0  0  0  8  6  6  5  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  4  4  3  2  0  0  0  0     0  0  0  0  0  3  4  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  2 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  5  6  5  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  0  0  0  0  0  0     0  0  0  0  3  7 10  8  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  4  7  6  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  2  0  0  0  0  0  0     0  0  0  0  3  3  3  0  0  0  0  0  0
     0  0  0  0  3  8  8  5  3  0  0  0  0     0  0  0  3  9 18 17 10  5  0  0  0  0     0  0  0  3  9 15 12  7  3  0  0  0  0
     0  0  0  3 13 32 36 24 11  4  0  0  0     0  0  0  8 37 71 69 42 19  8  0  0  0     0  0  0  7 28 46 40 23 11  4  0  0  0
     0  0  0  3 19 48 56 40 19  8  0  0  0     0  0  0  9 50 99100 65 31 13  0  0  0     0  0  0  8 33 56 52 33 15  6  0  0  0
     0  0  0  3 11 30 38 28 14  6  0  0  0     0  0  0  5 28 61 67 46 23  9  0  0  0     0  0  0  4 17 31 32 22 11  4  0  0  0
     0  0  0  0  4  9 13 11  6  5  0  0  0     0  0  0  2  8 18 23 17  9  4  0  0  0     0  0  0  0  5  9 12  9  5  2  0  0  0
     0  0  0  0  0  3  4  4  5  0  0  0  0     0  0  0  0  0  4  4  4  3  0  0  0  0     0  0  0  0  0  2  3  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  8  7  6  5  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  5 12 15 11  6  3  0  0  0  0     0  0  0  0 11  8  5  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  5 12 15 12  7  4  2  0  0  0     0  0  0  0  7  4  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4  5  7  7  5  4  0  0  0  0     0  0  0  0  0  3  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  3  4  4  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  3 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  2  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  7  9  6  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  2  2  0  0  0  0  0     0  0  0  0  4 12 17 13  6  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  2  0  0  0  0  0     0  0  0  0  3 10 15 12  7  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  2  4  7  6  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  4  5  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  2  3  3  3  0  0  0  0  0     0  0  0  0  4  6  5  4  0  0  0  0  0
     0  0  0  0  3  8  8  5  2  0  0  0  0     0  0  0  2  9 17 16  9  4  0  0  0  0     0  0  0  3 11 19 16  8  4  0  0  0  0
     0  0  0  2 12 31 35 22 10  4  0  0  0     0  0  0  6 28 58 57 33 15  6  0  0  0     0  0  0  6 26 46 40 22 10  4  0  0  0
     0  0  0  3 21 55 69 49 23 10  0  0  0     0  0  0  8 42 92100 65 30 13  0  0  0     0  0  0  7 30 56 54 33 15  6  0  0  0
     0  0  0  3 16 46 63 48 25 10  0  0  0     0  0  0  5 30 74 90 63 32 14  0  0  0     0  0  0  4 18 39 44 29 14  6  0  0  0
     0  0  0  0  6 18 28 23 13  6  0  0  0     0  0  0  2 11 30 41 32 17  7  0  0  0     0  0  0  2  7 16 21 16  8  4  0  0  0
     0  0  0  0  3  5  7  6  5  0  0  0  0     0  0  0  0  3  7 10  8  4  3  0  0  0     0  0  0  0  2  4  5  4  3  0  0  0  0
     0  0  0  0  0  0  3  3  0  0  0  0  0     0  0  0  0  0  2  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  8  6  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0 12 11  8  7  5  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4 12 18 14  7  4  0  0  0  0     0  0  0  0  0  8  5  4  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4 11 16 14  8  4  2  0  0  0     0  0  0  0  4  4  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  9  9  6  4  3  0  0  0     0  0  0  0  0  3  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  4  5  5  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  4 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  2  4  5  3  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  0  0  0  0  0     0  0  0  0  4 11 15 11  5  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  2  3  4  3  3  0  0  0  0     0  0  0  2  6 18 24 18  9  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  0  4 13 18 14  7  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  2  4  5  4  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  3  3  3  0  0  0  0  0     0  0  0  0  3  7  6  4  3  0  0  0  0     0  0  0  0  4  6  5  4  3  0  0  0  0
     0  0  0  0  6 15 17 10  5  3  0  0  0     0  0  0  3 12 26 27 15  7  4  0  0  0     0  0  0  3 11 21 19 11  5  3  0  0  0
     0  0  0  3 15 41 51 34 16  6  0  0  0     0  0  0  5 28 64 71 44 20  8  0  0  0     0  0  0  5 22 43 43 25 11  4  0  0  0
     0  0  0  4 21 59 77 55 26 11  0  0  0     0  0  0  6 35 85100 66 31 13  0  0  0     0  0  0  5 24 51 55 34 16  6  0  0  0
     0  0  0  3 14 41 56 41 20  8  0  0  0     0  0  0  4 22 58 72 49 24 10  0  0  0     0  0  0  3 14 33 38 25 11  5  0  0  0
     0  0  0  0  4 12 17 13  7  3  0  0  0     0  0  0  0  6 17 23 17  8  3  0  0  0     0  0  0  0  4 10 13  9  4  2  0  0  0
     0  0  0  0  3  4  4  3  3  0  0  0  0     0  0  0  0  2  3  4  3  2  0  0  0  0     0  0  0  0  0  3  3  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  5  5  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  6  8  7  5  3  0  0  0  0     0  0  0  0  0  4  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  8 13 11  6  4  3  0  0  0     0  0  0  0  5  5  4  4  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  7 13 13  8  4  3  0  0  0     0  0  0  0  4  4  5  4  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  4  8  8  5  3  3  0  0  0     0  0  0  0  2  4  5  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  3  3  3  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  5 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  5  5  3  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  4  4  3  0  0  0  0  0     0  0  0  3  6 13 16 11  5  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  4  6  7  5  4  0  0  0  0     0  0  0  3  8 20 24 17  8  5  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  5  6  4  4  0  0  0  0     0  0  0  0  4 12 17 12  6  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  4  5  4  0  0  0  0  0     0  0  0  0  2  4  4  4  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  4  5  4  3  0  0  0  0  0     0  0  0  0  3  7  6  3  0  0  0  0  0     0  0  0  0  3  4  4  3  0  0  0  0  0
     0  0  0  2  7 16 18 11  5  3  0  0  0     0  0  0  3 11 26 26 15  7  3  0  0  0     0  0  0  3  8 16 15  8  4  0  0  0  0
     0  0  0  3 16 40 48 32 14  6  0  0  0     0  0  0  5 26 62 68 42 19  8  0  0  0     0  0  0  4 19 42 42 24 11  4  0  0  0
     0  0  0  3 20 54 67 47 22  9  0  0  0     0  0  0  5 33 84100 66 30 12  0  0  0     0  0  0  5 26 60 66 40 18  7  0  0  0
     0  0  0  2 10 33 46 32 16  7  0  0  0     0  0  0  3 19 57 74 49 23 10  0  0  0     0  0  0  3 17 45 54 33 15  6  0  0  0
     0  0  0  0  3  7 11  9  4  2  0  0  0     0  0  0  0  5 16 23 15  7  3  0  0  0     0  0  0  0  5 15 19 12  5  3  0  0  0
     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  0  2  3  2  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  4  4  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  4  5  5  4  3  0  0  0  0     0  0  0  0  3  7  6  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  6 12 11  7  5  3  0  0  0     0  0  0  0  5  6  6 12  4  0  0  0  0     0  0  0  0  0  3  0  2  0  0  0  0  0
     0  0  0  3  9 17 18 11  6  5  0  0  0     0  0  0  0  5 11 12  7  9  3  0  0  0     0  0  0  0  0  5  5  3  0  0  0  0  0
     0  0  0  3  6 13 16 10  5  4  0  0  0     0  0  0  0  6  8  9  9  7  0  0  0  0     0  0  0  0  0  4  4  2  0  0  0  0  0
     0  0  0  0  4  6  6  5  4  0  0  0  0     0  0  0  0  0  5 10  6  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  6 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  4  4  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  0  4  7  7  4  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  5  5  5  3  2  0  0  0  0     0  0  0  0  7 15 16 10  5  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  4  8  9  6  4  0  0  0  0     0  0  0  3  9 21 23 15  7  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  6 10 10  7  5  0  0  0  0     0  0  0  0  6 17 20 13  6  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0 11 10  7  0  0  0  0  0     0  0  0  0  6 10 10  6  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0 10  6  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  5  5  3  0  0  0  0  0     0  0  0  0  4  7  6  3  0  0  0  0  0     0  0  0  0  3  5  5  3  0  0  0  0  0
     0  0  0  3  9 19 18 10  5  0  0  0  0     0  0  0  4 15 32 30 17  7  3  0  0  0     0  0  0  3 12 26 25 14  6  0  0  0  0
     0  0  0  4 19 45 46 28 12  5  0  0  0     0  0  3  7 35 78 78 45 20  8  0  0  0     0  0  0  7 33 71 66 37 16  6  0  0  0
     0  0  0  4 23 55 59 37 17  7  0  0  0     0  0  2  9 44 98100 60 27 10  0  0  0     0  0  0  9 44 91 87 50 21  8  0  0  0
     0  0  0  4 15 39 42 26 12  5  0  0  0     0  0  0  5 28 68 72 42 18  7  0  0  0     0  0  0  5 29 64 64 35 15  6  0  0  0
     0  0  0  3  7 16 17 10  4  0  0  0  0     0  0  0  3 10 25 26 14  7  3  0  0  0     0  0  0  2 10 24 23 12  6  4  0  0  0
     0  0  0  0  5  5  4  3  0  0  0  0  0     0  0  0  0  3  4  4  3  0  0  0  0  0     0  0  0  0  0  4  4  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  4  4  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  6  9  9  8  4  0  0  0  0     0  0  0  0  6  6  8  8  3  0  0  0  0     0  0  0  0  0  3  0  0  0  0  0  0  0
     0  0  0  4 12 26 23 14  7  7  0  0  0     0  0  0  3  8 12 13 17  7  0  0  0  0     0  0  0  0  5  8  6  5  2  0  0  0  0
     0  0  0  5 18 36 33 19  9  6  0  0  0     0  0  0  5 11 14 15 14 10  8  0  0  0     0  0  0  0  9 16 14  8  3  0  0  0  0
     0  0  0  4 13 25 23 14  6  6  0  0  0     0  0  0  0  7 10 12  9  8  0  0  0  0     0  0  0  0  8 11 10  4  4  0  0  0  0
     0  0  0  0  6 11 10  5  5  0  0  0  0     0  0  0  0  4  9  7  5  4  0  0  0  0     0  0  0  0  0  4  5  3  0  0  0  0  0
     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  7 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  5  5  3  0  0  0  0  0     0  0  0  0  4  5  6  4  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  4  5  4  4  3  0  0  0  0     0  0  0  0  5 13 17 12  6  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  4  5  4  3  0  0  0  0     0  0  0  2  6 16 21 14  7  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  5  5  3  2  0  0  0  0     0  0  0  0  5 14 16 10  5  2  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  9  6  0  0  0  0  0  0     0  0  0  0  7 11 10  6  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0 14  8  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  0  2  3  3  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0
     0  0  0  0  3  9 12  9  4  2  0  0  0     0  0  0  2  5 15 20 14  6  3  0  0  0     0  0  0  0  4 14 17 11  5  3  0  0  0
     0  0  0  3 11 34 47 32 15  7  0  0  0     0  0  0  4 18 57 74 49 23  9  0  0  0     0  0  0  3 15 49 58 37 17  7  0  0  0
     0  0  0  3 17 50 64 44 21  9  0  0  0     0  0  0  4 28 81100 66 31 13  0  0  0     0  0  0  3 25 67 78 49 22  9  0  0  0
     0  0  0  3 13 39 45 28 13  5  0  0  0     0  0  0  3 21 57 65 38 17  7  0  0  0     0  0  0  3 17 43 47 27 12  5  0  0  0
     0  0  0  4  9 21 20 11  5  0  0  0  0     0  0  0  3 11 23 21 12  5  3  0  0  0     0  0  0  3  6 13 12  6  3  0  0  0  0
     0  0  0  0 11  8  6  3  0  0  0  0  0     0  0  0  0  6  5  4  3  0  0  0  0  0     0  0  0  0  3  3  3  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  4  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  7  8  6  4  0  0  0  0     0  0  0  0  0  3  5  5  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  6 18 19 13  7  4  0  0  0     0  0  0  0  5  9 10 10  6  0  0  0  0     0  0  0  0  0  4  4  3  0  0  0  0  0
     0  0  0  2 10 27 29 17  9  5  0  0  0     0  0  0  0  7 10 12  9  8  0  0  0  0     0  0  0  0  0  7  5  4  2  0  0  0  0
     0  0  0  2  7 17 17 10  5  4  0  0  0     0  0  0  0  5  7  8  8  5  0  0  0  0     0  0  0  0  0  4  4  0  0  0  0  0  0
     0  0  0  0  4  5  5  4  3  0  0  0  0     0  0  0  0  2  4  5  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  2  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  8 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  2  0  0  0  0  0     0  0  0  0  2  3  4  3  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  3  4  4  3  0  0  0  0     0  0  0  0  4 11 17 13  7  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  4  5  4  3  0  0  0  0     0  0  0  0  6 18 25 19  9  5  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  0  5 12 14 10  5  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  4  5  5  3  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  5  8  6  3  2  0  0  0     0  0  0  0  3  8 10  7  3  0  0  0  0     0  0  0  0  2  5  6  4  2  0  0  0  0
     0  0  0  2 10 33 47 35 17  7  0  0  0     0  0  0  3 16 47 61 42 19  8  0  0  0     0  0  0  2 11 29 34 22 10  4  0  0  0
     0  0  0  3 18 57 74 53 26 11  0  0  0     0  0  0  5 31 84100 66 31 14  0  0  0     0  0  0  4 22 54 60 38 17  7  0  0  0
     0  0  0  3 15 38 44 28 13  5  0  0  0     0  0  0  5 24 58 62 38 17  7  0  0  0     0  0  0  4 18 39 39 23 10  4  0  0  0
     0  0  0  4  7 12 11  6  3  0  0  0  0     0  0  0  4  9 17 14  8  4  0  0  0  0     0  0  0  4  6 11  9  5  3  0  0  0  0
     0  0  0  0  0  4  3  0  0  0  0  0  0     0  0  0  0  6  4  3  0  0  0  0  0  0     0  0  0  0  0  3  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  4  8  8  6  4  2  0  0  0     0  0  0  0  0  3  4  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  7 16 16 10  5  4  0  0  0     0  0  0  0  4  5  5  4  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  6 12 11  6  4  3  0  0  0     0  0  0  0  4  5  5  4  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  4  5  4  3  0  0  0  0  0     0  0  0  0  0  2  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  9 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  2  3  3  3  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  0  3  9 13 10  5  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  4  4  0  0  0  0  0     0  0  0  0  4 12 17 13  7  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  2  0  0  0  0  0  0     0  0  0  0  3  6  8  5  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  4  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  4  9 11  7  3  0  0  0  0     0  0  0  3  6 15 17 10  5  2  0  0  0     0  0  0  0  4  9  9  5  3  0  0  0  0
     0  0  0  3 14 40 51 36 17  7  0  0  0     0  0  0  6 30 70 78 50 23 10  0  0  0     0  0  0  4 20 40 40 24 11  4  0  0  0
     0  0  0  3 19 54 67 47 23  9  0  0  0     0  0  0  7 39 92100 64 30 12  0  0  0     0  0  0  6 27 54 52 31 14  6  0  0  0
     0  0  0  2  8 23 29 19  9  4  0  0  0     0  0  0  3 17 41 43 26 12  5  0  0  0     0  0  0  3 13 26 24 14  6  3  0  0  0
     0  0  0  0  2  4  4  3  2  0  0  0  0     0  0  0  0  3  7  7  4  2  0  0  0  0     0  0  0  0  3  5  4  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6  9  8  5  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  7 12 10  6  3  0  0  0  0     0  0  0  0  0  3  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  5  6  5  3  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 STANDARD DEVIATION (DEGREES) OF BEAM DIVERGENCE
 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1    0.025  0.027  0.029  0.027  0.025  0.025  0.024  0.023  0.023
     2    0.024  0.026  0.028  0.027  0.025  0.025  0.024  0.023  0.023
     3    0.025  0.026  0.028  0.027  0.025  0.025  0.024  0.023  0.023
     4    0.024  0.026  0.027  0.026  0.025  0.025  0.023  0.023  0.023
     5    0.024  0.025  0.027  0.026  0.025  0.025  0.023  0.023  0.023
     6    0.024  0.026  0.027  0.026  0.025  0.025  0.024  0.023  0.023
     7    0.024  0.025  0.026  0.026  0.025  0.026  0.024  0.023  0.023
     8    0.024  0.025  0.026  0.025  0.025  0.025  0.024  0.023  0.023
     9    0.025  0.025  0.026  0.025  0.025  0.025  0.023  0.023  0.023
    10    0.024  0.025  0.026  0.025  0.024  0.025  0.024  0.023  0.023
    11    0.024  0.025  0.026  0.025  0.025  0.025  0.024  0.023  0.023
    12    0.024  0.025  0.026  0.025  0.024  0.025  0.024  0.023  0.024
    13    0.023  0.025  0.026  0.025  0.024  0.025  0.024  0.023  0.023
    14    0.024  0.025  0.025  0.025  0.024  0.025  0.024  0.023  0.023
    15    0.024  0.025  0.026  0.025  0.025  0.025  0.023  0.023  0.023
    16    0.024  0.025  0.025  0.025  0.024  0.025  0.024  0.023  0.023
    17    0.024  0.025  0.025  0.025  0.024  0.024  0.024  0.024  0.023
    18    0.024  0.024  0.026  0.025  0.025  0.025  0.025  0.023  0.023
    19    0.024  0.024  0.025  0.025  0.024  0.025  0.024  0.023  0.023
    20    0.024  0.024  0.025  0.025  0.024  0.025  0.025  0.024  0.023
    21    0.024  0.025  0.025  0.025  0.024  0.025  0.025  0.024  0.023
    22    0.024  0.024  0.025  0.025  0.024  0.026  0.025  0.024  0.023
    23    0.024  0.025  0.025  0.025  0.025  0.025  0.025  0.024  0.023
    24    0.024  0.025  0.026  0.025  0.024  0.026  0.026  0.024  0.023
    25    0.024  0.025  0.026  0.025  0.025  0.026  0.025  0.024  0.024


 STANDARD DEVIATION (DEGREES) OF REFLECTING RANGE
 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1    0.036  0.031  0.037  0.039  0.036  0.033  0.035  0.036  0.032
     2    0.036  0.033  0.037  0.035  0.032  0.036  0.037  0.034  0.031
     3    0.035  0.031  0.035  0.034  0.035  0.037  0.036  0.035  0.031
     4    0.035  0.034  0.036  0.035  0.034  0.033  0.037  0.036  0.032
     5    0.036  0.034  0.035  0.034  0.033  0.034  0.037  0.037  0.033
     6    0.037  0.030  0.035  0.036  0.041  0.049  0.042  0.036  0.032
     7    0.035  0.030  0.033  0.034  0.035  0.045  0.042  0.036  0.031
     8    0.035  0.030  0.034  0.036  0.033  0.035  0.036  0.035  0.031
     9    0.035  0.027  0.033  0.034  0.034  0.037  0.037  0.035  0.029
    10    0.036  0.031  0.036  0.036  0.035  0.042  0.039  0.035  0.031
    11    0.041  0.030  0.034  0.036  0.043  0.062  0.055  0.041  0.033
    12    0.038  0.031  0.034  0.035  0.036  0.042  0.040  0.040  0.036
    13    0.035  0.028  0.032  0.036  0.037  0.040  0.038  0.037  0.031
    14    0.035  0.028  0.033  0.034  0.035  0.039  0.039  0.037  0.031
    15    0.036  0.028  0.033  0.039  0.038  0.037  0.038  0.037  0.032
    16    0.036  0.029  0.031  0.034  0.036  0.040  0.038  0.037  0.033
    17    0.045  0.030  0.036  0.043  0.064  0.078  0.069  0.042  0.033
    18    0.036  0.031  0.033  0.035  0.036  0.044  0.041  0.036  0.033
    19    0.037  0.032  0.034  0.035  0.035  0.038  0.041  0.040  0.033
    20    0.037  0.030  0.034  0.034  0.038  0.046  0.040  0.036  0.032
    21    0.036  0.030  0.033  0.034  0.036  0.039  0.038  0.039  0.035
    22    0.035  0.031  0.033  0.035  0.036  0.042  0.039  0.036  0.032
    23    0.039  0.032  0.034  0.035  0.043  0.057  0.046  0.037  0.033
    24    0.039  0.032  0.035  0.035  0.037  0.045  0.040  0.039  0.036
    25    0.039  0.030  0.034  0.039  0.043  0.043  0.042  0.042  0.037


 ***** SUGGESTED VALUES FOR INPUT PARAMETERS *****
 BEAM_DIVERGENCE=   0.24530  BEAM_DIVERGENCE_E.S.D.=   0.02453
 REFLECTING_RANGE=  0.25550  REFLECTING_RANGE_E.S.D.=  0.03650



 NUMBER OF STRONG REFLECTIONS USED FOR LEARNING REFERENCE PROFILES

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1       54     36     65     71     56     31     54     75     62
     2       48     33     60     68     54     31     55     65     53
     3       48     29     56     67     55     30     54     60     57
     4       51     30     61     66     52     31     60     70     55
     5       53     33     59     69     53     30     57     79     64
     6       46     33     63     67     52     26     56     64     51
     7       46     32     63     67     49     30     54     66     56
     8       57     31     58     72     59     34     57     70     60
     9       46     35     61     69     60     33     55     56     50
    10       50     31     61     71     52     25     55     70     58
    11       50     34     65     67     59     31     54     62     51
    12       53     31     53     68     58     36     61     72     64
    13       47     36     61     64     43     27     58     67     56
    14       50     33     62     74     65     33     56     62     57
    15       47     28     52     61     44     28     57     71     55
    16       54     37     61     71     64     35     57     65     58
    17       47     30     59     65     47     32     56     67     57
    18       48     30     58     71     62     32     59     66     56
    19       51     33     57     65     54     29     56     66     57
    20       48     29     55     60     43     29     55     68     55
    21       52     33     61     73     63     31     56     65     52
    22       52     30     63     67     56     30     56     63     52
    23       53     36     61     65     52     29     52     67     61
    24       51     31     63     72     57     32     53     66     52
    25       49     35     61     66     54     30     55     66     58
 margin    1252    808   1500   1697   1364    764   1397   1668   1408



 NUMBER OF GRID POINTS IN REFERENCE PROFILES WITH  SIGNAL >   2.% OF MAXIMUM

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1      156    149    182    186    157    148    141    130    118
     2      144    142    178    167    141    154    142    117    111
     3      144    145    172    171    149    157    134    125    120
     4      150    142    169    166    153    143    136    130    113
     5      153    156    164    164    148    147    141    128    119
     6      152    135    159    163    154    158    149    128    114
     7      146    133    150    159    153    171    148    122    111
     8      139    129    157    158    147    150    137    124    109
     9      153    126    148    152    140    155    134    125    107
    10      140    131    154    158    146    160    145    121    109
    11      169    128    146    155    161    174    151    131    110
    12      150    133    151    152    141    166    146    138    121
    13      139    124    148    154    149    166    141    128    111
    14      138    117    146    150    148    166    142    130    109
    15      138    128    147    156    151    156    137    133    115
    16      139    125    138    148    153    167    148    137    111
    17      170    127    148    158    161    183    168    145    116
    18      140    123    145    154    156    179    157    130    111
    19      148    125    142    153    142    158    151    138    114
    20      147    128    140    148    154    183    160    134    117
    21      145    127    144    148    148    169    156    137    122
    22      150    126    145    148    152    187    162    134    117
    23      154    136    142    148    164    179    168    140    117
    24      155    130    150    152    150    190    169    148    125
    25      154    135    148    161    163    178    158    144    130
 margin     200    176    205    219    235    264    237    188    149

 MINIMUM INTEGRATION VOLUME IS  107 GRID POINTS



 CORRELATIONS BETWEEN REFERENCE PROFILES WITH CENTRAL PROFILE #1

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1     1.00   0.95   0.89   0.95   0.97   0.98   0.99   0.98   0.95
     2     1.00   0.95   0.89   0.97   0.98   0.92   0.95   0.98   0.93
     3     1.00   0.96   0.90   0.98   0.98   0.96   0.98   0.97   0.94
     4     1.00   0.83   0.87   0.98   0.98   0.97   0.97   0.98   0.93
     5     1.00   0.93   0.90   0.98   0.98   0.97   0.98   0.97   0.95
     6     1.00   0.95   0.94   0.98   0.92   0.85   0.95   0.97   0.93
     7     1.00   0.97   0.94   0.98   0.97   0.88   0.91   0.97   0.95
     8     1.00   0.97   0.94   0.97   0.98   0.96   0.98   0.97   0.95
     9     1.00   0.96   0.92   0.99   0.98   0.95   0.97   0.96   0.93
    10     1.00   0.96   0.93   0.97   0.98   0.87   0.92   0.97   0.94
    11     1.00   0.93   0.89   0.97   0.95   0.75   0.75   0.93   0.95
    12     1.00   0.95   0.95   0.98   0.95   0.91   0.98   0.92   0.90
    13     1.00   0.95   0.96   0.98   0.98   0.95   0.97   0.98   0.94
    14     1.00   0.97   0.96   0.98   0.98   0.95   0.96   0.95   0.94
    15     1.00   0.95   0.95   0.96   0.96   0.97   0.98   0.98   0.96
    16     1.00   0.97   0.95   0.98   0.98   0.94   0.97   0.96   0.94
    17     1.00   0.91   0.93   0.97   0.88   0.64   0.81   0.97   0.91
    18     1.00   0.95   0.96   0.97   0.97   0.93   0.91   0.98   0.95
    19     1.00   0.86   0.93   0.98   0.96   0.96   0.95   0.97   0.96
    20     1.00   0.97   0.97   0.98   0.98   0.89   0.93   0.97   0.97
    21     1.00   0.95   0.96   0.99   0.99   0.97   0.96   0.96   0.94
    22     1.00   0.98   0.97   0.98   0.99   0.93   0.91   0.97   0.96
    23     1.00   0.96   0.95   0.98   0.96   0.79   0.86   0.97   0.96
    24     1.00   0.97   0.96   0.97   0.97   0.90   0.91   0.96   0.92
    25     1.00   0.93   0.94   0.97   0.96   0.96   0.93   0.94   0.95
 margin    1.00   0.98   0.96   0.99   0.99   0.96   0.97   0.98   0.96



 CORRELATIONS BETWEEN REFERENCE PROFILES WITH RUN-AVERAGED PROFILES

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1     0.99   0.97   0.96   0.97   0.95   0.94   0.97   0.99   1.00
     2     1.00   0.96   0.96   0.99   0.99   0.98   0.97   0.99   0.99
     3     0.99   0.98   0.96   0.97   0.97   0.99   0.99   1.00   0.99
     4     0.99   0.92   0.96   0.99   0.96   0.86   0.90   0.99   0.98
     5     0.99   0.97   0.97   0.99   0.99   0.94   0.95   0.98   1.00
     6     0.97   0.99   0.97   0.92   0.85   0.87   0.93   0.97   0.99
     7     0.99   0.99   0.99   1.00   0.99   0.95   0.94   0.95   0.96
     8     1.00   0.99   0.98   0.99   0.99   0.98   0.98   0.99   0.99
     9     1.00   0.99   0.99   0.99   0.99   0.99   0.99   0.99   0.98
    10     1.00   0.95   0.96   0.99   1.00   0.98   0.98   0.99   0.99
    11     0.98   0.98   0.99   0.99   0.92   0.81   0.77   0.88   0.95
    12     0.99   0.94   1.00   1.00   0.99   0.99   0.97   0.93   0.91
    13     1.00   0.99   1.00   0.99   0.98   0.99   0.99   1.00   0.99
    14     1.00   0.99   1.00   1.00   1.00   0.99   0.99   0.99   0.99
    15     0.98   0.99   0.98   0.94   0.90   0.91   0.96   0.99   0.99
    16     1.00   0.98   0.99   1.00   1.00   0.99   0.99   0.98   0.97
    17     0.88   0.98   0.91   0.82   0.63   0.51   0.65   0.83   0.95
    18     0.99   0.96   0.99   1.00   0.99   0.98   0.95   0.96   0.95
    19     0.99   0.95   0.99   0.99   0.92   0.88   0.90   0.96   0.99
    20     1.00   0.99   0.99   0.99   0.99   0.97   0.98   1.00   1.00
    21     0.99   0.99   0.99   0.99   0.99   0.98   0.98   0.99   0.99
    22     0.99   0.99   0.99   0.99   1.00   0.99   0.96   0.98   0.98
    23     0.98   0.98   0.99   0.99   0.93   0.85   0.89   0.97   1.00
    24     0.99   0.96   0.98   1.00   1.00   0.98   0.98   0.99   0.97
    25     0.97   0.98   0.98   0.96   0.89   0.86   0.86   0.91   0.95
 margin    1.00   1.00   1.00   1.00   1.00   1.00   1.00   1.00   1.00



      CORRELATION BETWEEN OBSERVED AND EXPECTED PROFILES

  I/Sigma  = mean of intensity/Sigma, where
  Sigma^2  =  4.0*[variance(I; from counting statistics)
                   + 0.00010*I^2]
  CORR     = mean correlation factor between observed and
             expected reflection profiles.
  E.S.D.   = estimated standard deviation of CORR.
  <I>      = mean LP-corrected reflection intensity,
             assuming unpolarized incident beam.
  NUMBER   = number of accepted reflections used to calculate
             I/Sigma, CORR, E.S.D., and <I>.

       I/Sigma        CORR      E.S.D.      <I>      NUMBER

       -1... 0      -0.120      0.010        -1.          2
        0... 1       0.201      0.085         6.         60
        1... 2       0.444      0.109        20.        102
        2... 3       0.597      0.084        41.        134
        3... 6       0.719      0.087       101.        617
        6... 9       0.783      0.075       235.        749
        9...12       0.797      0.071       424.        972
       12...15       0.806      0.075       683.       1046
       15...18       0.802      0.072       994.        989
       18...21       0.806      0.076      1378.        972
       21...24       0.813      0.075      1919.        893
       24...27       0.806      0.077      2413.        806
       27...30       0.808      0.078      3166.        709
       30...33       0.802      0.082      4094.        637
       33...36       0.799      0.085      5313.        632
       36...39       0.792      0.080      7018.        585
       39...42       0.783      0.089      9325.        639
       42...45       0.791      0.078     14120.        680
       45...48       0.783      0.088     25301.        978
       48...51       0.783      0.100     62328.        871

   13073 REFLECTIONS SAVED ON FILE "INTEGRATE.HKL"


 total elapsed wall-clock time for INTEGRATE      16.4 sec

 TASK     cpu time (sec)    elapsed wall-clock time (sec)
   1         264.1                    16.3
