
 ***** INTEGRATE ***** (VERSION Mar 15, 2019  BUILT=20190806)   6-Jan-2020

 INPUT PARAMETER VALUES
 ----------------------
 MAXIMUM_NUMBER_OF_PROCESSORS=   0  (0: automatic choice of OpenMP threads)
 number of OpenMP threads used  20
 MAXIMUM_NUMBER_OF_JOBS=   0  (0:automatic choice of forked main tasks)
 number of forked INTEGRATE tasks   1
 CLUSTER_NODES=
 SECONDS=     0
 NAME_TEMPLATE_OF_DATA_FRAMES=~/data/lysozyme/nitrate/ubatch4/lys_nitr_7_1_????.cbf CBF
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
 UNIT_CELL_CONSTANTS=    27.43    32.15    34.52  88.641 108.451 111.892
 ROTATION_AXIS= -1.00000 -0.00135  0.00038
 OSCILLATION_RANGE=  0.100000 DEGREES
 STARTING_ANGLE=   5.000 STARTING_FRAME=       1
 X-RAY_WAVELENGTH=   0.97700 ANGSTROM
 INCIDENT_BEAM_DIRECTION=  0.000773  0.005098  1.023528
 NUMBER OF DETECTOR SEGMENTS   1
 NX=  2463 NY=  2527 QX= 0.17200 QY= 0.17200
 ORGX=   1265.44  ORGY=   1237.37
 DETECTOR_DISTANCE=   231.314
 DIRECTION_OF_DETECTOR_X-AXIS=   1.00000   0.00000   0.00000
 DIRECTION_OF_DETECTOR_Y-AXIS=   0.00000   1.00000   0.00000
 SEGMENT=    1 2463    1 2527  ! x1,x2,y1,y2 (pixel units)
    SEGMENT_ORGX=      0.00   SEGMENT_ORGY=      0.00
    SEGMENT_DISTANCE=     0.00
    DIRECTION_OF_SEGMENT_X-AXIS=  1.000000  0.000000  0.000000
    DIRECTION_OF_SEGMENT_Y-AXIS=  0.000000  1.000000  0.000000
 RESOLUTION RANGE RECORDED BY DETECTOR (ANGSTROM)    44.608     1.212
 NUMBER OF TRUSTED DETECTOR PIXELS    5213554
 MEAN CONTENTS OF TRUSTED PIXELS IN BACKGROUND TABLE     84.221
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
     1   0  1.058  5208029    0    138      64     3  0.03762  0.05805
     2   0  1.029  5207336    0    149      71     2  0.03899  0.07962
     3   0  1.029  5206481    0    142      77     2  0.03679  0.04867
     4   0  1.030  5207154    0    138      71     4  0.03886  0.04996
     5   0  1.030  5208200    0    132      62     2  0.03892  0.08063
     6   0  1.030  5207979    0    113      64     2  0.03843  0.08832
     7   0  1.036  5208751    0    102      56     1  0.04103  0.12097
     8   0  1.034  5209291    0    103      48     2  0.04396  0.08395
     9   0  1.035  5208989    0     99      48     2  0.04292  0.06521
    10   0  1.035  5208341    0    102      60     0  0.04053  0.07970
    11   0  1.032  5208217    0    114      59     0  0.04411  0.08115
    12   0  1.034  5208002    0    111      60     4  0.04178  0.07364
    13   0  1.035  5207734    0    114      64     1  0.03724  0.08563
    14   0  1.038  5207587    0    121      67     1  0.03786  0.08255
    15   0  1.040  5207785    0    116      66     2  0.04054  0.07145
    16   0  1.038  5208257    0    129      59     6  0.03885  0.05600
    17   0  1.040  5208202    0    121      62     3  0.04104  0.06264
    18   0  1.042  5207584    0    131      69     3  0.04258  0.06618
    19   0  1.040  5207217    0    118      68     5  0.03734  0.05536
    20   0  1.041  5208882    0    113      55     5  0.04247  0.07488

   475 OUT OF    476 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.47
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.07
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.425    32.135    34.492  88.600 108.443 111.913
 REC. CELL PARAMETERS   0.041655  0.033733  0.030737  83.892  70.654  67.333
 COORDINATES OF UNIT CELL A-AXIS   -22.756   -11.159   -10.475
 COORDINATES OF UNIT CELL B-AXIS    -6.410    26.325    17.279
 COORDINATES OF UNIT CELL C-AXIS    14.265    20.255   -23.999
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.000     0.017     0.002
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.034
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001348  0.000378
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000628  0.005202  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.30   1244.09
 DETECTOR ORIGIN (PIXELS) AT                     1265.48   1237.26
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.24
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     475 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   51     29     59     71     59     31     45     65     57
 SIGMAB (degree)       0.032  0.032  0.031  0.031  0.031  0.032  0.031  0.032  0.031
 SIGMAR (degree)       0.026  0.024  0.026  0.027  0.025  0.021  0.025  0.027  0.026


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    475 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  1  4  9 12  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  1  6 15 20 16  9  4  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  4 10 12  8  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  7  9  6  3  2  0  0  0     0  0  0  1  7 13 14  9  4  2  1  0  0     0  0  0  1  5  7  7  5  2  1  0  0  0
     0  0  0  3 15 33 43 31 16  7  3  1  0     0  0  1  7 31 59 65 43 22 10  4  1  0     0  0  1  5 20 34 34 22 10  4  2  0  0
     0  0  0  5 24 55 69 52 28 12  5  2  0     0  0  1 11 47 92100 71 37 15  7  2  0     0  0  1  8 31 53 52 35 18  7  3  1  0
     0  0  0  4 15 36 42 29 16  7  3  1  0     0  0  1  7 28 60 63 42 22  9  4  1  0     0  0  1  5 18 34 34 22 11  4  2  1  0
     0  0  0  1  4  9  9  5  3  1  0  0  0     0  0  0  2  8 16 15  7  4  2  0  0  0     0  0  0  1  5  9  9  4  2  1  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0
     0  0  0  2  4  7  7  4  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  0  1  0  0  0  0  0  0
     0  0  0  2  8 12 11  7  3  1  1  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  0  1  0  0  0  0  0  0
     0  0  0  1  4  8  7  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  1  0  1  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    624 REFLECTION INTENSITIES INTEGRATED


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
    21   0  1.042  5209869    0    113      42     2  0.04154  0.10790
    22   0  1.043  5208429    0    123      57     0  0.04430  0.10847
    23   0  1.042  5206627    0    122      74     6  0.03953  0.05688
    24   0  1.042  5207895    0    126      63     2  0.03932  0.06314
    25   0  1.048  5208148    0    116      59     1  0.04179  0.07742
    26   0  1.047  5208361    0    117      59     0  0.04164  0.08852
    27   0  1.046  5207505    0    120      68     0  0.04338  0.09452
    28   0  1.046  5207705    0    128      65     1  0.04303  0.09958
    29   0  1.047  5207281    0    127      73     7  0.03702  0.05169
    30   0  1.048  5208172    0    133      63     3  0.04386  0.07703
    31   0  1.048  5208649    0    132      55     0  0.04298  0.07684
    32   0  1.047  5207712    0    134      61     3  0.03944  0.05977
    33   0  1.048  5207294    0    133      71     2  0.04046  0.07338
    34   0  1.048  5207037    0    129      73     2  0.04059  0.09008
    35   0  1.047  5208086    0    130      62     4  0.03970  0.08195
    36   0  1.047  5208229    0    120      63     3  0.04218  0.08722
    37   0  1.048  5207763    0    120      65     7  0.04049  0.05315
    38   0  1.049  5208642    0    107      56     0  0.04432  0.09943
    39   0  1.050  5208232    0    112      61     2  0.04037  0.08443
    40   0  1.050  5208187    0    111      60     4  0.03944  0.05978

   477 OUT OF    478 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.50
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.06
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.425    32.140    34.497  88.608 108.439 111.909
 REC. CELL PARAMETERS   0.041653  0.033728  0.030731  83.887  70.656  67.335
 COORDINATES OF UNIT CELL A-AXIS   -22.757   -11.160   -10.473
 COORDINATES OF UNIT CELL B-AXIS    -6.413    26.326    17.284
 COORDINATES OF UNIT CELL C-AXIS    14.263    20.258   -24.005
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.004     0.012     0.001
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.034
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001348  0.000378
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000697  0.005174  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.32   1244.10
 DETECTOR ORIGIN (PIXELS) AT                     1265.40   1237.30
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.26
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     477 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   50     30     56     67     52     29     62     67     57
 SIGMAB (degree)       0.032  0.031  0.031  0.031  0.031  0.031  0.030  0.031  0.031
 SIGMAR (degree)       0.026  0.023  0.026  0.027  0.025  0.023  0.026  0.028  0.026


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    477 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  1  3  9 11  8  5  2  1  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  1  6 14 19 15  8  4  2  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  1  4  9 12  8  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  8  9  6  3  1  0  0  0     0  0  0  2  6 13 15  9  4  2  0  0  0     0  0  0  1  4  8  8  4  2  1  0  0  0
     0  0  0  3 13 34 39 29 15  7  3  1  0     0  0  1  6 27 57 59 41 20  9  4  1  0     0  0  1  5 17 32 31 21 10  4  2  1  0
     0  0  0  4 23 54 67 50 25 11  5  1  0     0  0  1 10 47 90100 69 33 14  6  2  0     0  0  1  8 31 53 54 35 16  7  3  1  0
     0  0  0  3 14 34 41 28 13  6  2  0  0     0  0  1  6 27 57 63 39 18  8  3  1  0     0  0  1  5 17 33 35 21  9  4  2  0  0
     0  0  0  0  3  6  7  4  2  1  0  0  0     0  0  0  1  7 12 11  7  3  1  0  0  0     0  0  0  1  4  7  6  4  2  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  7  7  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  8 12 11  7  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  4  8  8  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    623 REFLECTION INTENSITIES INTEGRATED


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
    41   0  1.053  5208287    0    111      58     1  0.03943  0.07114
    42   0  1.056  5209592    0    111      45     1  0.03799  0.08988
    43   0  1.053  5208790    0     94      54     3  0.03969  0.07819
    44   0  1.051  5208936    0    108      56     4  0.03998  0.07872
    45   0  1.053  5208148    0    120      67     4  0.04225  0.08780
    46   0  1.053  5208071    0    124      64     1  0.04056  0.09097
    47   0  1.050  5208381    0    132      57     1  0.04193  0.06550
    48   0  1.055  5208434    0    113      56     2  0.03819  0.04926
    49   0  1.053  5208528    0    116      57     4  0.03594  0.06335
    50   0  1.055  5208673    0    116      57     1  0.04426  0.09625
    51   0  1.055  5207518    0    121      67     4  0.03902  0.06967
    52   0  1.056  5207597    0    133      71     2  0.04036  0.08596
    53   0  1.057  5207017    0    129      74     3  0.04066  0.07465
    54   0  1.057  5206816    0    134      73     4  0.03926  0.06083
    55   0  1.059  5208050    0    135      66     1  0.04087  0.07450
    56   0  1.059  5207564    0    130      66     1  0.03843  0.06816
    57   0  1.058  5207723    0    131      66     5  0.03818  0.04804
    58   0  1.061  5207839    0    129      67     3  0.03886  0.05279
    59   0  1.058  5208273    0    129      60     4  0.03894  0.04917
    60   0  1.059  5208734    0    130      55     1  0.03881  0.06854

   474 OUT OF    475 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.55
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.04
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.424    32.149    34.514  88.643 108.449 111.906
 REC. CELL PARAMETERS   0.041660  0.033721  0.030721  83.844  70.630  67.325
 COORDINATES OF UNIT CELL A-AXIS   -22.759   -11.165   -10.461
 COORDINATES OF UNIT CELL B-AXIS    -6.414    26.329    17.297
 COORDINATES OF UNIT CELL C-AXIS    14.267    20.259   -24.026
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.005    -0.001    -0.005
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.034
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001348  0.000378
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000976  0.005135  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.43   1244.04
 DETECTOR ORIGIN (PIXELS) AT                     1265.15   1237.30
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.31
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     474 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   46     33     58     64     57     28     50     69     57
 SIGMAB (degree)       0.032  0.032  0.032  0.031  0.031  0.032  0.031  0.032  0.031
 SIGMAR (degree)       0.027  0.024  0.027  0.027  0.026  0.023  0.026  0.028  0.027


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    474 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  3  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  1  4 10 13 10  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  1  6 15 21 15  8  4  1  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  1  4 10 13 10  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  8  9  6  3  1  0  0  0     0  0  0  2  7 14 14  8  4  2  1  0  0     0  0  0  1  5  9  9  4  2  1  0  0  0
     0  0  0  2 16 35 42 30 15  7  3  1  0     0  0  1  6 32 61 62 40 21  9  4  1  0     0  0  1  5 23 42 35 22 10  5  2  0  0
     0  0  0  5 23 53 67 50 27 13  5  2  0     0  0  1 11 47 91100 70 36 16  6  2  0     0  0  1  9 34 60 56 38 18  8  3  1  0
     0  0  1  3 15 35 46 32 15  6  2  1  0     0  0  1  7 28 58 66 44 19  9  3  1  0     0  0  1  5 21 37 38 25 11  5  2  0  0
     0  0  0  1  3  8  9  6  3  1  0  0  0     0  0  0  2  7 14 14  8  4  2  1  0  0     0  0  0  1  6 10  9  5  2  1  0  0  0
     0  0  0  0  0  1  0  1  0  0  0  0  0     0  0  0  0  0  2  1  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  3  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 12  9  5  2  1  0  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3 11 16 13  9  4  2  1  0 -1     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  2  6 10 10  6  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  3  3  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    627 REFLECTION INTENSITIES INTEGRATED


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
    61   0  1.058  5207999    0    133      63     0  0.03745  0.06773
    62   0  1.059  5208045    0    139      68     3  0.03901  0.07180
    63   0  1.057  5208073    0    144      64     3  0.03756  0.06707
    64   0  1.058  5207549    0    142      73     5  0.03764  0.06075
    65   0  1.059  5206995    0    145      77     3  0.03805  0.05340
    66   0  1.058  5207837    0    138      64     2  0.03927  0.05669
    67   0  1.060  5209001    0    124      52     2  0.03905  0.05702
    68   0  1.062  5208779    0    111      54     3  0.04196  0.06503
    69   0  1.061  5208665    0    121      59     6  0.04192  0.06477
    70   0  1.060  5208564    0    126      62     5  0.04232  0.09209
    71   0  1.061  5207910    0    128      65     6  0.04480  0.06846
    72   0  1.064  5207279    0    130      72     7  0.03951  0.04904
    73   0  1.066  5208000    0    139      65     2  0.04155  0.07156
    74   0  1.066  5207425    0    119      71     3  0.03917  0.06205
    75   0  1.068  5208893    0    119      53     0  0.04317  0.09224
    76   0  1.067  5208973    0    117      54     1  0.04671  0.10359
    77   0  1.063  5207887    0    109      68     1  0.04193  0.09605
    78   0  1.067  5207807    0    116      70     7  0.03959  0.07096
    79   0  1.065  5207952    0    117      66     3  0.04335  0.09885
    80   0  1.066  5208025    0    123      66     2  0.04109  0.08743

   472 OUT OF    472 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.51
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.11
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.427    32.147    34.504  88.624 108.450 111.910
 REC. CELL PARAMETERS   0.041656  0.033722  0.030729  83.863  70.637  67.327
 COORDINATES OF UNIT CELL A-AXIS   -22.759   -11.165   -10.469
 COORDINATES OF UNIT CELL B-AXIS    -6.415    26.331    17.290
 COORDINATES OF UNIT CELL C-AXIS    14.268    20.255   -24.014
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.000     0.008     0.001
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.034
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001348  0.000378
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000807  0.005074  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.37   1244.10
 DETECTOR ORIGIN (PIXELS) AT                     1265.31   1237.44
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.30
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     472 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   49     34     57     64     53     28     56     63     57
 SIGMAB (degree)       0.032  0.032  0.032  0.031  0.031  0.031  0.030  0.031  0.031
 SIGMAR (degree)       0.026  0.024  0.026  0.028  0.026  0.023  0.025  0.027  0.026


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    472 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  1  3 10 14 10  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  2  1  1  0  0  0     0  0  0  1  6 16 20 16  9  4  2  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  1  4 11 12 10  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0 -1 -1  0     0  0  0  0  1  1  1  1  0  0 -1 -1  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  1  4  7  9  6  3  1  0  0  0     0  0  0  2  7 12 13  9  4  2  0  0  0     0  0  0  1  4  7  7  5  2  1  0  0  0
     0  0  0  3 14 36 44 29 17  7  3  1  0     0  0  1  6 28 60 65 40 22  9  4  1  0     0  0  1  5 19 37 36 21 11  4  2  1  0
     0  0  0  4 24 60 71 51 28 12  5  1  0     0  0  1  9 47 96100 67 35 15  6  2  0     0  0  1  7 32 58 54 34 17  7  3  1  0
     0  0  0  3 17 38 43 31 15  7  3  1  0     0  0  0  6 29 60 60 41 20  8  4  1  0     0  0  0  4 19 36 32 20 10  4  2  0  0
     0  0  0  0  3  7  8  6  3  1  0  0  0     0  0  0  1  6 11 12  7  4  2  0  0  0     0  0  0  1  4  6  7  3  2  1  0  0  0
     0  0  0  0  0  1  1  1  0  0 -1  0  0     0  0  0  0  1  2  2  2  0  0 -1 -1  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0 -1  0 -1 -1 -1  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0 -1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  9  8  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  8 13 11  7  3  2  1  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  1  5  8  7  4  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  1
     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    633 REFLECTION INTENSITIES INTEGRATED


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
    81   0  1.066  5207221    0    118      68     1  0.04192  0.09152
    82   0  1.063  5207460    0    119      76     2  0.04335  0.09992
    83   0  1.065  5208064    0    113      60     0  0.04778  0.10233
    84   0  1.062  5207879    0    108      64     0  0.04175  0.09520
    85   0  1.064  5207434    0    117      70     2  0.04195  0.08716
    86   0  1.060  5207893    0    121      64     1  0.04308  0.08457
    87   0  1.065  5207496    0    117      67     1  0.04235  0.08443
    88   0  1.066  5207804    0    110      66     1  0.04040  0.08790
    89   0  1.067  5208294    0    119      55     1  0.04261  0.08460
    90   0  1.066  5207496    0    119      66     2  0.04095  0.08786
    91   0  1.068  5207462    0    125      72     7  0.04095  0.06487
    92   0  1.067  5207843    0    123      62     5  0.04104  0.05418
    93   0  1.068  5208521    0    123      58     3  0.04134  0.08576
    94   0  1.072  5208176    0    117      65     3  0.04348  0.10705
    95   0  1.069  5208172    0    111      62     3  0.03959  0.08056
    96   0  1.071  5208478    0    109      62     4  0.04101  0.07012
    97   0  1.070  5208982    0    111      50     2  0.04025  0.07297
    98   0  1.073  5209138    0    130      50     4  0.04180  0.05784
    99   0  1.074  5208397    0    130      58     2  0.03880  0.05474
   100   0  1.072  5207749    0    141      64     3  0.03612  0.04967

   455 OUT OF    456 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.52
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.422    32.142    34.492  88.618 108.448 111.910
 REC. CELL PARAMETERS   0.041661  0.033727  0.030738  83.870  70.641  67.329
 COORDINATES OF UNIT CELL A-AXIS   -22.756   -11.160   -10.471
 COORDINATES OF UNIT CELL B-AXIS    -6.413    26.325    17.291
 COORDINATES OF UNIT CELL C-AXIS    14.264    20.254   -24.000
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.008     0.012    -0.002
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.034
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001348  0.000378
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000832  0.005087  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.31   1244.06
 DETECTOR ORIGIN (PIXELS) AT                     1265.22   1237.38
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.25
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     455 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     28     49     64     55     27     52     66     57
 SIGMAB (degree)       0.032  0.032  0.032  0.032  0.031  0.031  0.030  0.031  0.031
 SIGMAR (degree)       0.027  0.024  0.027  0.028  0.025  0.023  0.025  0.027  0.026


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    455 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  1  4 11 13 10  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  1  7 16 20 15  8  4  2  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  4 10 12  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  3  8  9  6  3  1  0  0  0     0  0  0  1  6 13 14  9  5  1  0  0  0     0  0  0  1  4  8  7  4  2  0  0  0  0
     0  0  0  3 16 39 44 32 15  7  3  1  0     0  0  1  6 30 60 62 41 20  9  4  1  0     0  0  1  5 21 34 33 20 10  4  2  1  0
     0  0  0  5 28 59 71 51 25 12  5  2  0     0  0  1 10 49 93100 68 33 16  6  2  0     0  0  1  8 32 54 54 35 16  8  3  1  0
     0  0  0  4 17 37 41 30 15  6  2  1  0     0  0  1  7 28 58 60 41 20  8  3  1  0     0  0  0  6 18 34 32 21 10  4  2  0  0
     0  0  0  1  3  6  7  5  2  1  0  0  0     0  0  0  2  6 10 12  7  3  1  0  0  0     0  0  0  1  4  7  7  4  2  1  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  1  2  0  0  0 -1 -1  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  7  7  3  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  8 14 12  7  3  1  0  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  2  4  9  8  5  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  1  0  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    623 REFLECTION INTENSITIES INTEGRATED


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
   101   0  1.070  5207188    0    138      73     5  0.03860  0.05299
   102   0  1.070  5207859    0    139      68     4  0.03679  0.05546
   103   0  1.072  5208330    0    135      60     3  0.03742  0.05379
   104   0  1.072  5209099    0    122      48     0  0.03967  0.06442
   105   0  1.073  5208629    0    110      60     1  0.03793  0.08466
   106   0  1.072  5208827    0    114      55     3  0.04057  0.07830
   107   0  1.073  5208513    0    119      56     3  0.04145  0.08012
   108   0  1.073  5208275    0    121      58     6  0.03787  0.04841
   109   0  1.072  5208317    0    132      57     6  0.03866  0.04272
   110   0  1.073  5208661    0    124      59     2  0.04218  0.08666
   111   0  1.072  5207764    0    119      67     2  0.03943  0.07895
   112   0  1.070  5208710    0    120      55     4  0.03892  0.05902
   113   0  1.073  5209018    0    104      46     3  0.04011  0.09055
   114   0  1.073  5209714    0    106      43     1  0.04682  0.13436
   115   0  1.073  5208444    0    106      57     6  0.04231  0.06825
   116   0  1.072  5208441    0    120      58     6  0.03930  0.06113
   117   0  1.074  5207962    0    124      66     2  0.04414  0.09327
   118   0  1.073  5207733    0    130      71     2  0.04195  0.09403
   119   0  1.074  5207953    0    145      67     5  0.04062  0.06393
   120   0  1.073  5206963    0    141      75     3  0.03778  0.07301

   462 OUT OF    463 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.58
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.425    32.145    34.506  88.641 108.419 111.913
 REC. CELL PARAMETERS   0.041653  0.033725  0.030722  83.857  70.663  67.322
 COORDINATES OF UNIT CELL A-AXIS   -22.756   -11.158   -10.478
 COORDINATES OF UNIT CELL B-AXIS    -6.414    26.325    17.297
 COORDINATES OF UNIT CELL C-AXIS    14.265    20.256   -24.017
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.006     0.015    -0.000
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.034
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001348  0.000378
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000753  0.005051  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.29   1244.03
 DETECTOR ORIGIN (PIXELS) AT                     1265.30   1237.39
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.26
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     462 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   45     31     57     65     49     26     52     67     56
 SIGMAB (degree)       0.032  0.032  0.032  0.032  0.031  0.032  0.031  0.031  0.031
 SIGMAR (degree)       0.027  0.028  0.027  0.027  0.025  0.023  0.026  0.028  0.027


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    462 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  1  3  8  9  8  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0     0  0  0  1  5 11 14 11  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  1  3  7  9  8  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0 -1  0     0  0  0  0  0  1  1  0  0  0  0 -1  0
     0  0  0  1  3  8  8  5  3  1  0  0  0     0  0  0  2  6 15 15  9  4  2  0  0  0     0  0  0  2  6 10 10  5  2  1  0  0  0
     0  0  0  3 13 31 34 27 14  6  2  1  0     0  0  1  8 30 62 66 45 23 10  4  1  0     0  0  1  7 26 47 47 29 15  6  2  1  0
     0  0  0  4 19 42 54 41 22 11  4  1  0     0  0  1 11 45 89100 72 36 17  6  2  0     0  0  1 10 40 68 68 46 22 10  4  1  0
     0  0  0  3 12 27 34 25 12  6  2  0  0     0  0  0  7 26 56 64 44 20  9  3  1  0     0  0  1  7 22 43 44 28 13  6  2  1  0
     0  0  0  0  2  5  7  4  2  0  0  0  0     0  0  0  1  6 11 14  7  3  1  0  0  0     0  0  0  2  5 10 11  5  2  1  0  0  0
     0  0  0  0  0  0  1  0  0  0  0 -1  0     0  0  0  0  0  1  2  0  0 -1 -1 -1  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0 -1 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  1  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  3  3  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  3  8 13 13  8  4  1  0  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  1  4 13 20 18 12  6  2  1  0  0     0  0  0  1  2  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  2  8 13 12  7  3  1  0  0  0     0  0  0  1  1  2  2  1  0  0  0  0  0     1  0  0  1  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  3  4  1  1  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    644 REFLECTION INTENSITIES INTEGRATED


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
   121   0  1.071  5206214    0    139      83     1  0.04079  0.07867
   122   0  1.071  5206160    0    146      80     4  0.03953  0.07114
   123   0  1.072  5206485    0    132      78     8  0.04001  0.06359
   124   0  1.074  5207643    0    119      63     4  0.04338  0.07579
   125   0  1.073  5208690    0    115      49     4  0.04050  0.07085
   126   0  1.069  5207803    0    111      61     5  0.04127  0.06197
   127   0  1.072  5206726    0    123      76     3  0.04383  0.06342
   128   0  1.070  5207441    0    125      66     2  0.03824  0.06803
   129   0  1.075  5207361    0    116      65     2  0.03939  0.06690
   130   0  1.073  5207865    0    118      56     2  0.04025  0.05592
   131   0  1.072  5208623    0    114      54     1  0.03958  0.08550
   132   0  1.075  5208407    0    109      56     1  0.04120  0.07350
   133   0  1.074  5208304    0    125      55     1  0.03845  0.06497
   134   0  1.072  5207970    0    130      57     2  0.04189  0.07441
   135   0  1.073  5206932    0    146      73     3  0.03792  0.06734
   136   0  1.074  5206534    0    145      72     5  0.03679  0.04611
   137   0  1.073  5207023    0    141      67     4  0.04045  0.05645
   138   0  1.074  5207868    0    133      59     1  0.03915  0.09753
   139   0  1.075  5208217    0    124      54     4  0.04301  0.08628
   140   0  1.075  5207611    0    122      63     4  0.04005  0.09061

   470 OUT OF    471 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.48
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.07
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.429    32.151    34.490  88.617 108.447 111.913
 REC. CELL PARAMETERS   0.041652  0.033718  0.030740  83.871  70.642  67.327
 COORDINATES OF UNIT CELL A-AXIS   -22.760   -11.161   -10.478
 COORDINATES OF UNIT CELL B-AXIS    -6.415    26.330    17.298
 COORDINATES OF UNIT CELL C-AXIS    14.266    20.256   -23.995
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.013     0.016    -0.002
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.034
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001348  0.000378
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000702  0.004981  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.30   1244.06
 DETECTOR ORIGIN (PIXELS) AT                     1265.38   1237.51
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.30
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     470 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   49     32     56     69     59     28     48     60     58
 SIGMAB (degree)       0.032  0.032  0.032  0.031  0.031  0.032  0.030  0.031  0.031
 SIGMAR (degree)       0.026  0.023  0.027  0.028  0.026  0.022  0.025  0.027  0.026


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    470 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  2  2  1  0  0  0  0     0  0  0  1  4 11 14 11  5  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  3  3  3  2  1  0  0  0     0  0  0  1  7 18 22 18  9  4  2  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  1  4 12 14 11  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  3  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  1  1  0  0  0  0 -1  0     0  0  0  0  1  1  2  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  1  4  8 11  7  3  1  1  0  0     0  0  0  2  6 13 15  9  4  2  1  0  0     0  0  0  1  4  7  8  4  2  1  0  0  0
     0  0  0  3 16 40 48 34 16  8  3  1  0     0  0  1  7 30 62 66 42 20 10  4  1  0     0  0  1  5 17 33 32 20 10  5  2  1  0
     0  0  0  6 29 63 74 56 28 13  5  2  0     0  0  1 13 49 95100 70 34 16  6  2  0     0  0  1  9 28 50 49 33 16  8  3  1  0
     0  0  0  3 16 40 45 33 16  6  3  1  0     0  0  0  6 27 57 61 40 20  8  3  1  0     0  0  0  4 16 31 30 19  9  3  1  0  0
     0  0  0  1  4 10  9  6  2  1  0  0  0     0  0  0  1  6 15 14  8  3  1  0  0  0     0  0  0  1  4  8  7  4  2  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0 -1  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  1  0 -1  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  1  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4  6  6  3  2  1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 10  9  6  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0 -1     0  0  0  0  0  0  0  0  0  0  0  0 -1
     0  0  0  1  4  7  6  3  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0
     0 -1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0 -1  0  0  0  0  0  0

    628 REFLECTION INTENSITIES INTEGRATED


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
   141   0  1.072  5207154    0    129      64     5  0.04118  0.07829
   142   0  1.075  5206900    0    131      70     5  0.04141  0.06550
   143   0  1.079  5207282    0    125      65     4  0.04006  0.07133
   144   0  1.074  5207250    0    117      64     5  0.03919  0.06425
   145   0  1.073  5207838    0    118      60     5  0.04277  0.07267
   146   0  1.073  5208059    0    118      57     1  0.04386  0.09606
   147   0  1.074  5207648    0    123      59     0  0.03785  0.10271
   148   0  1.076  5207933    0    122      57     1  0.04164  0.10331
   149   0  1.072  5207396    0    129      64     2  0.04226  0.10179
   150   0  1.073  5206364    0    133      73     8  0.03761  0.06195
   151   0  1.077  5206886    0    141      70     7  0.03770  0.06767
   152   0  1.073  5207105    0    144      69     4  0.04187  0.07469
   153   0  1.076  5206579    0    142      77     7  0.04215  0.06434
   154   0  1.074  5206258    0    133      78     8  0.04147  0.06692
   155   0  1.076  5207671    0    148      65     4  0.04537  0.06881
   156   0  1.077  5207228    0    140      67     5  0.03955  0.06292
   157   0  1.078  5206651    0    133      74     5  0.03757  0.05321
   158   0  1.078  5208262    0    139      59     5  0.04020  0.06840
   159   0  1.079  5208256    0    128      57     3  0.04202  0.08687
   160   0  1.080  5208230    0    117      60     6  0.03950  0.05127

   495 OUT OF    495 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.48
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.11
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.427    32.149    34.506  88.631 108.438 111.907
 REC. CELL PARAMETERS   0.041651  0.033719  0.030725  83.862  70.648  67.329
 COORDINATES OF UNIT CELL A-AXIS   -22.759   -11.162   -10.472
 COORDINATES OF UNIT CELL B-AXIS    -6.416    26.330    17.296
 COORDINATES OF UNIT CELL C-AXIS    14.266    20.259   -24.016
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.006     0.009    -0.000
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.034
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001348  0.000378
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000794  0.005161  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.36   1244.02
 DETECTOR ORIGIN (PIXELS) AT                     1265.31   1237.24
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.30
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     495 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   52     30     55     71     55     34     61     68     58
 SIGMAB (degree)       0.032  0.032  0.031  0.031  0.031  0.031  0.031  0.031  0.031
 SIGMAR (degree)       0.026  0.024  0.026  0.027  0.026  0.022  0.025  0.027  0.026


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    495 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  3  8 10  8  4  2  1  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  1  5 13 18 13  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  3  9 11  8  4  2  1  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0 -1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0 -1  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  6  8  6  2  1  0  0  0     0  0  0  1  5 11 13  8  3  1  0  0  0     0  0  0  1  4  7  8  5  2  1  0  0  0
     0  0  0  2 13 31 38 28 13  6  3  1  0     0  0  0  5 28 56 63 42 19  9  4  1  0     0  0  0  5 20 35 36 21 10  5  2  0  0
     0  0  0  4 23 52 63 47 24 11  4  1  0     0  0  1  9 46 92100 70 34 16  6  2  0     0  0  1  7 32 59 58 39 19  8  3  1  0
     0  0  0  3 13 34 39 28 15  6  2  1  0     0  0  1  5 26 59 62 41 21  8  3  1  0     0  0  1  4 19 38 37 24 11  4  2  0  0
     0  0  0  0  3  7  7  4  2  1  0  0  0     0  0  0  1  6 13 12  6  4  1  0  0  0     0  0  0  1  5  8  7  4  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  6  9  8  4  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  9 15 13  9  4  2  1  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  9  9  6  2  1  1  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    650 REFLECTION INTENSITIES INTEGRATED


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
   161   0  1.079  5208190    0    114      61     6  0.04051  0.05447
   162   0  1.079  5208308    0    120      61     3  0.03779  0.06560
   163   0  1.079  5208719    0    124      59     1  0.03798  0.07838
   164   0  1.079  5208113    0    112      62     1  0.03789  0.05708
   165   0  1.080  5208261    0    105      61     5  0.04037  0.07519
   166   0  1.081  5209280    0     98      49     2  0.04370  0.10197
   167   0  1.081  5209087    0    100      50     1  0.04110  0.11907
   168   0  1.081  5208612    0    104      53     4  0.03855  0.08519
   169   0  1.084  5208225    0    108      60     4  0.03787  0.07853
   170   0  1.086  5208476    0    110      57     6  0.03996  0.06422
   171   0  1.081  5207965    0    124      65     4  0.04174  0.07973
   172   0  1.085  5207944    0    122      64     6  0.03949  0.06889
   173   0  1.086  5207773    0    119      66     6  0.03952  0.06545
   174   0  1.084  5208894    0    115      56     4  0.03890  0.07867
   175   0  1.082  5208675    0    113      53     2  0.03946  0.07706
   176   0  1.089  5208663    0    118      60     3  0.04471  0.08123
   177   0  1.087  5208470    0    119      59     5  0.04210  0.05670
   178   0  1.085  5208681    0    112      57     3  0.04066  0.07830
   179   0  1.087  5208731    0    128      55     2  0.03808  0.08863
   180   0  1.087  5208209    0    123      63     5  0.04092  0.07286

   439 OUT OF    439 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.54
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.08
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.422    32.151    34.498  88.617 108.486 111.903
 REC. CELL PARAMETERS   0.041669  0.033717  0.030741  83.859  70.602  67.333
 COORDINATES OF UNIT CELL A-AXIS   -22.758   -11.171   -10.453
 COORDINATES OF UNIT CELL B-AXIS    -6.417    26.337    17.287
 COORDINATES OF UNIT CELL C-AXIS    14.267    20.252   -24.009
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.002    -0.005     0.001
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.034
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001348  0.000378
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.001109  0.005162  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.38   1244.12
 DETECTOR ORIGIN (PIXELS) AT                     1264.93   1237.33
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.29
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     439 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   48     32     60     61     44     22     45     63     56
 SIGMAB (degree)       0.032  0.032  0.032  0.031  0.031  0.031  0.030  0.031  0.031
 SIGMAR (degree)       0.029  0.024  0.027  0.028  0.027  0.027  0.027  0.027  0.026


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    439 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  3  4  3  1  0  0  0  0
     0  0  0  0  1  0  1  0  1  1  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  1  5 14 18 12  7  3  1  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  3  5  4  2  1  0  0  0     0  0  0  2 10 24 30 23 12  6  2  1  0
     0  0  0  0  1  1  0  1  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  1  6 15 18 13  7  3  1  0  0
     0  0  0  0  1  1  0  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  2  3  3  3  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0 -1 -1 -1  0     0  0  0  0  0  0  0  0  0 -1 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  1  4 10 11  7  4  1  0  0  0     0  0  0  1  6 14 14  8  4  1  0  0  0     0  0  0  1  4  9  8  4  2  1  0  0  0
     0  0  1  4 19 47 54 36 19  8  3  1  0     0  0  1  7 29 62 63 41 21  8  3  1  0     0  0  1  5 22 40 36 23 11  4  2  0  0
     0  0  1  6 35 76 86 63 31 14  5  2  1     0  0  1 11 49 95100 70 34 15  6  2  1     0  0  1  8 33 58 56 36 17  7  3  1  0
     0  0  0  4 20 47 52 36 18  8  3  1  0     0  0  1  6 26 56 59 41 19  9  3  1  0     0  0  0  4 17 33 32 21 10  4  2  0  0
     0  0  0  1  4  9 10  7  4  1  0  0  0     0  0  0  1  6 11 11  8  3  1  0  0  0     0  0  0  1  4  6  6  5  2  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 10  8  4  2  1  0  0  0     0  0  0  1  1  1  1  0  0  0  0  0  0     0  0  0  0  1  0  1  0  0  0  0  0  0
     0  0  0  2  8 13 12  7  3  1  0  0  0     0  0  0  1  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0
     0  0  0  1  4  8  7  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0 -1 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    605 REFLECTION INTENSITIES INTEGRATED


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
   181   0  1.088  5207754    0    118      70     3  0.03913  0.07605
   182   0  1.083  5208325    0    119      55     3  0.04101  0.07130
   183   0  1.084  5208428    0    122      56     2  0.04024  0.07280
   184   0  1.085  5207683    0    118      61     0  0.03682  0.06847
   185   0  1.084  5208062    0    120      63     0  0.04150  0.07741
   186   0  1.084  5208168    0    126      59     2  0.03944  0.06388
   187   0  1.090  5207640    0    116      67     6  0.04143  0.06749
   188   0  1.087  5207973    0    124      64     2  0.04185  0.10301
   189   0  1.087  5208274    0    127      63     4  0.04160  0.08501
   190   0  1.086  5208434    0    118      59     4  0.04112  0.08177
   191   0  1.087  5208562    0    114      53     6  0.03715  0.05579
   192   0  1.087  5208039    0    114      61     4  0.03854  0.05571
   193   0  1.086  5208548    0    120      58     2  0.04306  0.07756
   194   0  1.085  5208577    0    131      59     6  0.03863  0.07051
   195   0  1.084  5208453    0    125      60     5  0.03790  0.06992
   196   0  1.086  5208047    0    121      63     6  0.03798  0.04994
   197   0  1.087  5208787    0    123      55     3  0.04015  0.06095
   198   0  1.089  5209233    0    110      50     5  0.03622  0.03998
   199   0  1.087  5209694    0    105      47     3  0.03866  0.04650
   200   0  1.091  5210253    0    101      38     2  0.04398  0.07794

   476 OUT OF    476 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.48
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.07
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.422    32.143    34.491  88.617 108.456 111.906
 REC. CELL PARAMETERS   0.041662  0.033724  0.030741  83.871  70.634  67.334
 COORDINATES OF UNIT CELL A-AXIS   -22.756   -11.163   -10.465
 COORDINATES OF UNIT CELL B-AXIS    -6.414    26.330    17.286
 COORDINATES OF UNIT CELL C-AXIS    14.264    20.249   -24.003
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.001     0.006    -0.000
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.034
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001348  0.000378
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000829  0.005101  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.34   1244.11
 DETECTOR ORIGIN (PIXELS) AT                     1265.25   1237.40
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.24
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     476 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   49     28     54     69     56     30     54     66     55
 SIGMAB (degree)       0.032  0.032  0.031  0.032  0.031  0.031  0.031  0.031  0.031
 SIGMAR (degree)       0.026  0.024  0.026  0.027  0.025  0.022  0.024  0.027  0.026


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    476 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  1  4  9 12  9  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0     0  0  0  1  6 15 18 14  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  0  4 10 11  7  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  1  4  7  9  5  3  1  0  0  0     0  0  0  2  7 13 16  8  4  2  0  0  0     0  0  0  1  5  8  9  5  2  1  0  0  0
     0  0  0  4 16 35 44 31 15  6  3  1  0     0  0  1  7 31 62 67 45 21  9  4  1  0     0  0  1  5 21 37 36 24 11  5  2  0  0
     0  0  0  5 26 56 68 49 25 12  4  2  0     0  0  1 10 47 93100 68 33 15  6  2  0     0  0  1  7 31 54 55 36 17  8  3  1  0
     0  0  0  3 17 38 43 27 14  6  3  1  0     0  0  0  6 29 59 64 38 18  8  3  1  0     0  0  0  5 18 34 35 21 10  4  2  0  0
     0  0  0  1  3  7  9  4  3  1  0  0  0     0  0  0  1  6 13 13  7  4  1  0  0  0     0  0  0  1  4  7  8  4  2  1  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  8  7  5  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  8 12 12  7  3  1  1  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  8  8  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    608 REFLECTION INTENSITIES INTEGRATED


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
   201   0  1.089  5209276    0    117      51     5  0.03891  0.05971
   202   0  1.086  5208474    0    118      61     4  0.04170  0.06418
   203   0  1.087  5208154    0    135      64     3  0.03895  0.06949
   204   0  1.082  5207746    0    135      70     7  0.03613  0.05338
   205   0  1.086  5208556    0    138      60     6  0.03768  0.05262
   206   0  1.086  5208714    0    130      58     2  0.03871  0.08550
   207   0  1.086  5209035    0    121      49     4  0.03743  0.05270
   208   0  1.086  5209125    0    117      49     3  0.03773  0.05348
   209   0  1.087  5208328    0    123      57     2  0.04079  0.06802
   210   0  1.087  5207864    0    129      62     2  0.03944  0.07090
   211   0  1.086  5207442    0    131      71     7  0.04012  0.06021
   212   0  1.085  5208012    0    126      65     2  0.03869  0.07799
   213   0  1.086  5207956    0    125      63     7  0.03616  0.05667
   214   0  1.087  5208808    0    130      55     1  0.03957  0.07800
   215   0  1.084  5207209    0    125      73     5  0.03929  0.05991
   216   0  1.084  5207673    0    137      71     7  0.03788  0.05322
   217   0  1.087  5206757    0    133      75     6  0.03894  0.04616
   218   0  1.086  5208225    0    129      63     5  0.03895  0.05817
   219   0  1.088  5208490    0    130      61     2  0.03926  0.06655
   220   0  1.086  5207985    0    112      70     4  0.03771  0.07810

   490 OUT OF    491 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.53
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.06
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.425    32.146    34.492  88.629 108.424 111.906
 REC. CELL PARAMETERS   0.041649  0.033721  0.030734  83.871  70.664  67.333
 COORDINATES OF UNIT CELL A-AXIS   -22.757   -11.157   -10.479
 COORDINATES OF UNIT CELL B-AXIS    -6.417    26.329    17.292
 COORDINATES OF UNIT CELL C-AXIS    14.263    20.250   -24.005
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.003     0.014     0.002
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.034
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001348  0.000378
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000612  0.005085  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.30   1244.03
 DETECTOR ORIGIN (PIXELS) AT                     1265.49   1237.35
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.27
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     490 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   51     33     53     66     60     32     57     71     58
 SIGMAB (degree)       0.032  0.032  0.032  0.032  0.031  0.031  0.031  0.031  0.031
 SIGMAR (degree)       0.026  0.023  0.026  0.027  0.025  0.023  0.026  0.028  0.026


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    490 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  1  3  8 11  8  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  1  6 13 18 13  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  1  4  9 12  8  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  3  6  8  5  3  1  0  0  0     0  0  0  1  7 13 14  8  4  1  0  0  0    -1  0  0  1  6  9  8  5  2  1  0  0  0
     0  0  0  3 14 32 40 27 14  6  2  1  0     0  0  1  7 29 57 65 42 21  9  3  1  0     0  0  1  6 22 38 39 24 12  5  2  0  0
     0  0  0  5 24 50 62 44 22 11  4  2  0     0  0  1 10 49 92100 66 33 16  6  2  0     0  0  1  9 36 61 61 39 19  9  3  1  0
     0  0  0  3 14 34 40 25 14  6  2  1  0     0  0  0  6 29 61 63 39 21  8  3  1  0     0  0  1  6 21 39 37 22 12  5  1  0  0
     0  0  0  1  3  6  9  5  3  1  0  0  0     0  0  0  1  7 11 14  9  5  1  0  0  0     0  0  0  1  5  7  9  4  3  1  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0 -1  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  8  8  5  2  1  1  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  2  9 14 13  9  4  2  1  0  1     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  9  9  5  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    628 REFLECTION INTENSITIES INTEGRATED


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
   221   0  1.087  5208282    0    115      62     1  0.03923  0.12972
   222   0  1.083  5207929    0    120      65     2  0.03914  0.10703
   223   0  1.088  5207852    0    118      67     7  0.04265  0.06842
   224   0  1.088  5208590    0    114      59     0  0.04589  0.12348
   225   0  1.089  5208127    0    108      58     4  0.03737  0.06978
   226   0  1.087  5208340    0    117      59     3  0.04166  0.08206
   227   0  1.087  5207910    0    112      62     7  0.04152  0.05109
   228   0  1.089  5208895    0    112      51     1  0.04018  0.08024
   229   0  1.087  5209276    0    102      49     0  0.04151  0.11034
   230   0  1.089  5209003    0    109      52     4  0.03696  0.07249
   231   0  1.089  5209112    0    121      49     5  0.03956  0.06027
   232   0  1.088  5207918    0    119      63     7  0.03904  0.05896
   233   0  1.089  5208919    0    120      57     6  0.04320  0.08278
   234   0  1.090  5208348    0    126      60     6  0.04099  0.07993
   235   0  1.088  5207272    0    119      70     6  0.04284  0.07234
   236   0  1.086  5207676    0    120      64     0  0.04252  0.10110
   237   0  1.088  5207858    0    125      63     2  0.04562  0.08772
   238   0  1.086  5208074    0    122      59     2  0.04193  0.08254
   239   0  1.087  5207760    0    125      67     1  0.04129  0.08027
   240   0  1.090  5207332    0    119      66     2  0.04169  0.08343

   463 OUT OF    463 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.49
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.423    32.146    34.492  88.619 108.436 111.905
 REC. CELL PARAMETERS   0.041655  0.033721  0.030736  83.877  70.655  67.336
 COORDINATES OF UNIT CELL A-AXIS   -22.755   -11.162   -10.472
 COORDINATES OF UNIT CELL B-AXIS    -6.417    26.338    17.277
 COORDINATES OF UNIT CELL C-AXIS    14.266    20.241   -24.010
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.020     0.011     0.009
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.034
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001348  0.000378
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000613  0.005114  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.30   1244.13
 DETECTOR ORIGIN (PIXELS) AT                     1265.49   1237.41
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.27
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     463 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     31     61     67     55     26     47     65     55
 SIGMAB (degree)       0.032  0.031  0.031  0.031  0.031  0.031  0.031  0.031  0.031
 SIGMAR (degree)       0.028  0.025  0.027  0.028  0.026  0.024  0.026  0.027  0.027


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    463 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  1  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  2  2  1  0  0  0  0     0  0  0  1  4 11 15 11  6  2  1  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  3  3  3  2  1  0  0  0     0  0  0  2  7 17 21 16  8  4  2  1  1
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  1  4 10 12  8  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0 -1  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  1  3  7  9  6  3  1  0  0  0     0  0  0  1  5 11 13  8  4  2  0  0  0     0  0  0  1  4  6  7  5  2  1  0  0  0
     0  0  0  3 15 36 47 32 15  7  3  1  0     0  0  1  6 27 57 67 42 20  9  4  1  0     0  0  0  5 19 35 38 22 11  5  2  0  0
     0  0  1  5 28 60 70 51 25 12  5  2  0     0  0  1 11 47 96100 68 32 15  6  2  0     0  0  1  9 32 60 59 37 18  8  3  1  0
     0  0  0  3 14 36 41 26 14  6  2  1  0     0  0  0  6 25 57 58 38 18  8  3  1  0     0  0  0  4 18 35 35 23 10  4  2  0  0
     0  0  0  1  3  7  7  4  3  1  0  0  0     0  0  0  1  6 13 12  6  3  1  0  0  0     0  0  0  1  5  9  7  3  2  1  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  9  9  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  9 15 14  8  4  2  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  1  5  8  8  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    607 REFLECTION INTENSITIES INTEGRATED


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
   241   0  1.088  5207393    0    119      64     3  0.04159  0.07317
   242   0  1.088  5207972    0    118      61     2  0.04154  0.08157
   243   0  1.088  5207724    0    121      59     3  0.04058  0.07643
   244   0  1.088  5207279    0    115      63     4  0.03906  0.05601
   245   0  1.089  5207475    0    112      63     3  0.04007  0.07056
   246   0  1.089  5208122    0    107      63     3  0.04299  0.10805
   247   0  1.091  5208207    0    111      59     5  0.04129  0.07158
   248   0  1.089  5208012    0    116      58     3  0.04257  0.06800
   249   0  1.090  5207455    0    124      68     2  0.04139  0.07683
   250   0  1.091  5207199    0    132      74     1  0.04403  0.08433
   251   0  1.091  5207251    0    127      71     0  0.03841  0.09120
   252   0  1.088  5207567    0    119      66     3  0.04045  0.07581
   253   0  1.093  5208567    0    125      59     1  0.04870  0.12310
   254   0  1.091  5208962    0    126      49     0  0.04238  0.10593
   255   0  1.091  5206878    0    128      74     5  0.03973  0.06931
   256   0  1.092  5207593    0    131      71     5  0.04119  0.08131
   257   0  1.094  5207840    0    128      61     3  0.04202  0.07405
   258   0  1.094  5207928    0    128      61     1  0.04054  0.09827
   259   0  1.094  5208399    0    128      56     0  0.04098  0.09405
   260   0  1.092  5208370    0    125      57     6  0.04273  0.05680

   467 OUT OF    467 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.53
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.421    32.146    34.490  88.635 108.483 111.904
 REC. CELL PARAMETERS   0.041673  0.033723  0.030749  83.839  70.597  67.326
 COORDINATES OF UNIT CELL A-AXIS   -22.757   -11.173   -10.449
 COORDINATES OF UNIT CELL B-AXIS    -6.416    26.335    17.283
 COORDINATES OF UNIT CELL C-AXIS    14.264    20.237   -24.011
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.014    -0.013     0.005
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.034
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001348  0.000378
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.001071  0.005132  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.44   1244.08
 DETECTOR ORIGIN (PIXELS) AT                     1265.03   1237.34
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.25
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     467 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   49     32     59     67     51     29     52     66     54
 SIGMAB (degree)       0.032  0.032  0.032  0.031  0.031  0.031  0.030  0.031  0.031
 SIGMAR (degree)       0.029  0.025  0.028  0.029  0.027  0.026  0.026  0.027  0.026


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    467 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  3  3  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  3  1  1  0  0  0     0  0  0  1  6 14 18 15  8  3  1  0  0
     0  0  0  1  1  0  1  0  0  0  0  0  0     0  0  0  0  1  4  5  4  2  1  0  0  0     0  0  0  1  9 23 30 23 13  6  2  1  0
    -1  0  0  1  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  3  2  1  0  0  0     0  0  0  1  6 13 19 14  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  4  4  3  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0 -1  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  4  9 10  8  4  1  0  0  0     0  0  0  1  6 15 13  9  4  1  0  0  0     0  0  0  1  4  9  7  5  2  1  0  0  0
     0  0  0  3 19 44 54 37 19  9  3  1  0     0  0  1  6 29 63 70 43 22 10  4  1  0     0  0  0  4 19 36 35 21 11  4  2  1  0
     0  0  1  6 32 70 82 59 30 14  6  2  0     0  0  1 11 45 94100 66 32 15  6  2  0     0  0  1  7 27 52 50 32 15  7  3  1  0
     0  0  0  4 18 39 47 36 17  8  3  1  0     0  0  1  5 28 56 57 40 19  8  3  1  0     0  0  0  4 16 31 30 18  9  4  1  0  0
     0  0  0  1  4 11 11  7  4  1  0  0  0     0  0  0  2  6 14 14  9  5  2  0  0  0     0  0  0  1  4  8  7  4  2  1  0  0  0
     0  0  0  0  1  1  1  1  1  0  0 -1  0     0  0  0  0  1  1  2  1  1  0  0 -1  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0
     0  0  0  1  5  9  8  5  2  1  0  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  1  3  8 12 11  7  3  1  0  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  1  5  8  7  4  2  1  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  1  0  0  0  0     0  0  0  0  1  0  1  1  0  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    624 REFLECTION INTENSITIES INTEGRATED


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
   261   0  1.090  5207365    0    131      71     4  0.03913  0.06337
   262   0  1.090  5206934    0    127      76     5  0.03784  0.05301
   263   0  1.090  5206905    0    130      78     8  0.03763  0.06121
   264   0  1.090  5208080    0    130      65     7  0.04298  0.07134
   265   0  1.093  5207625    0    119      68     7  0.04226  0.08759
   266   0  1.093  5208679    0    120      60     6  0.04285  0.08947
   267   0  1.092  5208178    0    124      62     7  0.04034  0.07787
   268   0  1.087  5207503    0    119      73     8  0.04052  0.07709
   269   0  1.089  5207677    0    127      70     2  0.04142  0.09188
   270   0  1.091  5207286    0    135      68     7  0.03895  0.06761
   271   0  1.090  5207794    0    127      65     7  0.04146  0.08930
   272   0  1.088  5207935    0    130      62     4  0.04420  0.08774
   273   0  1.091  5208041    0    118      61     2  0.04461  0.10076
   274   0  1.092  5208348    0    106      57     3  0.04211  0.10163
   275   0  1.087  5208475    0    108      54     3  0.04145  0.06805
   276   0  1.089  5208519    0    125      55     2  0.04052  0.08282
   277   0  1.087  5207645    0    134      63     4  0.04027  0.05356
   278   0  1.085  5207717    0    131      65     4  0.03755  0.06304
   279   0  1.087  5208237    0    135      61     6  0.03836  0.05462
   280   0  1.087  5208457    0    117      62     7  0.04135  0.07179

   467 OUT OF    469 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.50
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.13
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.421    32.147    34.490  88.640 108.462 111.906
 REC. CELL PARAMETERS   0.041667  0.033723  0.030744  83.841  70.618  67.325
 COORDINATES OF UNIT CELL A-AXIS   -22.755   -11.164   -10.464
 COORDINATES OF UNIT CELL B-AXIS    -6.417    26.327    17.296
 COORDINATES OF UNIT CELL C-AXIS    14.266    20.246   -24.003
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.005     0.005     0.001
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.034
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001348  0.000378
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000892  0.005109  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.39   1244.04
 DETECTOR ORIGIN (PIXELS) AT                     1265.22   1237.32
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.26
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     467 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     32     56     65     49     28     53     67     56
 SIGMAB (degree)       0.032  0.031  0.032  0.032  0.031  0.031  0.030  0.031  0.031
 SIGMAR (degree)       0.027  0.024  0.026  0.027  0.026  0.022  0.026  0.028  0.027


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    467 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  1  0  0  0     0  0  0  1  4 11 16 13  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  4  3  2  1  0  0  0     0  0  0  1  7 19 25 20 11  5  2  1  1
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  1  5 12 16 11  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  3  4  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7 10  6  3  2  0  0  0     0  0  0  1  6 12 15  8  4  2  0  0  0     0  0  0  1  3  7  8  4  2  1  0  0  0
     0  0  0  3 18 38 49 36 19  8  3  1  0     0  0  1  6 31 59 66 44 22  9  3  1  0     0  0  0  4 19 35 35 22 11  4  1  0  0
     0  0  1  6 28 63 75 54 29 12  5  2  0     0  0  1 12 48 95100 66 33 15  6  2  0     0  0  1  8 29 53 54 34 16  7  3  1  0
     0  0  0  4 16 41 47 31 14  7  3  1  0     0  0  1  7 25 59 62 39 19  9  3  1  1     0  0  0  5 15 34 33 20 10  4  1  0  0
     0  0  0  1  4  9 10  6  3  1  0  0  0     0  0  0  2  7 13 14  8  3  1  0  0  0     0  0  0  1  5  8  8  5  2  1  0  0  0
     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0 -1  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  8  7  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 11 11  7  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  8  7  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    617 REFLECTION INTENSITIES INTEGRATED


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
   281   0  1.089  5209154    0    102      48     2  0.03682  0.08347
   282   0  1.088  5209112    0    118      54     6  0.04255  0.07564
   283   0  1.088  5208900    0    111      52     6  0.03903  0.05854
   284   0  1.087  5208724    0    107      52     6  0.03698  0.04834
   285   0  1.089  5209089    0    108      48     5  0.04094  0.06964
   286   0  1.087  5208669    0    113      53     3  0.03886  0.07160
   287   0  1.084  5208346    0    120      60     2  0.03906  0.06907
   288   0  1.085  5208950    0    123      54     1  0.03920  0.10320
   289   0  1.085  5208821    0    120      54     2  0.03932  0.06252
   290   0  1.085  5208378    0    118      63     3  0.03844  0.07328
   291   0  1.084  5209141    0    117      52     5  0.03935  0.06179
   292   0  1.084  5208692    0    131      58     2  0.04048  0.07567
   293   0  1.082  5207576    0    126      70     2  0.03603  0.06024
   294   0  1.082  5208296    0    139      63     1  0.04309  0.07969
   295   0  1.083  5208172    0    143      64     5  0.03830  0.06687
   296   0  1.084  5207812    0    133      69     7  0.03718  0.05249
   297   0  1.084  5208166    0    135      62     5  0.03720  0.05711
   298   0  1.082  5208305    0    135      57     2  0.03797  0.05778
   299   0  1.082  5208823    0    121      54     6  0.03834  0.05412
   300   0  1.084  5208202    0    120      62     4  0.03927  0.06956

   464 OUT OF    466 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.48
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.06
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.428    32.152    34.506  88.620 108.436 111.898
 REC. CELL PARAMETERS   0.041645  0.033713  0.030723  83.879  70.656  67.343
 COORDINATES OF UNIT CELL A-AXIS   -22.761   -11.161   -10.472
 COORDINATES OF UNIT CELL B-AXIS    -6.418    26.341    17.283
 COORDINATES OF UNIT CELL C-AXIS    14.271    20.251   -24.019
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.014     0.013     0.005
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.034
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001348  0.000378
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000593  0.005195  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.30   1244.11
 DETECTOR ORIGIN (PIXELS) AT                     1265.52   1237.28
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.33
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     464 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   46     32     58     67     52     26     50     66     60
 SIGMAB (degree)       0.031  0.031  0.032  0.032  0.031  0.032  0.030  0.031  0.031
 SIGMAR (degree)       0.027  0.024  0.026  0.027  0.026  0.023  0.026  0.027  0.026


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    464 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  2  2  1  0  0  0  0     0  0  0  0  3  8 12  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0     0  0  0  1  6 14 19 16  8  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  1  4  9 12  8  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  0  1  0  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  1  1  1  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0 -1 -1  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0 -1  0     0  0  0  0  1  1  1  0  0  0  0 -1  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  3  8 10  6  3  1  0  0  0     0  0  0  1  6 14 16  8  4  2  1  0  0     0  0  0  1  5  9  9  4  2  1  0  0  0
     0  0  0  3 13 32 41 31 16  7  3  1  0     0  0  0  6 28 58 65 44 21  9  4  1  0     0  0  1  4 20 35 35 22 10  5  2  0  0
     0  0  0  5 23 52 67 51 26 11  4  2  0     0  0  1 11 46 91100 69 35 15  5  2  0     0  0  1  8 32 55 55 35 17  7  3  1  0
     0  0  0  3 15 34 42 28 14  6  2  1  0     0  0  1  7 30 59 63 40 19  9  3  1  0     0  0  1  5 22 36 35 21 10  5  1  0  0
     0  0  0  1  3  8  8  4  3  1  0  0  0     0  0  0  2  6 13 13  7  4  1  0  0  0     0  0  0  1  4  9  7  4  2  1  0  0  0
     0  0  0  0  0  0  1  1  0  0  0 -1  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0 -1  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  9  8  5  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  8 13 13  8  4  2  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6  9  8  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
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
   301   0  1.083  5208441    0    115      60     5  0.03902  0.06247
   302   0  1.080  5208668    0    126      62     4  0.04069  0.07913
   303   0  1.080  5207969    0    130      65     6  0.04238  0.06516
   304   0  1.080  5207225    0    122      76     8  0.03704  0.05118
   305   0  1.079  5207837    0    120      70     3  0.03958  0.09060
   306   0  1.083  5208305    0    113      64     3  0.04291  0.09472
   307   0  1.083  5208860    0    110      55     2  0.03656  0.06379
   308   0  1.082  5208896    0    124      50     2  0.04117  0.06948
   309   0  1.084  5208915    0    120      56     4  0.03869  0.06465
   310   0  1.083  5208527    0    117      59     6  0.04103  0.06278
   311   0  1.086  5208735    0    116      58     6  0.04188  0.08569
   312   0  1.083  5208384    0    110      63     4  0.04351  0.11031
   313   0  1.082  5208386    0    128      63     7  0.03876  0.05533
   314   0  1.082  5207538    0    143      71     7  0.03764  0.05383
   315   0  1.084  5207912    0    131      68     5  0.03940  0.06399
   316   0  1.084  5208234    0    133      69     4  0.03981  0.07691
   317   0  1.085  5209095    0    129      50     5  0.03776  0.07884
   318   0  1.084  5208659    0    123      58     5  0.03971  0.07683
   319   0  1.083  5208002    0    130      66     3  0.03888  0.08657
   320   0  1.083  5207704    0    122      66     7  0.03772  0.05272

   474 OUT OF    475 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.51
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.06
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.424    32.148    34.491  88.639 108.425 111.897
 REC. CELL PARAMETERS   0.041651  0.033718  0.030736  83.862  70.659  67.339
 COORDINATES OF UNIT CELL A-AXIS   -22.757   -11.156   -10.475
 COORDINATES OF UNIT CELL B-AXIS    -6.419    26.332    17.289
 COORDINATES OF UNIT CELL C-AXIS    14.266    20.240   -24.009
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.009     0.013     0.003
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.034
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001348  0.000378
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000582  0.005072  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.31   1244.02
 DETECTOR ORIGIN (PIXELS) AT                     1265.54   1237.36
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.27
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     474 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   49     31     59     70     53     27     53     65     55
 SIGMAB (degree)       0.032  0.031  0.032  0.032  0.031  0.031  0.030  0.030  0.030
 SIGMAR (degree)       0.026  0.023  0.026  0.027  0.025  0.023  0.026  0.027  0.025


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    474 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  1  3  8 12  8  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  1  6 15 20 14  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  1  3  9 12  8  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  7  9  6  3  1  0  0  0     0  0  0  1  5 11 13  8  4  2  0  0  0     0  0  0  1  4  6  7  4  2  1  0  0  0
     0  0  0  3 12 32 42 27 14  7  2  1  0     0  0  0  5 25 55 64 39 20  9  3  1  0     0  0  0  4 16 31 34 21 10  5  2  0  0
     0  0  0  4 23 55 68 47 24 11  4  2  0     0  0  1 10 43 91100 66 33 15  6  2  1     0  0  1  7 28 53 52 34 17  7  3  1  0
     0  0  0  3 15 33 40 27 14  6  2  1  0     0  0  1  6 28 56 59 38 18  8  3  1  0     0  0  1  5 19 32 32 20  9  4  1  0  0
     0  0  0  0  3  6 10  6  4  1  0  0  0     0  0  0  1  7 12 15  9  5  1  0  0  0     0  0  0  1  5  7  8  5  2  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  1  2  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0 -1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  7  4  2  1  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 12 11  7  4  2  1  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  7  7  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    629 REFLECTION INTENSITIES INTEGRATED


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
   321   0  1.087  5208334    0    123      61     2  0.03963  0.07163
   322   0  1.086  5208099    0    123      64     3  0.04093  0.08161
   323   0  1.083  5207691    0    118      61     7  0.03773  0.04930
   324   0  1.085  5208017    0    126      62     5  0.04109  0.06977
   325   0  1.087  5208264    0    133      64     5  0.04391  0.07724
   326   0  1.089  5207452    0    148      73     7  0.03712  0.05940
   327   0  1.084  5207277    0    155      76     3  0.03875  0.07099
   328   0  1.087  5207428    0    144      73     4  0.03969  0.06681
   329   0  1.087  5208560    0    136      58     4  0.03892  0.07766
   330   0  1.087  5209583    0    115      50     5  0.03594  0.06516
   331   0  1.087  5209147    0    105      50     4  0.03765  0.06054
   332   0  1.087  5209907    0    100      42     5  0.04022  0.07212
   333   0  1.084  5209721    0    105      48     5  0.04211  0.07432
   334   0  1.086  5208936    0    112      52     6  0.03704  0.05081
   335   0  1.088  5208651    0    106      61     7  0.03848  0.05683
   336   0  1.086  5208454    0    113      60     5  0.04173  0.07828
   337   0  1.086  5208002    0    119      67     7  0.04116  0.05934
   338   0  1.088  5208024    0    125      67     3  0.03784  0.07213
   339   0  1.084  5207675    0    129      68     3  0.03752  0.06594
   340   0  1.085  5207886    0    126      68     7  0.04135  0.06119

   462 OUT OF    463 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.56
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.07
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.424    32.147    34.485  88.628 108.427 111.893
 REC. CELL PARAMETERS   0.041648  0.033717  0.030741  83.875  70.662  67.346
 COORDINATES OF UNIT CELL A-AXIS   -22.757   -11.154   -10.476
 COORDINATES OF UNIT CELL B-AXIS    -6.419    26.334    17.284
 COORDINATES OF UNIT CELL C-AXIS    14.265    20.238   -24.002
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.011     0.016     0.004
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.034
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001348  0.000378
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000476  0.005072  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.29   1244.10
 DETECTOR ORIGIN (PIXELS) AT                     1265.67   1237.43
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.26
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     462 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     30     51     64     53     34     55     64     55
 SIGMAB (degree)       0.032  0.032  0.032  0.031  0.031  0.031  0.031  0.031  0.031
 SIGMAR (degree)       0.028  0.024  0.026  0.027  0.026  0.024  0.025  0.027  0.025


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    462 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  1  4  9 12  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0 -1     0  0  0  0  1  2  3  3  1  0  0  0 -1     0  0  0  1  6 15 21 15  8  4  2  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0     0  0  0  0  4  9 13 10  5  2  1  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  8  8  6  3  1  0  0  0     0  0  0  1  6 14 12  8  4  2  0  0  0     0  0  0  1  4  7  5  3  2  1  0  0  0
     0  0  0  2 15 34 42 30 16  6  3  1  0     0  0  0  5 27 59 65 42 22  9  4  1  0     0  0  0  3 17 33 33 20 10  4  1  0  0
     0  0  0  3 22 54 69 49 26 12  5  2  0     0  0  1  8 43 91100 66 34 15  6  2  1     0  0  1  6 26 49 47 30 16  7  2  1  0
     0  0  0  2 15 33 42 30 16  7  3  1  0     0  0  0  5 28 54 60 41 20  9  3  1  0     0  0  0  4 17 29 29 20  9  4  1  0  0
     0  0  0  1  3  8  8  5  3  1  0  0  0     0  0  0  2  7 14 13  8  4  1  0  0  0     0  0  0  1  5  9  7  4  2  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  7  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5 10  9  6  2  1  0  0  0     0  0  0  1  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4  6  5  4  1  0  0  0  0     0  0  0  1  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  1  0  0  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0 -1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    607 REFLECTION INTENSITIES INTEGRATED


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
   341   0  1.087  5209020    0    128      52     6  0.04274  0.06842
   342   0  1.086  5208348    0    131      62     3  0.04085  0.08608
   343   0  1.084  5208169    0    121      65     1  0.03909  0.10178
   344   0  1.083  5208840    0    131      51     0  0.03958  0.12378
   345   0  1.082  5207959    0    122      63     5  0.04016  0.06992
   346   0  1.085  5207355    0    131      71     8  0.04081  0.05560
   347   0  1.086  5207816    0    129      62     7  0.03950  0.05585
   348   0  1.083  5207954    0    122      65     2  0.04531  0.12748
   349   0  1.084  5208055    0    116      61     0  0.04237  0.12619
   350   0  1.083  5208172    0    104      60     3  0.04280  0.08980
   351   0  1.081  5208309    0    113      58     2  0.03875  0.07492
   352   0  1.080  5207933    0    136      67     1  0.04666  0.09878
   353   0  1.080  5206606    0    133      76     3  0.04143  0.08081
   354   0  1.082  5206765    0    135      75     8  0.03966  0.06481
   355   0  1.082  5208098    0    131      64     6  0.03993  0.07099
   356   0  1.081  5208376    0    118      57     5  0.04137  0.05723
   357   0  1.085  5207977    0    116      58     3  0.04320  0.05628
   358   0  1.084  5208479    0    124      59     2  0.04153  0.08077
   359   0  1.085  5207389    0    120      67     4  0.03988  0.06284
   360   0  1.081  5208337    0    132      61     4  0.04407  0.07308

   475 OUT OF    475 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.49
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.07
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.419    32.143    34.484  88.642 108.446 111.901
 REC. CELL PARAMETERS   0.041665  0.033725  0.030747  83.848  70.634  67.331
 COORDINATES OF UNIT CELL A-AXIS   -22.754   -11.160   -10.465
 COORDINATES OF UNIT CELL B-AXIS    -6.416    26.330    17.285
 COORDINATES OF UNIT CELL C-AXIS    14.265    20.233   -24.006
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.013     0.006     0.004
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.034
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001348  0.000378
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000773  0.005091  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.34   1244.08
 DETECTOR ORIGIN (PIXELS) AT                     1265.32   1237.39
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.23
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     475 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   53     30     57     67     53     29     57     66     58
 SIGMAB (degree)       0.032  0.032  0.032  0.032  0.031  0.031  0.030  0.031  0.031
 SIGMAR (degree)       0.027  0.023  0.026  0.028  0.025  0.023  0.025  0.027  0.026


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    475 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  1  3 10 13 11  5  3  1  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  1  6 16 19 16  8  4  1  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  2  2  2  1  0  0  0  0     0  0  0  1  4 11 13 10  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  3  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  0  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  8  9  6  3  2  0  0  0     0  0  0  2  8 15 14  9  4  2  0  0  0     0  0  0  1  5  9  7  5  2  1  0  0  0
     0  0  0  3 15 37 46 35 17  8  3  1  0     0  0  1  6 28 61 68 46 22  9  4  1  0     0  0  1  4 17 36 36 22 11  4  2  0  0
     0  0  1  5 24 60 69 52 26 12  5  2  0     0  0  1 10 45 96100 68 32 15  6  2  0     0  0  1  7 28 55 54 35 17  7  3  1  0
     0  0  0  4 18 41 45 30 14  7  3  1  0     0  0  1  7 32 63 63 39 18  8  3  1  0     0  0  0  4 19 36 35 21  9  4  2  1  0
     0  0  0  1  5  9  9  7  3  1  0  0  0     0  0  0  2  9 15 13 10  4  2  0  0  0     0  0  0  1  5  8  7  5  2  1  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0 -1 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0 -1  0  0  0  0  0  0  0  0  0     0  0  0 -1  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  2  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  8  7  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 12 11  7  3  1  1  0  0     0  0  1  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  8  7  4  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    644 REFLECTION INTENSITIES INTEGRATED


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
   361   0  1.082  5207479    0    137      61     6  0.03916  0.05536
   362   0  1.082  5207334    0    133      65     2  0.04386  0.08850
   363   0  1.086  5207979    0    134      60     4  0.04190  0.07436
   364   0  1.084  5208635    0    129      55     2  0.04428  0.11065
   365   0  1.086  5208442    0    123      56     0  0.04097  0.11009
   366   0  1.084  5207840    0    123      63     2  0.04075  0.06809
   367   0  1.085  5208038    0    120      60     2  0.04264  0.07722
   368   0  1.086  5208226    0    122      61     4  0.04024  0.06937
   369   0  1.087  5208820    0    119      51     2  0.03996  0.07377
   370   0  1.086  5207697    0    101      65     4  0.04000  0.06924
   371   0  1.086  5208269    0    109      60     4  0.04058  0.08352
   372   0  1.085  5208151    0    109      58     5  0.04136  0.08460
   373   0  1.084  5208226    0    107      63     7  0.03976  0.07513
   374   0  1.087  5208056    0    121      64     7  0.04183  0.07664
   375   0  1.086  5207232    0    128      74     6  0.04212  0.07301
   376   0  1.085  5207603    0    137      64     4  0.04124  0.07052
   377   0  1.085  5207323    0    137      71     6  0.04072  0.07326
   378   0  1.085  5207804    0    141      62     4  0.03968  0.08116
   379   0  1.085  5207963    0    143      65     6  0.03961  0.06651
   380   0  1.085  5207554    0    143      63     6  0.04071  0.06254

   463 OUT OF    463 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.50
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.08
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.421    32.146    34.487  88.624 108.446 111.895
 REC. CELL PARAMETERS   0.041659  0.033719  0.030743  83.870  70.643  67.343
 COORDINATES OF UNIT CELL A-AXIS   -22.756   -11.160   -10.465
 COORDINATES OF UNIT CELL B-AXIS    -6.417    26.339    17.276
 COORDINATES OF UNIT CELL C-AXIS    14.267    20.233   -24.009
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.024     0.006     0.006
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.034
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001348  0.000378
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000638  0.005241  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.34   1244.10
 DETECTOR ORIGIN (PIXELS) AT                     1265.50   1237.22
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.25
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     463 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   46     28     60     70     55     30     52     65     48
 SIGMAB (degree)       0.032  0.031  0.032  0.032  0.031  0.031  0.030  0.031  0.031
 SIGMAR (degree)       0.026  0.023  0.026  0.028  0.025  0.022  0.026  0.028  0.027


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    463 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  1  1  1  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  1  4  9 13 10  5  2  1  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  1  6 15 20 15  9  4  2  1  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  1  4 10 13  9  5  2  1  0  0
     0  0  0  0  0  1  0  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0 -1 -1  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0 -1  0  0 -1  0  0  0     0  0  0  0  0  0 -1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  6  9  6  3  1  0  0  0     0  0  0  1  6 11 14  8  4  1  0  0  0     0  0  0  1  4  6  7  4  2  1  0  0  0
     0  0  0  3 13 34 44 32 16  7  3  1  0     0  0  1  6 26 59 63 41 20  8  3  1  0     0  0  1  5 17 34 30 19 10  4  1  0  0
     0  0  1  4 23 56 70 48 27 12  5  2  0     0  0  1  9 43 93100 65 35 15  6  2  0     1  0  1  7 27 51 50 32 16  7  3  1  0
     0  0  0  3 15 35 45 30 17  6  2  1  0     0  0  1  6 28 60 68 42 21  9  3  1  0     0  0  0  5 18 32 36 22 10  4  1  0  0
     0  0  0  1  4  9 10  5  2  1  0  0  0     0  0  0  2  8 16 16  9  3  2  0  0  0     0  0  0  1  5  9  8  5  2  1  0  0  0
     0  0  0  0  0  1  0  0  0  0  0 -1  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0
     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0 -1  0 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0 -1  0  0  0  0 -1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  1  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  1  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  8  7  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 12 12  7  4  2  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  8  8  5  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    645 REFLECTION INTENSITIES INTEGRATED


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
   381   0  1.086  5207869    0    136      62     7  0.03861  0.04524
   382   0  1.086  5207976    0    137      64     6  0.04075  0.05642
   383   0  1.087  5207638    0    136      71     6  0.03974  0.06022
   384   0  1.088  5207795    0    134      67     3  0.04035  0.07336
   385   0  1.087  5208025    0    132      61     3  0.04264  0.07257
   386   0  1.088  5208576    0    127      57     5  0.04084  0.06958
   387   0  1.088  5208105    0    128      70     7  0.04152  0.07531
   388   0  1.087  5207529    0    127      74     5  0.04167  0.08283
   389   0  1.087  5208305    0    124      60     5  0.04126  0.08600
   390   0  1.087  5208740    0    130      58     4  0.03990  0.09382
   391   0  1.086  5207524    0    119      77     8  0.04030  0.07422
   392   0  1.085  5208167    0    129      68     2  0.04201  0.11202
   393   0  1.086  5208860    0    118      52     3  0.03801  0.08195
   394   0  1.089  5207738    0    115      66     5  0.04062  0.07075
   395   0  1.090  5209099    0    113      54     3  0.04433  0.08710
   396   0  1.089  5208164    0    113      68     4  0.04339  0.08541
   397   0  1.086  5208412    0    118      60     6  0.04096  0.06284
   398   0  1.086  5208121    0    125      61     6  0.03901  0.04506
   399   0  1.087  5208251    0    127      61     5  0.04162  0.04887
   400   0  1.088  5208022    0    130      62     6  0.04137  0.06425

   482 OUT OF    483 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.54
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.14
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.424    32.149    34.483  88.654 108.442 111.893
 REC. CELL PARAMETERS   0.041655  0.033718  0.030747  83.839  70.635  67.335
 COORDINATES OF UNIT CELL A-AXIS   -22.759   -11.159   -10.468
 COORDINATES OF UNIT CELL B-AXIS    -6.420    26.331    17.292
 COORDINATES OF UNIT CELL C-AXIS    14.266    20.231   -24.005
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.010     0.005     0.003
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.034
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001348  0.000378
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000776  0.005049  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.31   1244.05
 DETECTOR ORIGIN (PIXELS) AT                     1265.29   1237.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.28
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     482 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   51     32     57     69     56     30     54     65     58
 SIGMAB (degree)       0.032  0.032  0.032  0.031  0.031  0.031  0.031  0.031  0.031
 SIGMAR (degree)       0.027  0.024  0.027  0.028  0.026  0.024  0.026  0.027  0.026


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    482 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0 -1  0  0  0 -1 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  1  0  0  0  0     0  0  0  0  1  2  3  3  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  1  5 12 18 14  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  2  3  4  4  2  1  1  0  0     0  0  0  1  9 22 29 23 12  6  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  3  3  2  1  0  0  0  0     0  0  0  1  6 15 18 13  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  3  4  3  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0 -1 -1  0    -1  0  0  0  1  1  1  1  0  0 -1 -1  0     0  0  0  0  0  1  0  1  0  0  0  0  0
     0  0  0  1  4  9 11  7  4  1  0  0  0     0  0  0  1  6 13 14  9  4  1  0  0  0     0  0  0  1  4  7  6  4  2  1  0  0  0
     0  0  0  3 18 42 56 38 19  8  3  1  0     0  0  1  5 28 57 68 43 21 10  3  1  0     0  0  1  4 15 27 29 18  9  4  1  0  0
     0  0  1  5 33 73 86 63 33 14  6  2  0     0  0  1 10 49 95100 68 34 15  6  2  0     0  0  1  7 24 43 43 28 14  6  2  0  0
     0  0  0  5 22 48 54 38 18  8  3  1  0     0  0  0  8 33 61 62 42 19  8  3  1  0     0  0  1  5 18 29 28 18  8  4  1  0  0
     0  0  0  1  5 10 11  8  3  1  0  0  0     0  0  0  2  7 14 14  9  3  1  0  0  0     0  0  0  1  4  7  7  4  1  1  0  0  0
     0  0  0  0  1  1  1  1  0  0  0 -1  0     0  0  0  0  1  2  2  1  0  0  0 -1  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  1  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  1  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  1  4  7  6  4  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  2  6  9  8  6  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  2  4  6  6  4  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    640 REFLECTION INTENSITIES INTEGRATED


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
   401   0  1.087  5208086    0    120      65     2  0.04222  0.07264
   402   0  1.084  5208251    0    114      60     3  0.03984  0.06538
   403   0  1.085  5209257    0    115      46     4  0.03770  0.06150
   404   0  1.086  5207616    0    134      65     1  0.03891  0.06617
   405   0  1.086  5207040    0    144      77     1  0.03757  0.06413
   406   0  1.084  5207000    0    146      74     3  0.04114  0.06376
   407   0  1.087  5208274    0    132      61     3  0.04151  0.07088
   408   0  1.090  5207961    0    129      67     3  0.04171  0.06743
   409   0  1.090  5208146    0    118      63     7  0.04087  0.05221
   410   0  1.089  5209034    0    123      57     4  0.04481  0.07591
   411   0  1.087  5208385    0    128      59     3  0.04220  0.07453
   412   0  1.085  5207644    0    121      69     6  0.04020  0.06070
   413   0  1.087  5207897    0    124      68     3  0.04447  0.09606
   414   0  1.085  5208893    0    119      53     1  0.04201  0.09524
   415   0  1.087  5208167    0    120      60     0  0.04130  0.08751
   416   0  1.084  5207713    0    132      67     2  0.04302  0.09549
   417   0  1.086  5208277    0    130      57     2  0.04149  0.07263
   418   0  1.085  5207788    0    122      68     0  0.04152  0.09577
   419   0  1.085  5208304    0    124      61     1  0.04465  0.12578
   420   0  1.087  5207721    0    138      74     8  0.04105  0.06870

   466 OUT OF    468 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.49
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.03
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.422    32.144    34.490  88.643 108.446 111.885
 REC. CELL PARAMETERS   0.041656  0.033720  0.030741  83.853  70.636  67.347
 COORDINATES OF UNIT CELL A-AXIS   -22.759   -11.157   -10.463
 COORDINATES OF UNIT CELL B-AXIS    -6.418    26.330    17.285
 COORDINATES OF UNIT CELL C-AXIS    14.266    20.237   -24.010
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.010     0.004     0.000
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.034
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001348  0.000378
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000714  0.005091  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.36   1244.08
 DETECTOR ORIGIN (PIXELS) AT                     1265.42   1237.39
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.25
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     466 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   48     31     63     68     50     29     52     65     55
 SIGMAB (degree)       0.032  0.032  0.032  0.032  0.031  0.032  0.031  0.031  0.030
 SIGMAR (degree)       0.026  0.023  0.026  0.028  0.025  0.022  0.025  0.027  0.026


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    466 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  3  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  1  3 10 14 10  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  1  5 14 19 15  8  4  1  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  1  4 10 12  9  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  3  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0 -1  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0 -1  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  4  9 10  7  3  1  0  0  0     0  0  0  1  7 15 16 10  4  2  0  0  0     0  0  0  1  4  9  8  5  2  1  0  0  0
     0  0  0  3 14 37 48 33 17  8  3  1  0     0  0  1  7 28 62 68 45 22 10  4  1  0     0  0  1  5 18 35 34 21  9  4  2  0  0
     0  0  0  4 21 54 67 49 26 12  5  2  0     0  0  1 10 41 90100 67 33 14  6  2  0     0  0  1  7 26 51 51 32 15  6  3  0  0
     0  0  0  3 16 34 41 29 14  6  3  1  0     0  0  1  6 29 57 62 38 19  8  3  1  0     0  0  0  4 18 32 31 18  9  3  1  0  0
     0  0  0  1  4  9  9  7  4  1  0  0  0     0  0  0  2  8 14 15 10  5  2  0  0  0     0  0  0  2  6  8  8  5  2  1  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     1  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4  8  7  4  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 12 11  6  3  1  1  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  1  5  7  7  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  1  2  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    642 REFLECTION INTENSITIES INTEGRATED


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
   421   0  1.087  5207888    0    130      69     5  0.04122  0.06156
   422   0  1.085  5208795    0    121      57     2  0.04134  0.07657
   423   0  1.087  5209372    0    118      49     3  0.04021  0.06049
   424   0  1.086  5208670    0    132      59     3  0.03785  0.08675
   425   0  1.082  5207800    0    129      67     2  0.03930  0.06743
   426   0  1.085  5208330    0    139      56     4  0.03803  0.04849
   427   0  1.085  5207995    0    125      66     6  0.03988  0.05465
   428   0  1.085  5208997    0    121      51     2  0.03726  0.08652
   429   0  1.086  5209281    0    118      51     6  0.03780  0.06172
   430   0  1.084  5208886    0    114      52     6  0.03607  0.04768
   431   0  1.084  5208307    0    117      56     4  0.03897  0.05789
   432   0  1.084  5208779    0    115      56     5  0.04181  0.08413
   433   0  1.081  5208821    0    117      57     2  0.04440  0.13443
   434   0  1.084  5208451    0    117      62     0  0.04427  0.12503
   435   0  1.083  5208516    0    122      56     0  0.03986  0.09416
   436   0  1.084  5208094    0    137      58     6  0.03952  0.05985
   437   0  1.079  5207551    0    126      73     2  0.03816  0.07527
   438   0  1.081  5207878    0    126      71     6  0.04321  0.07630
   439   0  1.081  5208328    0    116      62     5  0.04037  0.08921
   440   0  1.080  5208552    0    126      59     6  0.04239  0.07851

   454 OUT OF    458 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.50
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.13
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.420    32.146    34.505  88.639 108.449 111.890
 REC. CELL PARAMETERS   0.041660  0.033719  0.030729  83.854  70.634  67.343
 COORDINATES OF UNIT CELL A-AXIS   -22.758   -11.163   -10.457
 COORDINATES OF UNIT CELL B-AXIS    -6.419    26.333    17.284
 COORDINATES OF UNIT CELL C-AXIS    14.264    20.246   -24.025
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.012    -0.004     0.003
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.034
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001348  0.000378
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000951  0.005116  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.34   1244.06
 DETECTOR ORIGIN (PIXELS) AT                     1265.09   1237.34
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.28
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     454 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   48     31     54     60     57     29     50     60     53
 SIGMAB (degree)       0.031  0.032  0.032  0.031  0.031  0.031  0.031  0.031  0.031
 SIGMAR (degree)       0.027  0.023  0.026  0.027  0.025  0.022  0.025  0.027  0.026


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    454 STRONG REFLECTIONS *****


     0  0  1  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  2  2  1  0  0  0  0     0  0  0  0  3  8 11  9  4  2  1  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  2  2  3  1  1  0  0  0     0  0  0  1  6 14 16 13  7  3  1  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     1  0  0  1  3  8 10  7  4  2  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0 -1  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  1  1  1  0  0 -1 -1  0     0  0  0  0  0  1  1  1  0  0  0 -1  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  3  6  9  5  3  1  0  0  0     0  0  0  1  7 13 15  7  4  2  0  0  0     0  0  0  1  6 10 10  5  3  1  0  0  0
     0  0  0  2 14 32 41 30 14  6  3  1  0     0  0  1  5 28 57 67 43 20  9  4  1  0     0  0  1  5 22 38 39 22 11  5  2  1  0
     0  0  1  4 25 56 63 45 23 11  4  1  0     0  0  1  9 47 94100 66 32 15  6  2  0     0  0  1  8 34 62 60 37 18  8  3  1  0
     0  0  0  3 14 34 38 26 13  6  2  0  0     0  0  1  6 27 58 59 38 19  9  3  1  0     0  0  1  5 20 40 37 23 11  5  2  0  0
     0  0  0  1  3  6  7  6  2  1  0  0  0     0  0  0  2  7 11 12  9  3  1  0  0  0     0  0  0  2  6  9  8  6  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0 -1 -1  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1 -1  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0 -1  0 -1  0  0  0  0  0     0  0  0  0  0  0  0 -1  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  6  9  9  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  9 14 14  8  4  2  1  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5 10  8  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    635 REFLECTION INTENSITIES INTEGRATED


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
   441   0  1.078  5208273    0    128      67     7  0.03902  0.07878
   442   0  1.076  5207491    0    127      73     8  0.03843  0.05528
   443   0  1.078  5208085    0    131      64     6  0.04315  0.06151
   444   0  1.077  5208661    0    132      59     0  0.04341  0.09510
   445   0  1.081  5207489    0    133      67     1  0.03868  0.05576
   446   0  1.077  5207821    0    136      70     0  0.03850  0.07046
   447   0  1.078  5208348    0    134      64     0  0.04142  0.09766
   448   0  1.077  5208020    0    125      63     2  0.04087  0.07932
   449   0  1.075  5207912    0    131      67     4  0.04235  0.08425
   450   0  1.074  5207404    0    124      70     7  0.04119  0.05596
   451   0  1.075  5207641    0    121      67     7  0.04046  0.07257
   452   0  1.073  5207960    0    123      61     4  0.04478  0.09536
   453   0  1.074  5207542    0    117      68     6  0.04217  0.07269
   454   0  1.074  5208433    0    115      56     6  0.04021  0.07960
   455   0  1.075  5209291    0    104      52     3  0.04391  0.09916
   456   0  1.075  5209166    0     92      51     0  0.04048  0.11155
   457   0  1.074  5208905    0    109      55     0  0.04218  0.11525
   458   0  1.074  5208382    0    104      62     6  0.04198  0.06983
   459   0  1.072  5209071    0    112      51     2  0.04082  0.10631
   460   0  1.073  5208442    0    110      55     2  0.04295  0.09801

   465 OUT OF    465 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.54
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.422    32.144    34.498  88.652 108.460 111.894
 REC. CELL PARAMETERS   0.041664  0.033723  0.030738  83.833  70.616  67.334
 COORDINATES OF UNIT CELL A-AXIS   -22.759   -11.166   -10.455
 COORDINATES OF UNIT CELL B-AXIS    -6.418    26.327    17.289
 COORDINATES OF UNIT CELL C-AXIS    14.264    20.241   -24.019
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.006    -0.009     0.002
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.034
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001348  0.000378
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.001050  0.005074  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.36   1244.06
 DETECTOR ORIGIN (PIXELS) AT                     1264.98   1237.40
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.25
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     465 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     30     51     63     50     29     57     69     62
 SIGMAB (degree)       0.032  0.032  0.032  0.032  0.031  0.031  0.030  0.031  0.031
 SIGMAR (degree)       0.028  0.023  0.026  0.028  0.026  0.023  0.027  0.028  0.026


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    465 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  2  4  3  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  0  4 10 14 12  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  4  3  2  1  0  0  0     0  0  0  1  7 16 23 18 10  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  1  4 10 14 10  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  3  3  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0 -1 -1  0  0  0  0     0  0  0  0  0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0
     0  0  0  0  0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  3  8 11  8  4  2  0  0  0     0  0  0  1  6 15 16 10  4  2  1  0  0     0  0  0  1  5 10  9  5  2  1  0  0  0
     0  0  0  2 15 38 46 37 17  7  3  1  0     0  0  1  6 29 62 65 47 22  9  4  1  0     0  0  0  5 18 36 36 25 11  5  2  0  0
     0  0  0  5 26 60 74 54 29 13  5  2  0     0  0  1 10 47 94100 68 35 16  6  2  0     0  0  1  8 30 55 56 36 18  8  3  1  0
     0  0  0  3 16 37 43 31 15  7  3  1  0     0  0  0  6 28 58 61 41 18  9  3  1  0     0  0  0  4 19 35 33 22  9  5  2  0  0
     0  0  0  1  5  9 10  6  3  1  0  0  0     0  0  0  1  7 15 14  9  4  2  0  0  0     0  0  0  1  5  9  9  5  2  1  0  0  0
     0  0  0  0  1  1  1  0  0  0  0 -1  0     0  0  0  0  1  1  1  0  0  0  0 -1  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0 -1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  3  3  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  9  9  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  2  8 13 13  7  4  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  9  7  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    603 REFLECTION INTENSITIES INTEGRATED


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
   461   0  1.075  5208700    0    110      58     3  0.04558  0.09235
   462   0  1.073  5207982    0    118      66     2  0.04104  0.09578
   463   0  1.070  5207343    0    129      69     4  0.03878  0.06315
   464   0  1.070  5207368    0    129      72     6  0.04152  0.06331
   465   0  1.067  5208233    0    125      62     5  0.04304  0.08159
   466   0  1.067  5208606    0    122      59     3  0.04218  0.10230
   467   0  1.068  5208335    0    123      61     2  0.04289  0.10739
   468   0  1.069  5207862    0    124      73     1  0.03739  0.08692
   469   0  1.065  5208232    0    126      64     2  0.03626  0.07037
   470   0  1.065  5208235    0    127      65     6  0.03731  0.06864
   471   0  1.067  5208022    0    126      69     0  0.04002  0.10628
   472   0  1.067  5208118    0    115      63     0  0.04227  0.09987
   473   0  1.066  5208252    0    115      58     6  0.04160  0.06840
   474   0  1.068  5208859    0    109      53     4  0.04359  0.07785
   475   0  1.068  5209091    0     93      48     0  0.04213  0.09588
   476   0  1.069  5209566    0     98      45     2  0.04079  0.08227
   477   0  1.065  5209552    0    103      46     2  0.04264  0.08237
   478   0  1.068  5209049    0    106      49     5  0.04046  0.06023
   479   0  1.067  5208848    0    103      54     2  0.03577  0.07221
   480   0  1.067  5208901    0    101      54     6  0.03843  0.05637

   466 OUT OF    466 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.50
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.13
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.424    32.154    34.490  88.641 108.429 111.884
 REC. CELL PARAMETERS   0.041647  0.033709  0.030737  83.863  70.656  67.351
 COORDINATES OF UNIT CELL A-AXIS   -22.762   -11.152   -10.470
 COORDINATES OF UNIT CELL B-AXIS    -6.417    26.342    17.286
 COORDINATES OF UNIT CELL C-AXIS    14.269    20.232   -24.013
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.022     0.011    -0.004
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.034
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001348  0.000378
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000479  0.005239  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.29   1244.06
 DETECTOR ORIGIN (PIXELS) AT                     1265.66   1237.18
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.32
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     466 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     28     60     71     49     24     53     68     54
 SIGMAB (degree)       0.032  0.031  0.032  0.032  0.031  0.032  0.031  0.031  0.030
 SIGMAR (degree)       0.028  0.025  0.026  0.027  0.025  0.023  0.026  0.028  0.027


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    466 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  1  1  2  4  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  2  2  1  0  0  0  0     0  0  0  1  4 12 16 13  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  1     0  0  0  0  1  3  4  4  2  1  0  0  0     0  0  0  1  7 20 28 22 12  6  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  1  5 13 17 13  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  3  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0 -1 -1  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     1  0  0  0  0  0  1  1  0  0 -1  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  3  8 12  7  4  2  1  0  0     0  0  0  1  5 14 17 10  4  2  1  0  0     0  0  0  1  4  8  8  5  2  1  0  0  0
     0  0  0  3 17 44 52 36 20  8  4  1  0     0  0  1  7 30 66 66 43 23  9  4  1  0     0  0  1  4 17 35 32 20 10  4  2  0  0
     0  0  0  5 27 65 80 59 31 14  6  2  0     0  0  1 11 47 96100 69 34 15  6  2  1     0  0  1  7 27 50 48 32 15  7  3  1  0
     0  0  0  3 19 41 50 34 17  7  3  1  0     0  0  1  7 31 59 63 41 19  9  3  1  0     0  0  0  4 18 30 31 20  9  4  2  0  0
     0  0  0  1  5  9 10  5  3  1  0  0  0     0  0  0  2  8 14 13  7  4  2  0  0  0     0  0  0  2  5  8  7  4  2  1  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0 -1 -1  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0    -1  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  8  7  4  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  2  6 11 10  7  3  1  0  0  0     0  0  0  1  1  1  1  1  0  0  0  0  0     0  0  0  1  1  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  6  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    594 REFLECTION INTENSITIES INTEGRATED


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
   481   0  1.065  5208573    0    112      57     4  0.04091  0.07873
   482   0  1.066  5208241    0    111      58     4  0.03602  0.06341
   483   0  1.067  5208538    0    122      55     5  0.03718  0.04716
   484   0  1.065  5207988    0    133      65     2  0.04208  0.09692
   485   0  1.061  5207791    0    120      62     4  0.03813  0.07963
   486   0  1.063  5207959    0    119      64     4  0.04070  0.08846
   487   0  1.064  5208714    0    124      57     4  0.04075  0.08396
   488   0  1.064  5208266    0    111      64     5  0.04143  0.06973
   489   0  1.063  5208079    0    118      63     4  0.03983  0.08468
   490   0  1.064  5209123    0    118      51     4  0.04568  0.07756
   491   0  1.061  5208311    0    117      55     4  0.04219  0.07888
   492   0  1.063  5207345    0    129      72     8  0.03914  0.05870
   493   0  1.063  5208419    0    136      63     3  0.04375  0.08079
   494   0  1.062  5208238    0    126      61     4  0.03795  0.07868
   495   0  1.061  5208017    0    124      64     4  0.03993  0.09307
   496   0  1.063  5208261    0    125      60     0  0.04256  0.11348
   497   0  1.064  5207668    0    122      64     3  0.03996  0.07457
   498   0  1.062  5207856    0    134      65     4  0.04163  0.06553
   499   0  1.061  5208103    0    131      60     3  0.03909  0.06991
   500   0  1.063  5208165    0    121      59     1  0.04022  0.09805

   482 OUT OF    482 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.53
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.415    32.138    34.493  88.650 108.426 111.879
 REC. CELL PARAMETERS   0.041660  0.033725  0.030735  83.856  70.656  67.354
 COORDINATES OF UNIT CELL A-AXIS   -22.755   -11.148   -10.463
 COORDINATES OF UNIT CELL B-AXIS    -6.415    26.327    17.279
 COORDINATES OF UNIT CELL C-AXIS    14.266    20.232   -24.018
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.020     0.008    -0.004
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.034
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001348  0.000378
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000639  0.005023  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.28   1244.02
 DETECTOR ORIGIN (PIXELS) AT                     1265.44   1237.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.21
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     482 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   49     28     57     65     59     30     56     71     57
 SIGMAB (degree)       0.032  0.031  0.031  0.031  0.031  0.032  0.031  0.031  0.030
 SIGMAR (degree)       0.027  0.024  0.026  0.027  0.025  0.023  0.026  0.027  0.026


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    482 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  3  7  9  8  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  1  0  0  0     0  0  0  1  5 13 16 13  7  3  1  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  1  3  8 10  7  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  1  1  0  0 -1 -1  0     0  0  0  0  0  0  1  0  0  0  0 -1  0
     0  0  0  1  2  6  7  5  2  1  0  0  0     0  0  0  1  5 11 13  9  4  1  1  0  0     0  0  0  1  5  8  8  5  2  1  0  0  0
     0  0  0  2 11 29 37 27 14  6  2  1  0     0  0  0  5 25 57 63 44 22 10  4  1  0     0  0  1  5 19 39 41 26 12  5  2  1  0
     0  0  0  4 20 50 58 45 24 11  4  1  0     0  0  1  9 42 92100 70 36 16  6  2  0     0  0  1  8 32 62 62 40 20  9  3  1  0
     0  0  0  3 12 32 36 25 14  6  2  1  0     0  0  1  6 25 57 62 40 21  8  3  1  0     0  0  1  5 20 38 39 22 11  5  2  0  0
     0  0  0  1  3  8  7  4  2  1  0  0  0     0  0  0  2  7 14 11  7  3  1  0  0  0     0  0  0  2  5  9  7  4  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5 10 10  6  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  8 16 15  9  5  2  1  0  0     0  0  0  0  1  2  1  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 10 10  6  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    629 REFLECTION INTENSITIES INTEGRATED


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
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  2  3  3  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  2  0  0  0  0  0     0  0  0  0  5 12 14 10  5  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  3  3  3  0  0  0  0  0     0  0  0  2  8 17 21 15  8  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  0  5 11 13  9  5  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  5 10 11  7  3  0  0  0  0     0  0  0  3  8 16 17 10  4  0  0  0  0     0  0  0  2  6 10  9  5  2  0  0  0  0
     0  0  0  4 20 40 47 32 15  6  0  0  0     0  0  0  9 35 65 69 44 20  8  0  0  0     0  0  0  6 23 39 38 23 10  4  0  0  0
     0  0  0  7 29 59 69 50 24 10  0  0  0     0  0  2 13 51 94100 69 32 13  0  0  0     0  0  0  9 33 56 54 36 16  6  0  0  0
     0  0  0  4 18 39 45 31 15  6  0  0  0     0  0  0  7 31 63 67 44 20  8  0  0  0     0  0  0  5 20 36 37 23 10  4  0  0  0
     0  0  0  0  4 10 11  7  3  0  0  0  0     0  0  0  2  9 17 18 10  4  0  0  0  0     0  0  0  0  6 10 10  5  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  7 10  9  5  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  9 14 13  8  3  0  0  0  0     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  5  9  9  5  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  2 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  6  7  5  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  8 10  8  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  4  6  4  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  0  0  0  0  0  0     0  0  0  0  0  3  3  0  0  0  0  0  0
     0  0  0  0  3  6  8  5  3  0  0  0  0     0  0  0  3  8 17 17 10  5  0  0  0  0     0  0  0  3  8 13 11  6  3  0  0  0  0
     0  0  0  3 11 26 32 22 11  5  0  0  0     0  0  0  7 33 69 71 44 22  9  0  0  0     0  0  0  8 32 53 46 27 13  5  0  0  0
     0  0  0  3 15 35 44 33 17  7  0  0  0     0  0  2  9 46 93100 67 32 14  5  0  0     0  0  2 11 45 75 68 41 19  8  0  0  0
     0  0  0  2  8 21 27 19 10  4  0  0  0     0  0  0  4 25 55 60 39 19  8  0  0  0     0  0  0  6 26 44 41 24 11  5  0  0  0
     0  0  0  0  3  5  6  5  3  0  0  0  0     0  0  0  0  6 13 14 10  4  0  0  0  0     0  0  0  2  7 11 10  6  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  3  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  7 10  8  5  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4 10 15 12  7  3  2  0  0  0     0  0  0  0  0  2  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  6  9  8  5  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  3  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  3 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  3  8 10  7  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  0  6 14 17 13  7  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  2  2  0  0  0  0  0     0  0  0  0  4 10 14 10  5  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  2  3  4  3  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  6  7  5  2  0  0  0  0     0  0  0  0  5 10 11  6  3  0  0  0  0     0  0  0  0  4  7  6  3  2  0  0  0  0
     0  0  0  2 12 28 34 23 11  5  0  0  0     0  0  0  5 23 48 52 33 16  7  0  0  0     0  0  0  4 17 32 31 18  9  4  0  0  0
     0  0  0  4 22 51 63 45 23 11  0  0  0     0  0  0  9 43 90100 68 34 15  6  0  0     0  0  0  8 34 62 62 39 19  8  0  0  0
     0  0  0  3 16 40 49 35 19  8  0  0  0     0  0  0  7 33 71 80 54 28 12  0  0  0     0  0  0  7 27 50 50 33 16  7  0  0  0
     0  0  0  0  5 11 13 10  5  3  0  0  0     0  0  0  3 11 21 23 15  8  3  0  0  0     0  0  0  3  9 15 15 10  4  3  0  0  0
     0  0  0  0  0  2  3  0  0  0  0  0  0     0  0  0  0  3  3  3  2  0  0  0  0  0     0  0  0  0  0  3  3  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  2  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4  8  7  4  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  9 16 15  9  4  2  0  0  0     0  0  0  0  0  2  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  8 14 12  8  4  2  0  0  0     0  0  0  0  0  2  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  4  4  4  3  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  4 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  2  4  5  4  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  2  0  0  0  0     0  0  0  0  6 15 20 15  8  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  2  3  4  4  3  0  0  0  0     0  0  0  0  9 21 27 21 11  5  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  2  0  0  0  0     0  0  0  0  5 13 16 11  6  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  5 12 15 10  5  2  0  0  0     0  0  0  2  8 17 19 11  5  3  0  0  0     0  0  0  0  5 10 10  6  3  0  0  0  0
     0  0  0  4 20 47 58 41 21  9  0  0  0     0  0  0  6 31 65 72 47 24 10  0  0  0     0  0  0  4 19 37 38 24 12  5  0  0  0
     0  0  0  6 29 67 80 58 30 14  0  0  0     0  0  0 10 45 92100 68 34 16  7  0  0     0  0  0  7 28 53 53 35 17  8  0  0  0
     0  0  0  3 16 39 47 32 16  7  0  0  0     0  0  0  6 25 54 60 38 19  8  0  0  0     0  0  0  4 16 31 32 20 10  4  0  0  0
     0  0  0  0  4  8  9  6  3  0  0  0  0     0  0  0  0  6 11 12  8  4  0  0  0  0     0  0  0  0  4  7  7  4  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  5 10  9  6  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  8 14 13  8  4  2  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  5  8  8  5  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  5 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  2  0  0  0  0     0  0  0  0  4 10 14 11  6  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  4  4  3  0  0  0  0     0  0  0  2  7 17 23 18 10  5  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  0  4 10 13  9  5  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  8 10  7  3  0  0  0  0     0  0  0  3  6 14 15  9  4  2  0  0  0     0  0  0  0  4  8  7  4  2  0  0  0  0
     0  0  0  3 17 38 48 35 18  8  0  0  0     0  0  0  6 30 61 67 45 22 10  0  0  0     0  0  0  4 17 32 33 21 11  5  0  0  0
     0  0  0  6 27 62 75 54 28 13  6  0  0     0  0  2 11 46 93100 67 34 16  6  0  0     0  0  0  7 26 48 48 31 16  7  0  0  0
     0  0  0  4 15 35 42 28 14  6  0  0  0     0  0  0  6 24 49 53 33 16  7  0  0  0     0  0  0  4 13 25 25 15  7  3  0  0  0
     0  0  0  0  3  6  6  4  2  0  0  0  0     0  0  0  0  5  9  8  5  2  0  0  0  0     0  0  0  0  3  4  4  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  4  7  7  4  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 10  9  6  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  5  5  3  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  6 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  0  0  0  0  0  0  0     0  0  0  0  3  6  8  6  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  2  4 10 12  9  5  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  0  0  0  0  0  0     0  0  0  0  3  6  7  5  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  6  7  5  3  0  0  0  0     0  0  0  2  6 14 14  8  4  2  0  0  0     0  0  0  0  4  7  6  4  2  0  0  0  0
     0  0  0  3 13 31 36 24 12  5  0  0  0     0  0  0  8 35 65 64 39 19  8  0  0  0     0  0  0  5 18 31 29 17  8  4  0  0  0
     0  0  0  5 24 49 55 37 19  8  0  0  0     0  0  3 15 57100 96 60 29 13  0  0  0     0  0  2  9 29 46 43 26 13  5  0  0  0
     0  0  0  4 14 30 32 21 10  4  0  0  0     0  0  3  9 33 59 56 34 16  7  0  0  0     0  0  0  5 16 27 25 15  7  3  0  0  0
     0  0  0  0  3  7  6  4  3  0  0  0  0     0  0  0  3  8 14 12  7  3  0  0  0  0     0  0  0  0  4  6  5  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  5  4  3  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  7  6  4  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  4  4  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  7 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  2  0  0  0  0  0     0  0  0  0  3  7 10  8  4  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  2  0  0  0  0     0  0  0  0  5 15 20 16  8  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  2  0  0  0  0  0     0  0  0  0  4 11 14 10  5  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  4  6  4  3  0  0  0  0     0  0  0  0  3  7  9  6  3  0  0  0  0     0  0  0  0  3  4  4  3  2  0  0  0  0
     0  0  0  0  9 28 38 27 14  6  0  0  0     0  0  0  3 19 46 55 36 18  8  0  0  0     0  0  0  3 11 24 27 17  8  4  0  0  0
     0  0  0  3 23 58 71 52 27 12  0  0  0     0  0  0  8 43 92100 68 34 15  6  0  0     0  0  0  6 25 49 50 33 16  7  0  0  0
     0  0  0  3 18 42 47 32 16  7  0  0  0     0  0  0  6 32 65 66 43 20  9  0  0  0     0  0  0  4 19 36 34 22 10  4  0  0  0
     0  0  0  0  4  9  9  5  3  0  0  0  0     0  0  0  2  8 15 14  7  3  0  0  0  0     0  0  0  0  5  8  7  4  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  5  5  3  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 11 11  7  3  2  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  5  8  8  5  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  8 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  2  0  0  0  0     0  0  0  0  3 10 14 12  6  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  4  4  3  0  0  0  0     0  0  0  0  8 20 27 21 12  6  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  3  0  0  0  0     0  0  0  0  6 14 17 13  7  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  2  3  4  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  5  7  5  3  0  0  0  0     0  0  0  0  4  8 10  6  3  0  0  0  0     0  0  0  0  3  5  5  3  2  0  0  0  0
     0  0  0  3 13 34 44 33 17  8  0  0  0     0  0  0  5 21 49 57 40 20  9  0  0  0     0  0  0  4 14 29 31 20 10  4  0  0  0
     0  0  0  5 27 65 80 60 32 15  8  0  0     0  0  0  9 43 91100 71 36 17  7  0  0     0  0  0  6 27 51 53 36 18  8  0  0  0
     0  0  0  4 21 46 53 37 18  8  0  0  0     0  0  0  7 32 64 66 44 21  9  0  0  0     0  0  0  5 19 36 36 23 11  5  0  0  0
     0  0  0  2  5 10 10  6  3  0  0  0  0     0  0  0  2  8 14 13  7  3  0  0  0  0     0  0  0  0  5  8  7  4  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4  8  8  5  3  0  0  0  0     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  8 13 13  9  4  3  0  0  0     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  6  9  9  5  3  0  0  0  0     0  0  0  0  0  3  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  9 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  0  0  0  0  0     0  0  0  0  4 11 14 11  6  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  4  3  0  0  0  0     0  0  0  0  6 15 19 15  8  5  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  2  0  0  0  0  0     0  0  0  0  3  7  9  7  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  4  9 11  7  4  0  0  0  0     0  0  0  3  9 17 18 11  5  3  0  0  0     0  0  0  3  7 12 11  6  3  0  0  0  0
     0  0  0  4 18 42 51 37 19  8  0  0  0     0  0  0  9 36 75 80 53 27 11  0  0  0     0  0  0  8 28 50 47 29 14  6  0  0  0
     0  0  0  4 24 55 66 49 25 11  0  0  0     0  0  2 10 47 94100 68 34 15  6  0  0     0  0  2  9 35 61 58 37 18  8  0  0  0
     0  0  0  3 10 25 29 20 10  5  0  0  0     0  0  0  4 19 42 44 28 14  6  0  0  0     0  0  0  4 14 27 26 15  7  3  0  0  0
     0  0  0  0  2  3  4  3  0  0  0  0  0     0  0  0  0  3  6  7  4  2  0  0  0  0     0  0  0  0  3  4  4  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  4  3  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  8 12 11  6  3  2  0  0  0     0  0  0  0  0  2  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  9 15 13  8  4  2  0  0  0     0  0  0  0  3  3  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  4  7  6  3  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 STANDARD DEVIATION (DEGREES) OF BEAM DIVERGENCE
 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1    0.023  0.023  0.022  0.022  0.022  0.023  0.022  0.023  0.023
     2    0.022  0.022  0.022  0.022  0.022  0.022  0.022  0.022  0.023
     3    0.023  0.023  0.023  0.022  0.022  0.023  0.022  0.023  0.023
     4    0.023  0.023  0.022  0.022  0.022  0.022  0.021  0.022  0.022
     5    0.022  0.023  0.023  0.023  0.022  0.022  0.021  0.022  0.022
     6    0.023  0.023  0.023  0.023  0.022  0.023  0.022  0.022  0.023
     7    0.023  0.023  0.023  0.022  0.022  0.023  0.022  0.022  0.023
     8    0.022  0.022  0.022  0.022  0.022  0.022  0.022  0.022  0.022
     9    0.022  0.023  0.023  0.022  0.022  0.023  0.021  0.022  0.022
    10    0.023  0.023  0.022  0.022  0.022  0.022  0.022  0.022  0.022
    11    0.023  0.023  0.023  0.023  0.022  0.022  0.022  0.022  0.022
    12    0.022  0.022  0.022  0.022  0.022  0.022  0.022  0.022  0.022
    13    0.022  0.023  0.023  0.022  0.022  0.022  0.022  0.022  0.023
    14    0.022  0.022  0.023  0.023  0.022  0.023  0.022  0.022  0.022
    15    0.022  0.022  0.023  0.023  0.022  0.023  0.022  0.022  0.022
    16    0.022  0.022  0.023  0.022  0.022  0.022  0.022  0.022  0.022
    17    0.022  0.022  0.023  0.022  0.022  0.023  0.022  0.022  0.022
    18    0.023  0.023  0.023  0.022  0.022  0.023  0.022  0.022  0.022
    19    0.022  0.022  0.023  0.023  0.022  0.022  0.022  0.022  0.022
    20    0.022  0.023  0.023  0.022  0.022  0.022  0.022  0.023  0.022
    21    0.023  0.022  0.023  0.022  0.022  0.023  0.022  0.022  0.022
    22    0.022  0.023  0.023  0.022  0.022  0.022  0.022  0.022  0.022
    23    0.023  0.023  0.023  0.022  0.022  0.022  0.022  0.022  0.022
    24    0.022  0.022  0.023  0.023  0.022  0.023  0.022  0.022  0.022
    25    0.022  0.022  0.022  0.022  0.022  0.023  0.022  0.022  0.022


 STANDARD DEVIATION (DEGREES) OF REFLECTING RANGE
 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1    0.024  0.023  0.025  0.026  0.024  0.020  0.024  0.027  0.026
     2    0.025  0.022  0.026  0.027  0.024  0.021  0.025  0.028  0.025
     3    0.027  0.029  0.028  0.027  0.025  0.021  0.025  0.028  0.027
     4    0.025  0.024  0.026  0.030  0.028  0.021  0.024  0.026  0.026
     5    0.026  0.024  0.027  0.027  0.024  0.021  0.026  0.028  0.025
     6    0.026  0.030  0.027  0.026  0.024  0.022  0.025  0.027  0.028
     7    0.025  0.021  0.026  0.028  0.025  0.021  0.025  0.028  0.025
     8    0.025  0.026  0.027  0.027  0.024  0.020  0.024  0.027  0.026
     9    0.030  0.027  0.027  0.031  0.033  0.031  0.030  0.028  0.027
    10    0.026  0.025  0.026  0.027  0.024  0.020  0.024  0.027  0.026
    11    0.025  0.022  0.025  0.026  0.024  0.023  0.026  0.028  0.026
    12    0.027  0.025  0.027  0.030  0.027  0.022  0.025  0.027  0.026
    13    0.028  0.024  0.027  0.033  0.033  0.024  0.025  0.027  0.026
    14    0.027  0.024  0.026  0.027  0.024  0.021  0.028  0.031  0.028
    15    0.025  0.024  0.026  0.027  0.025  0.021  0.025  0.026  0.026
    16    0.026  0.022  0.026  0.028  0.025  0.022  0.025  0.026  0.024
    17    0.025  0.021  0.024  0.027  0.024  0.021  0.024  0.027  0.024
    18    0.025  0.024  0.026  0.028  0.026  0.020  0.025  0.027  0.025
    19    0.025  0.023  0.025  0.028  0.024  0.019  0.024  0.028  0.026
    20    0.027  0.023  0.027  0.031  0.028  0.023  0.026  0.029  0.026
    21    0.025  0.022  0.025  0.027  0.024  0.020  0.024  0.027  0.025
    22    0.026  0.028  0.027  0.028  0.025  0.020  0.024  0.026  0.026
    23    0.026  0.026  0.026  0.029  0.027  0.022  0.026  0.027  0.026
    24    0.027  0.022  0.025  0.027  0.025  0.021  0.027  0.032  0.029
    25    0.026  0.025  0.026  0.026  0.025  0.023  0.025  0.026  0.026


 ***** SUGGESTED VALUES FOR INPUT PARAMETERS *****
 BEAM_DIVERGENCE=   0.22303  BEAM_DIVERGENCE_E.S.D.=   0.02230
 REFLECTING_RANGE=  0.17891  REFLECTING_RANGE_E.S.D.=  0.02556



 NUMBER OF STRONG REFLECTIONS USED FOR LEARNING REFERENCE PROFILES

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1       51     29     59     71     59     31     45     65     57
     2       50     30     56     67     52     29     62     67     57
     3       46     33     58     64     57     28     50     69     57
     4       49     34     57     64     53     28     56     63     57
     5       47     28     49     64     55     27     52     66     57
     6       45     31     57     65     49     26     52     67     56
     7       49     32     56     69     59     28     48     60     58
     8       52     30     55     71     55     34     61     68     58
     9       48     32     60     61     44     22     45     63     56
    10       49     28     54     69     56     30     54     66     55
    11       51     33     53     66     60     32     57     71     58
    12       47     31     61     67     55     26     47     65     55
    13       49     32     59     67     51     29     52     66     54
    14       47     32     56     65     49     28     53     67     56
    15       46     32     58     67     52     26     50     66     60
    16       49     31     59     70     53     27     53     65     55
    17       47     30     51     64     53     34     55     64     55
    18       53     30     57     67     53     29     57     66     58
    19       46     28     60     70     55     30     52     65     48
    20       51     32     57     69     56     30     54     65     58
    21       48     31     63     68     50     29     52     65     55
    22       48     31     54     60     57     29     50     60     53
    23       47     30     51     63     50     29     57     69     62
    24       47     28     60     71     49     24     53     68     54
    25       49     28     57     65     59     30     56     71     57
 margin    1211    765   1419   1664   1341    714   1325   1647   1404



 NUMBER OF GRID POINTS IN REFERENCE PROFILES WITH  SIGNAL >   2.% OF MAXIMUM

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1      134    121    126    131    121    110    120    135    129
     2      129    114    126    132    116    112    118    138    125
     3      139    125    144    133    120    121    124    142    133
     4      137    118    130    141    129    109    117    126    126
     5      133    120    139    137    119    113    119    129    122
     6      137    138    144    135    119    114    120    132    134
     7      135    120    137    140    125    112    115    131    123
     8      130    125    131    129    116    115    118    131    125
     9      148    126    137    138    133    129    128    134    122
    10      135    126    125    135    120    112    117    134    120
    11      131    117    126    136    123    112    124    132    120
    12      138    120    130    142    127    114    120    131    119
    13      145    130    152    146    131    122    123    128    127
    14      137    120    138    139    120    111    127    137    130
    15      134    115    136    136    128    116    122    127    122
    16      131    114    136    142    124    114    119    125    117
    17      135    116    133    136    119    112    118    129    117
    18      140    123    133    140    126    114    115    127    123
    19      131    113    131    141    122    110    119    132    124
    20      139    121    140    146    136    118    125    139    129
    21      134    117    140    137    123    116    122    135    120
    22      130    119    132    138    126    113    115    128    124
    23      144    125    141    148    133    117    124    126    122
    24      138    120    131    138    122    117    129    137    125
    25      130    116    129    136    121    117    118    124    115
 margin     160    147    170    160    150    141    145    161    155

 MINIMUM INTEGRATION VOLUME IS  109 GRID POINTS



 CORRELATIONS BETWEEN REFERENCE PROFILES WITH CENTRAL PROFILE #1

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1     1.00   0.98   0.97   0.98   0.98   0.97   0.97   0.97   0.95
     2     1.00   0.96   0.97   0.99   0.98   0.97   0.98   0.96   0.95
     3     1.00   0.85   0.91   0.98   0.96   0.96   0.97   0.98   0.95
     4     1.00   0.91   0.94   0.96   0.95   0.98   0.97   0.97   0.95
     5     1.00   0.93   0.92   0.99   0.98   0.98   0.97   0.98   0.95
     6     1.00   0.90   0.93   0.99   0.98   0.97   0.98   0.99   0.95
     7     1.00   0.96   0.97   0.97   0.97   0.97   0.98   0.97   0.95
     8     1.00   0.93   0.96   0.98   0.98   0.98   0.98   0.97   0.96
     9     1.00   0.85   0.95   0.98   0.93   0.90   0.95   0.97   0.91
    10     1.00   0.94   0.95   0.98   0.99   0.98   0.96   0.97   0.96
    11     1.00   0.97   0.98   0.99   0.98   0.94   0.98   0.96   0.95
    12     1.00   0.95   0.97   0.95   0.97   0.98   0.97   0.98   0.95
    13     1.00   0.91   0.96   0.96   0.90   0.93   0.98   0.99   0.90
    14     1.00   0.90   0.93   0.99   0.98   0.96   0.95   0.94   0.96
    15     1.00   0.95   0.96   0.99   0.99   0.97   0.96   0.96   0.96
    16     1.00   0.97   0.97   0.96   0.97   0.98   0.98   0.97   0.96
    17     1.00   0.99   0.97   0.97   0.97   0.95   0.97   0.95   0.95
    18     1.00   0.92   0.95   0.98   0.97   0.98   0.97   0.97   0.97
    19     1.00   0.95   0.98   0.97   0.97   0.96   0.96   0.94   0.95
    20     1.00   0.96   0.97   0.97   0.97   0.98   0.97   0.97   0.96
    21     1.00   0.97   0.94   0.98   0.97   0.97   0.97   0.98   0.94
    22     1.00   0.84   0.95   0.98   0.97   0.98   0.98   0.98   0.91
    23     1.00   0.83   0.93   0.98   0.97   0.97   0.97   0.97   0.92
    24     1.00   0.97   0.97   0.98   0.98   0.97   0.96   0.92   0.93
    25     1.00   0.92   0.97   0.98   0.99   0.98   0.98   0.98   0.96
 margin    1.00   0.95   0.97   0.99   0.99   0.98   0.98   0.98   0.97



 CORRELATIONS BETWEEN REFERENCE PROFILES WITH RUN-AVERAGED PROFILES

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1     0.99   0.99   1.00   1.00   0.99   0.98   0.99   1.00   0.99
     2     1.00   0.99   1.00   0.99   0.99   0.98   1.00   0.99   1.00
     3     0.99   0.95   0.98   0.99   0.99   0.99   1.00   0.99   0.99
     4     0.99   0.99   0.99   0.99   0.98   0.99   0.99   0.99   0.99
     5     1.00   1.00   0.99   0.99   1.00   0.99   0.99   0.99   0.99
     6     0.99   0.95   0.96   0.98   0.98   0.96   0.95   0.96   0.95
     7     1.00   0.98   0.99   1.00   1.00   1.00   0.99   0.99   0.99
     8     0.99   0.99   0.99   0.99   0.99   0.99   1.00   0.99   1.00
     9     0.98   0.98   0.99   0.98   0.92   0.82   0.96   0.99   0.99
    10     0.99   0.99   1.00   1.00   0.99   1.00   1.00   1.00   1.00
    11     0.99   0.99   0.99   0.98   0.96   0.92   0.98   0.99   0.99
    12     0.99   1.00   0.99   0.98   0.99   0.99   0.97   0.98   0.99
    13     0.99   1.00   0.98   0.96   0.92   0.91   0.99   0.99   0.98
    14     0.99   0.99   0.99   0.99   0.99   0.99   0.97   0.97   0.98
    15     1.00   1.00   0.99   0.99   0.99   1.00   0.99   0.99   1.00
    16     0.99   0.98   0.99   1.00   1.00   0.99   0.99   0.99   0.99
    17     1.00   0.94   0.99   0.99   0.98   0.97   0.99   0.99   0.97
    18     1.00   1.00   1.00   1.00   0.99   0.99   1.00   1.00   1.00
    19     0.99   0.99   0.99   1.00   1.00   0.99   0.99   0.99   0.99
    20     0.98   0.96   0.97   0.98   0.97   0.96   0.98   0.99   0.98
    21     0.99   0.98   0.99   1.00   0.99   0.99   1.00   0.99   0.99
    22     0.99   0.95   0.99   1.00   1.00   1.00   0.98   0.98   0.96
    23     0.99   0.97   0.99   1.00   0.99   0.97   0.99   0.99   0.99
    24     0.99   0.93   0.99   0.99   0.98   0.99   0.97   0.95   0.92
    25     0.99   0.98   0.99   1.00   0.99   0.98   0.97   0.97   0.97
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

       -1... 0      -0.120      0.054        -1.         12
        0... 1       0.273      0.110         4.        109
        1... 2       0.502      0.102        13.        228
        2... 3       0.647      0.087        27.        255
        3... 6       0.717      0.093        77.       1124
        6... 9       0.754      0.090       185.       1456
        9...12       0.764      0.091       337.       1494
       12...15       0.772      0.087       543.       1269
       15...18       0.778      0.087       786.       1054
       18...21       0.784      0.086      1046.        878
       21...24       0.794      0.086      1381.        675
       24...27       0.789      0.092      1787.        625
       27...30       0.782      0.087      2306.        562
       30...33       0.788      0.090      2938.        484
       33...36       0.786      0.094      3817.        465
       36...39       0.793      0.090      5025.        452
       39...42       0.790      0.094      6841.        461
       42...45       0.796      0.092     10410.        518
       45...48       0.790      0.102     17278.        631
       48...51       0.797      0.116     29730.        373

   13125 REFLECTIONS SAVED ON FILE "INTEGRATE.HKL"


 total elapsed wall-clock time for INTEGRATE      19.6 sec

 TASK     cpu time (sec)    elapsed wall-clock time (sec)
   1         296.6                    19.4
