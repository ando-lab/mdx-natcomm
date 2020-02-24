
 ***** INTEGRATE ***** (VERSION Mar 15, 2019  BUILT=20190806)   6-Jan-2020

 INPUT PARAMETER VALUES
 ----------------------
 MAXIMUM_NUMBER_OF_PROCESSORS=   0  (0: automatic choice of OpenMP threads)
 number of OpenMP threads used  20
 MAXIMUM_NUMBER_OF_JOBS=   0  (0:automatic choice of forked main tasks)
 number of forked INTEGRATE tasks   1
 CLUSTER_NODES=
 SECONDS=     0
 NAME_TEMPLATE_OF_DATA_FRAMES=~/data/lysozyme/nitrate/ubatch4/lys_nitr_10_1_????.cbf CBF
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
 UNIT_CELL_CONSTANTS=    27.43    32.13    34.50  88.676 108.456 111.864
 ROTATION_AXIS= -0.99999  0.00182  0.00312
 OSCILLATION_RANGE=  0.100000 DEGREES
 STARTING_ANGLE= 250.000 STARTING_FRAME=       1
 X-RAY_WAVELENGTH=   0.97700 ANGSTROM
 INCIDENT_BEAM_DIRECTION=  0.000779  0.005206  1.023528
 NUMBER OF DETECTOR SEGMENTS   1
 NX=  2463 NY=  2527 QX= 0.17200 QY= 0.17200
 ORGX=   1265.36  ORGY=   1237.28
 DETECTOR_DISTANCE=   231.202
 DIRECTION_OF_DETECTOR_X-AXIS=   1.00000   0.00000   0.00000
 DIRECTION_OF_DETECTOR_Y-AXIS=   0.00000   1.00000   0.00000
 SEGMENT=    1 2463    1 2527  ! x1,x2,y1,y2 (pixel units)
    SEGMENT_ORGX=      0.00   SEGMENT_ORGY=      0.00
    SEGMENT_DISTANCE=     0.00
    DIRECTION_OF_SEGMENT_X-AXIS=  1.000000  0.000000  0.000000
    DIRECTION_OF_SEGMENT_Y-AXIS=  0.000000  1.000000  0.000000
 RESOLUTION RANGE RECORDED BY DETECTOR (ANGSTROM)    44.470     1.211
 NUMBER OF TRUSTED DETECTOR PIXELS    5186964
 MEAN CONTENTS OF TRUSTED PIXELS IN BACKGROUND TABLE    139.489
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
     1   0  1.036  5180764    0    107      64     3  0.04629  0.08288
     2   0  0.997  5179991    0    120      73     8  0.04452  0.06162
     3   0  0.998  5180029    0    129      73     4  0.04255  0.08063
     4   0  0.994  5179819    0    131      74     3  0.04411  0.08742
     5   0  0.996  5179906    0    129      72     5  0.04260  0.06609
     6   0  0.999  5180041    0    132      74     5  0.04544  0.08126
     7   0  0.997  5180115    0    128      67     1  0.04407  0.08708
     8   0  0.999  5179195    0    136      83     3  0.04489  0.07621
     9   0  0.997  5178235    0    144      94    10  0.04056  0.05814
    10   0  0.998  5179085    0    154      88     6  0.04213  0.08194
    11   0  0.998  5178973    0    153      79     2  0.04595  0.09278
    12   0  0.999  5179193    0    141      80     2  0.04109  0.09784
    13   0  0.998  5179790    0    133      72     2  0.04523  0.08409
    14   0  0.999  5180371    0    127      72     4  0.04199  0.08853
    15   0  1.001  5180316    0    126      64     0  0.04608  0.10542
    16   0  0.997  5179822    0    147      74     2  0.04441  0.09349
    17   0  1.000  5179339    0    137      81     6  0.03979  0.06968
    18   0  0.999  5179844    0    134      76     1  0.04039  0.09030
    19   0  1.002  5181115    0    124      59     1  0.04421  0.08578
    20   0  1.001  5180726    0    114      65     6  0.04315  0.07149

   514 OUT OF    515 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.53
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.429    32.135    34.492  88.692 108.441 111.867
 REC. CELL PARAMETERS   0.041643  0.033728  0.030741  83.806  70.624  67.350
 COORDINATES OF UNIT CELL A-AXIS   -14.256   -14.237    18.612
 COORDINATES OF UNIT CELL B-AXIS    14.447    25.591    13.002
 COORDINATES OF UNIT CELL C-AXIS   -21.606    21.433   -16.233
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.018     0.017    -0.014
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.035
 LAB COORDINATES OF ROTATION AXIS -0.999993  0.001823  0.003120
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000770  0.005169  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.34   1244.01
 DETECTOR ORIGIN (PIXELS) AT                     1265.33   1237.22
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.25
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     514 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   51     36     63     70     59     34     58     69     63
 SIGMAB (degree)       0.033  0.034  0.034  0.033  0.032  0.032  0.032  0.033  0.033
 SIGMAR (degree)       0.028  0.024  0.027  0.028  0.026  0.023  0.026  0.028  0.026


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    514 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  3  3  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  1  5 12 16 12  7  3  1  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  3  3  2  1  0  0  0     0  0  0  2  7 17 23 19 12  6  2  1  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  1  4 11 13 11  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  3  4  3  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  4  9 10  6  3  1  0  0  0     0  0  0  2  6 12 13  8  4  2  0  0  0     0  0  0  1  3  6  7  4  2  1  0  0  0
     0  0  1  5 19 40 49 35 17  8  3  1  0     0  0  1  8 32 60 65 43 20  9  4  1  0     0  0  1  5 17 31 32 21  9  4  2  0  0
     0  0  1  7 29 60 73 58 31 13  5  2  0     0  0  2 12 51 94100 72 36 16  6  2  1     0  0  1  8 29 50 50 34 17  7  3  1  0
     0  0  1  5 17 38 43 32 17  7  3  1  0     0  0  1  8 29 58 61 40 20  9  3  1  0     0  0  1  5 15 30 31 20 10  4  1  0  0
     0  0  0  2  4  9  9  6  3  1  0  0  0     0  0  0  3  7 13 13  8  4  2  0  0  0     0  0  0  2  4  7  6  4  2  1  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0 -1  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  6  6  4  1  1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  0  1  0  0  0  0  0  0
     0  0  0  2  6 10  9  6  3  1  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  1  0  1  1  0  0  0  0  0
     0  0  0  1  3  6  6  4  2  1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0
     1  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0

    636 REFLECTION INTENSITIES INTEGRATED


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
    21   0  1.005  5181647    0    123      56     5  0.04766  0.08918
    22   0  1.006  5180471    0    125      69     3  0.04480  0.08465
    23   0  1.001  5179982    0    127      72     2  0.03856  0.07216
    24   0  1.000  5180364    0    124      67     2  0.04017  0.05738
    25   0  1.005  5180129    0    129      70     1  0.04475  0.07370
    26   0  1.002  5181459    0    125      59     1  0.04662  0.11130
    27   0  1.005  5181105    0    118      63     3  0.04568  0.08593
    28   0  1.001  5181349    0    122      57     2  0.04323  0.07295
    29   0  1.002  5180186    0    125      72     6  0.04144  0.04844
    30   0  1.000  5180450    0    129      73     5  0.04145  0.06791
    31   0  1.004  5180902    0    128      63     5  0.04099  0.06486
    32   0  1.004  5179829    0    116      75     2  0.04588  0.08356
    33   0  1.003  5179954    0    122      72     3  0.04464  0.08349
    34   0  1.005  5180554    0    127      65     4  0.04552  0.08092
    35   0  1.002  5180636    0    121      59     6  0.04109  0.06213
    36   0  1.006  5180261    0    116      71     1  0.04272  0.10715
    37   0  1.004  5180891    0    122      63     2  0.04309  0.08476
    38   0  1.004  5180915    0    111      64     0  0.04131  0.09563
    39   0  1.001  5180672    0    110      65     1  0.04172  0.08252
    40   0  1.003  5180748    0    117      66     0  0.04202  0.07755

   463 OUT OF    466 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.52
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.04
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.434    32.141    34.512  88.643 108.469 111.886
 REC. CELL PARAMETERS   0.041645  0.033724  0.030726  83.843  70.611  67.344
 COORDINATES OF UNIT CELL A-AXIS   -14.253   -14.249    18.612
 COORDINATES OF UNIT CELL B-AXIS    14.442    25.602    13.000
 COORDINATES OF UNIT CELL C-AXIS   -21.610    21.460   -16.235
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.014     0.015    -0.009
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.035
 LAB COORDINATES OF ROTATION AXIS -0.999993  0.001823  0.003120
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000910  0.005125  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.37   1244.05
 DETECTOR ORIGIN (PIXELS) AT                     1265.17   1237.32
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.25
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     463 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     32     57     67     57     27     48     63     54
 SIGMAB (degree)       0.033  0.034  0.034  0.033  0.031  0.032  0.032  0.033  0.033
 SIGMAR (degree)       0.027  0.025  0.027  0.028  0.026  0.023  0.026  0.028  0.027


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    463 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  2  7 10  7  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  1  4 11 15 12  8  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  3  6 10  8  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  1  2  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  1  3  6  7  6  3  1  0  0  0     0  0  0  2  7 12 13  9  4  2  0  0  0     0  0  0  2  5  9 10  7  3  1  0  0  0
     0  0  1  3 12 31 37 27 15  6  2  1  0     0  0  1  7 28 58 63 43 22  9  4  1  0     0  0  1  6 22 42 42 26 14  6  2  1  0
     0  0  1  5 21 47 58 45 25 11  4  1  0     0  0  1 11 46 96100 69 37 16  7  2  0     0  0  1  9 34 68 67 45 22  9  4  1  0
     0  0  0  2 13 28 36 28 14  7  2  1  0     0  0  1  7 28 58 61 44 20 10  3  1  0     0  0  1  6 22 41 41 27 13  6  2  0  0
     0  0  0  1  3  6  7  5  2  1  0  0  0     0  0  0  2  7 12 13  8  3  2  0  0  0     0  0  0  2  6 10  9  5  2  1  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  3  3  2  1  0  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 12 12  7  4  2  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0 -1
     0  0  1  3 10 20 19 12  5  2  1  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 12 11  8  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  3  3  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    613 REFLECTION INTENSITIES INTEGRATED


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
    41   0  0.999  5180356    0    112      70     0  0.03818  0.06962
    42   0  1.006  5181032    0    133      61     4  0.04269  0.06668
    43   0  1.002  5180102    0    128      71     5  0.04207  0.08331
    44   0  1.003  5179928    0    130      77     8  0.04572  0.07398
    45   0  1.005  5179990    0    132      81     6  0.04520  0.09473
    46   0  1.005  5180173    0    130      72     3  0.04471  0.08576
    47   0  1.007  5180471    0    131      69     3  0.04959  0.08304
    48   0  1.009  5179957    0    122      68     3  0.04474  0.08118
    49   0  1.010  5179863    0    122      77     2  0.04648  0.09520
    50   0  1.008  5180203    0    131      72     3  0.04153  0.08532
    51   0  1.008  5180071    0    133      73     5  0.04431  0.07373
    52   0  1.009  5179481    0    133      87     5  0.04110  0.09241
    53   0  1.010  5179946    0    128      74     2  0.04265  0.08000
    54   0  1.014  5181438    0    119      60     5  0.04078  0.07267
    55   0  1.014  5181450    0    105      58     4  0.03926  0.05167
    56   0  1.007  5180872    0    106      65     3  0.04689  0.07397
    57   0  1.011  5181129    0    121      65     6  0.04494  0.07113
    58   0  1.010  5180070    0    134      73     5  0.04093  0.07389
    59   0  1.011  5179352    0    140      78     2  0.04141  0.08552
    60   0  1.012  5179038    0    143      81     0  0.04321  0.09628

   473 OUT OF    473 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.51
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.432    32.136    34.498  88.672 108.461 111.877
 REC. CELL PARAMETERS   0.041644  0.033729  0.030739  83.816  70.609  67.344
 COORDINATES OF UNIT CELL A-AXIS   -14.252   -14.249    18.612
 COORDINATES OF UNIT CELL B-AXIS    14.444    25.592    13.004
 COORDINATES OF UNIT CELL C-AXIS   -21.607    21.445   -16.228
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.009     0.007    -0.005
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.035
 LAB COORDINATES OF ROTATION AXIS -0.999993  0.001823  0.003120
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000803  0.005044  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.32   1244.06
 DETECTOR ORIGIN (PIXELS) AT                     1265.26   1237.43
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.23
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     473 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     32     60     63     50     27     53     69     60
 SIGMAB (degree)       0.033  0.034  0.034  0.033  0.032  0.032  0.032  0.033  0.032
 SIGMAR (degree)       0.028  0.025  0.027  0.027  0.026  0.023  0.027  0.028  0.027


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    473 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  1  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  2  1  1  0  0  0  0     0  0  0  1  3  8 11  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  2  1  1  0  0  0     0  0  0  1  5 14 19 14  8  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  1  3  8 12  9  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  1  2  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  1  3  7  7  5  3  1  0  0  0     0  0  0  1  6 12 12  8  4  2  0  0  0     0  0  0  1  4  8  7  4  2  1  0  0  0
     0  0  0  4 13 32 40 31 16  6  2  1  0     0  0  1  8 27 57 62 45 22  8  3  1  0     0  0  1  6 18 36 36 25 12  4  2  0  0
     0  0  0  5 20 51 64 48 26 11  5  1  0     0  0  1 10 42 88100 68 34 14  6  2  0     0  0  1  8 30 56 57 37 18  8  3  1  0
     0  0  0  3 14 32 39 27 14  7  2  1  0     0  0  1  7 28 55 60 38 19  9  3  1  0     0  0  1  6 20 35 35 21 10  4  1  0  0
     0  0  0  1  3  7  7  5  3  1  0  0  0     0  0  0  3  7 13 12  8  4  1  0  0  0     0  0  0  2  5  8  7  4  2  1  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  9  8  5  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 14 13  9  4  1  1  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  9  8  5  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
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
    61   0  1.013  5178602    0    129      89     2  0.04383  0.09376
    62   0  1.012  5179689    0    120      73     3  0.04490  0.07923
    63   0  1.011  5181046    0    115      59     2  0.04888  0.09255
    64   0  1.012  5181071    0    123      57     4  0.04631  0.07428
    65   0  1.012  5180177    0    122      69     5  0.04401  0.07474
    66   0  1.016  5179775    0    123      76     8  0.04432  0.07279
    67   0  1.015  5180130    0    129      69     7  0.04292  0.06091
    68   0  1.015  5180098    0    138      70     4  0.04516  0.07165
    69   0  1.015  5179690    0    141      80     1  0.04330  0.08093
    70   0  1.016  5180391    0    127      70     2  0.03895  0.06383
    71   0  1.020  5181229    0    125      57     5  0.04193  0.05470
    72   0  1.017  5181327    0    122      59     3  0.04550  0.07737
    73   0  1.017  5180180    0    127      73     6  0.04414  0.06954
    74   0  1.019  5179813    0    137      77     0  0.04608  0.11949
    75   0  1.021  5180559    0    134      75     3  0.04279  0.08283
    76   0  1.018  5180486    0    134      71     2  0.04221  0.07822
    77   0  1.017  5180780    0    127      69     4  0.04035  0.05387
    78   0  1.016  5180346    0    121      72     0  0.03934  0.08265
    79   0  1.018  5180887    0    118      66     1  0.04109  0.08707
    80   0  1.017  5180695    0    129      66     1  0.04166  0.09098

   474 OUT OF    474 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.52
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.06
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.426    32.141    34.484  88.692 108.432 111.868
 REC. CELL PARAMETERS   0.041645  0.033722  0.030747  83.810  70.633  67.350
 COORDINATES OF UNIT CELL A-AXIS   -14.249   -14.242    18.611
 COORDINATES OF UNIT CELL B-AXIS    14.441    25.600    13.005
 COORDINATES OF UNIT CELL C-AXIS   -21.612    21.421   -16.222
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.029    -0.001    -0.012
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.035
 LAB COORDINATES OF ROTATION AXIS -0.999993  0.001823  0.003120
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000592  0.005249  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.27   1244.11
 DETECTOR ORIGIN (PIXELS) AT                     1265.49   1237.21
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.24
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     474 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   49     31     55     62     58     31     54     65     57
 SIGMAB (degree)       0.033  0.033  0.034  0.033  0.032  0.032  0.032  0.033  0.033
 SIGMAR (degree)       0.028  0.025  0.027  0.028  0.026  0.024  0.026  0.028  0.027


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    474 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  3  3  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  1  4 11 14 12  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  3  2  1  0  0  0     0  0  0  1  7 17 23 20 11  5  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  1  5 11 14 12  5  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  1  0  0  0  0  0  0  0  0  0  0  0     0  1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  4  8 11  7  4  1  0  0  0     0  0  0  1  6 13 14  8  4  2  0  0  0     0  0  0  1  3  7  7  4  2  1  0  0  0
     0  0  0  3 18 42 48 35 19  9  3  1  0     0  0  1  6 32 63 62 41 21 10  4  1  0     0  0  1  4 19 34 32 20 10  5  2  0  0
     0  0  1  6 30 64 76 60 31 14  5  2  1     0  0  2 12 51 96100 71 34 15  6  2  0     0  0  1  8 30 53 51 36 16  7  3  1  0
     0  0  1  5 19 38 48 33 16  8  2  1  0     0  0  1  9 33 58 63 40 19  9  3  1  0     0  0  1  6 20 32 32 20  9  4  1  0  0
     0  0  0  1  4  9 10  6  3  1  0  0  0     0  0  0  3  8 13 15  7  4  2  0  0  0     0  0  0  2  5  7  7  4  2  1  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  6  6  4  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  7 12 10  6  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  2  5  6  6  4  2  1  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  1  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    650 REFLECTION INTENSITIES INTEGRATED


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
    81   0  1.020  5180413    0    127      69     0  0.04178  0.09446
    82   0  1.019  5180177    0    126      71     5  0.04253  0.07598
    83   0  1.017  5181291    0    118      60     6  0.04220  0.06961
    84   0  1.020  5181740    0    112      55     2  0.04610  0.09020
    85   0  1.018  5180589    0    102      66     7  0.03897  0.05412
    86   0  1.019  5181169    0    113      62     2  0.04297  0.08901
    87   0  1.019  5181493    0    118      55     0  0.04221  0.10201
    88   0  1.017  5180366    0    130      70     1  0.04384  0.08496
    89   0  1.019  5179096    0    135      85     3  0.04193  0.07142
    90   0  1.018  5180289    0    140      68     0  0.04196  0.08278
    91   0  1.014  5179979    0    122      73     4  0.04070  0.06235
    92   0  1.020  5180927    0    123      65     2  0.04525  0.09026
    93   0  1.020  5181235    0    114      60     3  0.04097  0.07101
    94   0  1.023  5180696    0    120      68     1  0.04262  0.08729
    95   0  1.020  5180483    0    122      73     2  0.04203  0.09152
    96   0  1.021  5180337    0    124      72     7  0.04216  0.06323
    97   0  1.023  5180995    0    127      64     2  0.04598  0.11562
    98   0  1.022  5180600    0    127      72     3  0.04752  0.09595
    99   0  1.021  5179939    0    133      77     5  0.04289  0.07344
   100   0  1.018  5180049    0    139      70     7  0.04292  0.05815

   458 OUT OF    459 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.45
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.07
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.426    32.131    34.499  88.675 108.451 111.862
 REC. CELL PARAMETERS   0.041647  0.033730  0.030735  83.823  70.620  67.360
 COORDINATES OF UNIT CELL A-AXIS   -14.248   -14.242    18.610
 COORDINATES OF UNIT CELL B-AXIS    14.441    25.588    13.005
 COORDINATES OF UNIT CELL C-AXIS   -21.612    21.445   -16.225
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.010     0.004    -0.008
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.035
 LAB COORDINATES OF ROTATION AXIS -0.999993  0.001823  0.003120
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000916  0.005062  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.37   1244.12
 DETECTOR ORIGIN (PIXELS) AT                     1265.16   1237.47
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.23
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     458 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     29     55     63     51     30     60     67     52
 SIGMAB (degree)       0.032  0.033  0.033  0.033  0.032  0.032  0.031  0.032  0.032
 SIGMAR (degree)       0.029  0.024  0.027  0.028  0.026  0.024  0.026  0.028  0.027


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    458 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  1  0  0  0     0  0  0  1  3  9 12 10  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  1  0     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  1  6 16 20 16  9  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  4 10 12  9  5  2  1  0  1
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  3  7  8  7  3  1  0  0  0     0  0  0  1  6 11 12  8  3  2  1  0  0     0  0  0  1  4  7  7  5  2  1  0  0  0
     0  0  0  3 13 34 44 32 17  7  3  1  0     0  0  1  6 26 56 63 42 21  8  3  1  0     0  0  0  4 17 34 35 23 11  4  2  0  0
     0  0  0  5 25 60 72 52 28 12  5  2  0     0  0  1 11 46 94100 68 36 15  6  2  0     0  0  1  7 29 55 57 37 19  8  3  1  0
     0  0  0  3 15 35 41 31 15  6  2  1  0     0  0  1  7 28 57 62 41 20  8  3  1  0     0  0  1  5 18 34 36 21 10  4  1  0  0
     0  0  0  1  4  8  8  5  2  1  0  0  0     0  0  0  2  7 14 13  8  3  1  0  0  0     0  0  0  2  5  9  7  4  2  1  0  0  0
     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  8  8  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 13 13  8  4  2  1  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0
     0  0  0  2  5  8  9  5  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  2  2  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    626 REFLECTION INTENSITIES INTEGRATED


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
   101   0  1.022  5179822    0    140      77     4  0.04193  0.07449
   102   0  1.022  5179830    0    130      77     4  0.04182  0.06826
   103   0  1.022  5180818    0    117      64     7  0.04152  0.06183
   104   0  1.026  5181859    0    121      54     6  0.04597  0.06427
   105   0  1.021  5180862    0    125      68     3  0.04227  0.08517
   106   0  1.021  5181348    0    130      59     3  0.04026  0.06991
   107   0  1.022  5179933    0    128      77     1  0.04186  0.08428
   108   0  1.017  5180382    0    126      73     1  0.04080  0.09729
   109   0  1.023  5180571    0    128      72     2  0.04424  0.07709
   110   0  1.027  5180368    0    119      71     4  0.03949  0.05934
   111   0  1.026  5181215    0    115      65     4  0.04124  0.07101
   112   0  1.029  5181590    0    116      57     0  0.04192  0.10135
   113   0  1.029  5180477    0    119      70     5  0.04707  0.06544
   114   0  1.028  5180518    0    133      72     4  0.04066  0.06189
   115   0  1.025  5179916    0    148      78     3  0.04132  0.06718
   116   0  1.026  5179491    0    148      79     7  0.04040  0.05090
   117   0  1.026  5180260    0    150      75     3  0.04131  0.06705
   118   0  1.027  5180582    0    143      72     4  0.03843  0.05960
   119   0  1.028  5181356    0    136      66     3  0.04370  0.07781
   120   0  1.026  5180873    0    126      70     2  0.04068  0.08607

   484 OUT OF    485 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.53
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.03
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.437    32.138    34.494  88.682 108.455 111.883
 REC. CELL PARAMETERS   0.041639  0.033728  0.030742  83.805  70.610  67.335
 COORDINATES OF UNIT CELL A-AXIS   -14.248   -14.258    18.614
 COORDINATES OF UNIT CELL B-AXIS    14.443    25.596    13.006
 COORDINATES OF UNIT CELL C-AXIS   -21.615    21.435   -16.222
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.017    -0.012     0.001
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.035
 LAB COORDINATES OF ROTATION AXIS -0.999993  0.001823  0.003120
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000549  0.005038  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.32   1244.06
 DETECTOR ORIGIN (PIXELS) AT                     1265.59   1237.44
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.25
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     484 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     34     60     65     54     28     53     70     65
 SIGMAB (degree)       0.033  0.034  0.034  0.033  0.032  0.032  0.031  0.032  0.032
 SIGMAR (degree)       0.028  0.024  0.027  0.028  0.026  0.023  0.027  0.028  0.026


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    484 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  1  4  9 13 10  6  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  1  6 14 21 16  9  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  1  4 10 13  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0 -1  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0 -1  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  1  3  8 10  7  3  1  0  0  0     0  0  0  1  5 13 14  9  4  1  0  0  0     0  0  0  1  3  7  7  4  2  1  0  0  0
     0  0  0  4 15 35 43 32 17  8  3  1  0     0  0  1  7 27 59 62 42 21  9  3  1  0     0  0  0  4 17 33 33 20 10  5  2  0  0
     0  0  0  5 24 55 72 54 28 12  5  1  0     0  0  1 10 43 91100 70 36 16  5  2  0     0  0  1  7 27 52 51 35 17  8  3  1  0
     0  0  0  3 16 36 45 32 16  7  3  1  0     0  0  1  7 27 57 62 43 20  8  3  1  0     0  0  1  5 16 32 33 23  9  4  2  0  0
     0  0  0  1  4  8 12  7  3  1  0  0  0     0  0  0  2  7 14 16  8  4  2  0  0  0     0  0  0  1  5  7  8  4  2  1  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  1  1  2  1  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  7  4  2  1  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 11 11  7  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4  7  8  5  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
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
   121   0  1.027  5180966    0    114      67     0  0.04140  0.09638
   122   0  1.026  5181223    0    121      62     2  0.04225  0.06609
   123   0  1.029  5181715    0    117      59     3  0.03975  0.07807
   124   0  1.033  5182148    0    103      52     2  0.04171  0.07433
   125   0  1.034  5182610    0    110      48     0  0.04324  0.10220
   126   0  1.035  5181855    0    110      58     0  0.04339  0.09027
   127   0  1.034  5180233    0    111      70     4  0.04347  0.06374
   128   0  1.032  5180939    0    120      70     1  0.04124  0.09154
   129   0  1.031  5181437    0    121      60     4  0.04012  0.05944
   130   0  1.033  5181101    0    117      63     1  0.04268  0.08356
   131   0  1.033  5181552    0    116      64     0  0.04621  0.11667
   132   0  1.032  5181366    0    118      64     1  0.04325  0.08903
   133   0  1.032  5181130    0    128      60     4  0.04110  0.05307
   134   0  1.033  5180252    0    123      72     0  0.03819  0.06113
   135   0  1.034  5180256    0    127      74     0  0.04100  0.06838
   136   0  1.038  5181282    0    136      62     1  0.03853  0.06581
   137   0  1.033  5180465    0    122      69     3  0.04192  0.07422
   138   0  1.035  5180932    0    126      63     3  0.04234  0.07319
   139   0  1.038  5181153    0    112      62     5  0.04745  0.07688
   140   0  1.035  5181037    0    103      58     0  0.04464  0.09508

   455 OUT OF    457 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.54
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.05
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.429    32.141    34.497  88.655 108.457 111.882
 REC. CELL PARAMETERS   0.041648  0.033723  0.030738  83.836  70.620  67.345
 COORDINATES OF UNIT CELL A-AXIS   -14.246   -14.249    18.611
 COORDINATES OF UNIT CELL B-AXIS    14.441    25.603    13.001
 COORDINATES OF UNIT CELL C-AXIS   -21.611    21.443   -16.221
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.018    -0.001    -0.011
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.035
 LAB COORDINATES OF ROTATION AXIS -0.999993  0.001823  0.003120
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000734  0.005212  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.32   1244.12
 DETECTOR ORIGIN (PIXELS) AT                     1265.35   1237.28
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.23
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     455 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   45     31     56     69     53     29     51     63     49
 SIGMAB (degree)       0.033  0.033  0.034  0.033  0.032  0.032  0.032  0.032  0.033
 SIGMAR (degree)       0.027  0.025  0.027  0.028  0.026  0.024  0.026  0.027  0.027


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    455 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  1  2  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  1  0  0  0     0  0  0  0  2  6 10  8  4  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  1  0  0  0     0  0  0  1  5 10 15 13  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  2  7 10  8  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  3  7  8  6  3  1  0  0  0     0  0  0  1  6 13 15  9  5  2  1  0  0     0  0  0  1  4  9 10  6  3  1  0  0  0
     0  0  0  3 11 29 39 29 14  7  3  1  0     0  0  1  6 24 56 66 45 21  9  4  1  0     0  0  0  5 18 38 40 27 12  5  2  0  0
     0  0  0  5 21 45 59 46 25 11  4  1  0     0  0  1 10 42 89100 71 37 16  6  2  0     0  0  1  8 31 62 62 42 20  9  3  1  0
     0  0  0  2 11 28 36 27 14  6  2  1  0     0  0  0  6 25 53 63 42 21  8  3  1  0     0  0  0  5 19 36 39 24 12  5  2  0  0
     0  0  0  1  3  7  8  5  2  1  0  0  0     0  0  0  2  7 13 13  8  4  2  0  0  0     0  0  0  2  5  9  8  5  2  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  1  0  0  1  0 -1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0 -1  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0
     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5 10 10  7  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  2  8 15 15 10  5  2  1  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5 10  9  5  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  1  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    587 REFLECTION INTENSITIES INTEGRATED


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
   141   0  1.035  5181970    0    114      52     1  0.04139  0.08583
   142   0  1.032  5180023    0    130      75     1  0.04179  0.07699
   143   0  1.035  5180134    0    133      70     4  0.04133  0.05933
   144   0  1.036  5180321    0    122      72     4  0.03835  0.06641
   145   0  1.039  5182298    0    126      52     0  0.04112  0.09661
   146   0  1.039  5181995    0    122      53     2  0.03904  0.07802
   147   0  1.040  5180535    0    120      68     2  0.04172  0.08800
   148   0  1.043  5180193    0    141      76     5  0.04147  0.06078
   149   0  1.039  5180954    0    142      65     4  0.03906  0.05642
   150   0  1.044  5181096    0    130      65     1  0.04307  0.06973
   151   0  1.042  5181407    0    129      62     1  0.04207  0.07786
   152   0  1.043  5181121    0    131      63     5  0.04086  0.05791
   153   0  1.045  5180311    0    137      77     7  0.03906  0.05105
   154   0  1.046  5180432    0    134      76     6  0.04183  0.07879
   155   0  1.047  5180923    0    118      68     4  0.04190  0.08414
   156   0  1.046  5181420    0    121      59     3  0.04321  0.06213
   157   0  1.046  5180966    0    123      65     1  0.04083  0.08124
   158   0  1.044  5180288    0    118      72     1  0.04169  0.08004
   159   0  1.046  5180549    0    118      71     0  0.04402  0.10161
   160   0  1.050  5180430    0    122      69     2  0.04297  0.09644

   496 OUT OF    496 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.48
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.04
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.434    32.138    34.496  88.676 108.458 111.878
 REC. CELL PARAMETERS   0.041642  0.033728  0.030740  83.813  70.611  67.342
 COORDINATES OF UNIT CELL A-AXIS   -14.246   -14.254    18.614
 COORDINATES OF UNIT CELL B-AXIS    14.443    25.594    13.005
 COORDINATES OF UNIT CELL C-AXIS   -21.614    21.440   -16.221
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.013    -0.010    -0.004
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.035
 LAB COORDINATES OF ROTATION AXIS -0.999993  0.001823  0.003120
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000669  0.005183  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.28   1244.07
 DETECTOR ORIGIN (PIXELS) AT                     1265.41   1237.26
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.26
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     496 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   49     32     59     71     56     29     57     70     64
 SIGMAB (degree)       0.034  0.034  0.034  0.033  0.032  0.032  0.031  0.033  0.032
 SIGMAR (degree)       0.027  0.024  0.026  0.027  0.026  0.023  0.026  0.028  0.026


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    496 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  1  4  9 13  9  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0     0  0  0  1  7 15 20 16  9  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  1  5  9 13  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0


     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  4  8  9  6  3  1  0  0  0     0  0  0  1  6 12 13  9  4  2  0  0  0     0  0  0  1  4  7  7  4  2  1  0  0  0
     0  0  0  3 15 36 45 32 18  8  3  1  0     0  0  1  7 27 59 64 43 23 10  4  1  0     0  0  0  4 17 33 34 22 12  5  2  0  0
     0  0  0  6 27 56 69 53 28 13  4  2  0     0  0  1 11 47 92100 72 36 16  6  2  0     0  0  1  7 29 51 52 37 18  8  3  1  0
     0  0  0  4 17 35 45 31 15  7  3  1  0     0  0  1  7 31 57 64 41 20  9  3  1 -1     0  0  1  5 19 33 34 21 10  4  2  0 -1
     0  0  0  1  4  9  9  7  3  1  0  0  0     0  0  0  2  7 14 13  9  4  2  0  0  0     0  0  0  2  5  9  7  5  2  1  0  0  0
     0  0  0  0  1  1  1  1  0  0  0 -1  0     0  0  0  0  1  2  1  1  0  0  0 -1  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  7  5  2  1  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 11 11  7  3  1  1  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  7  7  4  2  1  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    638 REFLECTION INTENSITIES INTEGRATED


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
   161   0  1.051  5180536    0    123      69     5  0.04128  0.06221
   162   0  1.049  5180761    0    129      68     3  0.04038  0.06338
   163   0  1.048  5180137    0    116      75     6  0.04181  0.05734
   164   0  1.050  5181318    0    113      61     1  0.04121  0.08735
   165   0  1.052  5182378    0    110      52     2  0.04307  0.07946
   166   0  1.055  5181260    0    103      62     5  0.04200  0.06330
   167   0  1.055  5180997    0    104      67     3  0.04181  0.08775
   168   0  1.053  5181550    0    119      59     2  0.04355  0.08062
   169   0  1.050  5181112    0    124      60     3  0.03663  0.05111
   170   0  1.054  5179841    0    142      75     0  0.03987  0.06526
   171   0  1.055  5180578    0    140      69     6  0.03965  0.05966
   172   0  1.053  5181287    0    131      66     7  0.04124  0.05415
   173   0  1.054  5181239    0    130      69     4  0.03930  0.06739
   174   0  1.056  5180573    0    135      72     3  0.04069  0.06246
   175   0  1.056  5180856    0    145      66     3  0.04071  0.05657
   176   0  1.056  5180124    0    146      74     3  0.03979  0.05606
   177   0  1.058  5180647    0    148      63     4  0.03785  0.05671
   178   0  1.058  5180899    0    144      66     3  0.03971  0.06309
   179   0  1.056  5180758    0    141      64     0  0.03875  0.07035
   180   0  1.059  5180331    0    130      68     4  0.03920  0.06639

   479 OUT OF    481 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.46
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.07
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.433    32.140    34.495  88.660 108.468 111.887
 REC. CELL PARAMETERS   0.041648  0.033726  0.030742  83.823  70.605  67.337
 COORDINATES OF UNIT CELL A-AXIS   -14.243   -14.259    18.611
 COORDINATES OF UNIT CELL B-AXIS    14.440    25.600    13.005
 COORDINATES OF UNIT CELL C-AXIS   -21.613    21.443   -16.219
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.015    -0.013    -0.002
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.035
 LAB COORDINATES OF ROTATION AXIS -0.999993  0.001823  0.003120
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000649  0.005122  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.28   1244.11
 DETECTOR ORIGIN (PIXELS) AT                     1265.43   1237.38
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.23
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     479 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   46     35     62     70     55     29     57     66     54
 SIGMAB (degree)       0.033  0.033  0.034  0.033  0.032  0.032  0.032  0.033  0.033
 SIGMAR (degree)       0.027  0.024  0.027  0.027  0.025  0.024  0.026  0.028  0.027


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    479 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  1  2  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  3  8 10  8  4  2  1  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  1  0  0  0     0  0  0  1  5 12 15 13  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  2  7  8  7  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  1  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  3  7  7  6  2  1  0  0  0     0  0  0  2  5 13 13  8  4  2  0  0  0     0  0  0  1  4  9  8  5  2  1  0  0  0
     0  0  0  3 14 30 38 28 14  6  3  1  0     0  0  1  6 29 59 64 43 20  9  4  1  0     0  0  1  5 22 41 41 25 12  5  2  0  0
     0  0  0  5 20 48 59 45 23 10  4  1  0     0  0  1 11 42 91100 70 36 15  6  2  1     0  0  1  9 33 64 62 42 21  9  3  1  1
     0  0  0  3 12 30 34 26 13  6  2  1  0     0  0  1  7 26 56 60 42 21  9  3  1  0     0  0  1  6 21 39 39 26 12  5  2  0  0
     0  0  0  1  3  6  6  5  3  1  0  0  0     0  0  0  2  7 13 12  9  4  2  0  0  0     0  0  0  2  5  8  8  5  3  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0 -1  0  0  0  0  0  0  0  0  0     0  0  0 -1  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  6 10  9  5  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  9 16 14  9  4  2  1  0  0     0  0  0  0  1  2  1  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 10  9  6  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  1  0     0  0  0  0  0  0  0  0  0  0  0  1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    641 REFLECTION INTENSITIES INTEGRATED


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
   181   0  1.060  5180681    0    131      69     6  0.04014  0.06520
   182   0  1.057  5180719    0    129      68     4  0.04191  0.06483
   183   0  1.056  5180341    0    133      65     4  0.04023  0.04978
   184   0  1.058  5181173    0    127      60     3  0.04047  0.06404
   185   0  1.055  5180718    0    115      60     2  0.04144  0.07387
   186   0  1.057  5181068    0    122      58     6  0.04234  0.06269
   187   0  1.058  5180814    0    121      61     4  0.04060  0.05620
   188   0  1.056  5180269    0    118      66     2  0.04312  0.07138
   189   0  1.054  5180265    0    134      68     2  0.04133  0.07056
   190   0  1.053  5179231    0    135      83     6  0.03868  0.06171
   191   0  1.055  5179713    0    133      77     8  0.04270  0.07095
   192   0  1.060  5180471    0    133      73     5  0.04271  0.08444
   193   0  1.059  5180482    0    116      68     4  0.04233  0.08724
   194   0  1.058  5180758    0    117      60     2  0.04351  0.07609
   195   0  1.056  5180861    0    110      60     1  0.04203  0.07280
   196   0  1.056  5181493    0    113      56     1  0.04324  0.07684
   197   0  1.058  5180970    0    118      59     1  0.04036  0.07831
   198   0  1.059  5180647    0    114      65     5  0.04200  0.06102
   199   0  1.055  5180925    0    125      66     1  0.04077  0.10659
   200   0  1.057  5181092    0    123      63     1  0.04241  0.09940

   462 OUT OF    464 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.51
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.07
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.433    32.132    34.500  88.668 108.472 111.881
 REC. CELL PARAMETERS   0.041647  0.033734  0.030739  83.815  70.598  67.340
 COORDINATES OF UNIT CELL A-AXIS   -14.243   -14.261    18.610
 COORDINATES OF UNIT CELL B-AXIS    14.443    25.587    13.008
 COORDINATES OF UNIT CELL C-AXIS   -21.613    21.451   -16.217
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.000    -0.014     0.005
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.035
 LAB COORDINATES OF ROTATION AXIS -0.999993  0.001823  0.003120
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000767  0.005019  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.29   1244.06
 DETECTOR ORIGIN (PIXELS) AT                     1265.28   1237.46
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.25
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     462 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   48     25     53     71     52     27     51     67     57
 SIGMAB (degree)       0.033  0.033  0.034  0.033  0.031  0.032  0.031  0.032  0.032
 SIGMAR (degree)       0.027  0.025  0.026  0.027  0.026  0.023  0.026  0.028  0.027


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    462 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  1  3  8 11  9  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  1  5 13 17 14  8  3  2  1  1
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  3  8 11  8  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  0  1  0 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  1  0  0  0     0  0  0  0  0  0  0  0  0  1  0  0  0


     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  3  7  9  6  3  1  0  0  0     0  0  0  2  7 14 15  9  4  1  0  0  0     0  0  0  1  5  9  8  5  2  1  0  0  0
     0  0  0  3 12 32 39 29 14  6  2  1  0     0  0  0  6 27 60 63 40 20  8  3  1  0     0  0  0  4 19 38 37 22 10  5  2  0  0
     0  0  0  5 23 52 63 45 24 10  4  1  0     0  0  1 10 48 93100 68 33 14  6  2  1     0  0  1  7 32 58 58 38 18  7  3  1  0
     0  0  0  3 14 32 38 26 13  5  2  0  1     0  0  1  7 30 59 61 39 18  7  3  1  1     0  0  1  5 20 37 35 21  9  4  2  0  0
     0  0  0  1  3  6  6  4  2  1  0  0  0     0  0  0  2  6 11 11  7  3  1  0  0  0     0  0  0  2  4  7  6  4  2  1  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  1  0  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  1  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  9  9  5  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  9 15 14  9  4  2  1  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  9  8  4  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    641 REFLECTION INTENSITIES INTEGRATED


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
   201   0  1.057  5180238    0    121      72     8  0.04121  0.05804
   202   0  1.055  5180710    0    119      69     3  0.04253  0.09681
   203   0  1.052  5180909    0    123      68     1  0.04545  0.09623
   204   0  1.050  5180474    0    118      66     3  0.04150  0.06225
   205   0  1.051  5180254    0    119      74     5  0.03996  0.05840
   206   0  1.053  5180295    0    124      69     1  0.04339  0.08246
   207   0  1.054  5181338    0    118      58     0  0.04261  0.09476
   208   0  1.053  5180410    0    122      69     4  0.04287  0.07472
   209   0  1.055  5180472    0    122      66     0  0.04242  0.09230
   210   0  1.054  5181068    0    116      60     4  0.04113  0.07982
   211   0  1.054  5180821    0    107      63     1  0.04230  0.09667
   212   0  1.055  5181168    0    119      59     1  0.04450  0.08495
   213   0  1.052  5179994    0    128      75     4  0.04222  0.07349
   214   0  1.057  5180005    0    124      68     2  0.03821  0.06387
   215   0  1.057  5180384    0    129      71     1  0.04486  0.09317
   216   0  1.053  5180737    0    124      72     1  0.04492  0.09508
   217   0  1.055  5180974    0    123      63     1  0.04347  0.08650
   218   0  1.054  5180695    0    132      68     2  0.04102  0.07311
   219   0  1.053  5181320    0    128      58     2  0.04426  0.06525
   220   0  1.055  5180998    0    127      59     2  0.04079  0.06195

   461 OUT OF    463 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.50
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.05
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.426    32.130    34.500  88.677 108.456 111.860
 REC. CELL PARAMETERS   0.041647  0.033731  0.030736  83.819  70.614  67.361
 COORDINATES OF UNIT CELL A-AXIS   -14.244   -14.245    18.612
 COORDINATES OF UNIT CELL B-AXIS    14.444    25.584    13.007
 COORDINATES OF UNIT CELL C-AXIS   -21.613    21.449   -16.219
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.001    -0.004    -0.008
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.035
 LAB COORDINATES OF ROTATION AXIS -0.999993  0.001823  0.003120
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000947  0.005261  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.36   1244.08
 DETECTOR ORIGIN (PIXELS) AT                     1265.11   1237.16
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.25
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     461 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   45     28     58     69     54     29     57     63     51
 SIGMAB (degree)       0.033  0.033  0.034  0.033  0.031  0.031  0.031  0.032  0.032
 SIGMAR (degree)       0.028  0.026  0.027  0.027  0.026  0.022  0.026  0.028  0.027


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    461 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  1  4  9 13 11  6  2  1  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  2  3  3  2  1  0  0  0     0  0  0  1  6 15 20 17  9  4  2  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  1  4 10 12  8  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  7  9  6  3  1  0  0  0     0  0  0  1  6 12 14  8  4  2  0  0  0     0  0  0  1  4  7  7  4  2  1  0  0  0
     0  0  0  3 15 36 44 33 18  7  3  1  0     0  0  1  6 27 59 62 43 21  9  3  1  0     0  0  0  4 17 33 33 23 11  4  1  0  0
     0  0  1  5 25 60 71 54 28 12  5  2  0     0  0  1 10 46 95100 69 34 14  6  2  0     0  0  1  7 29 54 52 35 17  7  3  1  0
     0  0  0  3 14 37 42 28 14  6  3  1  0     0  0  1  7 27 60 61 38 18  8  3  1  0     0  0  1  4 17 33 33 20  9  4  1  0  0
     0  0  0  1  4  6  7  6  3  1  0  0  0     0  0  0  2  7 11 12  9  4  2  0  0  0     0  0  0  1  5  7  7  5  2  1  0  0  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  1  1  0  0  0  0  0  0  0  0     0  0  0  1  1  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  8  7  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 12 12  7  3  1  1  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  8  8  5  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    629 REFLECTION INTENSITIES INTEGRATED


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
   221   0  1.054  5180674    0    141      67     3  0.04037  0.06197
   222   0  1.055  5180729    0    138      70     1  0.03814  0.07142
   223   0  1.051  5180161    0    127      71     3  0.03988  0.05328
   224   0  1.053  5181483    0    125      60     2  0.03943  0.06254
   225   0  1.050  5180732    0    124      64     2  0.04065  0.06645
   226   0  1.052  5180501    0    135      69     2  0.04125  0.07711
   227   0  1.054  5180522    0    129      75     6  0.04315  0.06570
   228   0  1.051  5181115    0    136      67     2  0.03895  0.07814
   229   0  1.053  5180844    0    134      71     3  0.04044  0.05948
   230   0  1.052  5180474    0    128      71     4  0.03904  0.05507
   231   0  1.050  5181190    0    129      63     4  0.04240  0.05607
   232   0  1.053  5181229    0    117      60     0  0.03993  0.06877
   233   0  1.056  5181161    0    118      65     2  0.04252  0.06649
   234   0  1.057  5180974    0    122      70     4  0.03791  0.06468
   235   0  1.059  5180683    0    118      68     3  0.03964  0.05459
   236   0  1.058  5181146    0    123      65     5  0.04377  0.06171
   237   0  1.058  5181607    0    112      61     2  0.04159  0.07322
   238   0  1.060  5181358    0    121      62     1  0.04043  0.07060
   239   0  1.057  5181748    0    114      57     3  0.04371  0.06719
   240   0  1.058  5181121    0    117      65     6  0.03869  0.05630

   471 OUT OF    471 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.53
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.428    32.124    34.491  88.672 108.472 111.885
 REC. CELL PARAMETERS   0.041657  0.033743  0.030747  83.809  70.596  67.335
 COORDINATES OF UNIT CELL A-AXIS   -14.240   -14.262    18.605
 COORDINATES OF UNIT CELL B-AXIS    14.439    25.580    13.006
 COORDINATES OF UNIT CELL C-AXIS   -21.609    21.445   -16.211
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.001    -0.019     0.010
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.035
 LAB COORDINATES OF ROTATION AXIS -0.999993  0.001823  0.003120
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000691  0.004948  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.29   1244.10
 DETECTOR ORIGIN (PIXELS) AT                     1265.38   1237.61
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.18
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     471 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   48     38     57     69     50     25     49     65     63
 SIGMAB (degree)       0.033  0.033  0.033  0.033  0.032  0.033  0.031  0.032  0.032
 SIGMAR (degree)       0.027  0.025  0.027  0.028  0.027  0.024  0.027  0.028  0.027


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    471 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  1  3  8 12  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0     0  0  0  1  5 13 17 14  8  4  1  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  1  4  8 11  8  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0 -1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  3  6  7  5  2  1  0  0  0     0  0  0  2  6 12 13  8  3  1  0  0  0     0  0  0  1  5  7  8  5  2  1  0  0  0
     0  0  0  3 13 31 40 29 15  7  2  1  0     0  0  1  6 27 57 65 43 23 10  4  1  0     0  0  1  5 21 41 39 24 13  5  2  0  0
     0  0  1  4 21 49 60 48 25 11  4  1  1     0  0  1 10 45 93100 72 36 16  6  2  0     0  0  1  8 34 64 62 42 19  9  3  1  0
     0  0  0  4 14 31 37 27 13  6  2  1  0     0  0  1  7 27 58 62 43 19  8  3  1  0     0  0  1  6 21 38 37 25 11  5  2  0  0
     0  0  0  2  4  8  7  5  3  1  0  0  0     0  0  0  3  7 14 13  8  4  2  0  0  0     0  0  0  2  5  9  8  4  2  1  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0 -1  0  1  0  0  0  0  0  0  0  0     0  0 -1  0  1  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  1  0     1  0  0  0  0  1  0  0  0  0  0  1  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  1  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  6 10  9  6  3  1  0  0  0     0  0  0  1  0  1  1  0  0  0  0  0  0     0  0  0  1  0  0  0  0  0  0  0  0  0
     0  0  0  2  9 16 16 10  4  2  1  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0
     0  0  0  2  6 10  9  6  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    631 REFLECTION INTENSITIES INTEGRATED


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
   241   0  1.059  5181178    0    133      67     4  0.04327  0.09191
   242   0  1.058  5180148    0    128      77     1  0.04412  0.10188
   243   0  1.059  5180329    0    138      72     2  0.04147  0.08344
   244   0  1.057  5179574    0    137      83     5  0.04185  0.07338
   245   0  1.057  5179592    0    138      77     6  0.03814  0.05638
   246   0  1.058  5179890    0    141      76     1  0.04023  0.09109
   247   0  1.063  5179926    0    126      75     2  0.04113  0.07428
   248   0  1.059  5180848    0    124      61     4  0.04087  0.05723
   249   0  1.057  5180171    0    124      72     1  0.04314  0.08470
   250   0  1.059  5180091    0    118      76     4  0.04266  0.06678
   251   0  1.061  5180506    0    121      70     4  0.04322  0.07416
   252   0  1.062  5181198    0    123      63     3  0.04403  0.08664
   253   0  1.059  5179977    0    124      76     3  0.04362  0.08257
   254   0  1.059  5180466    0    126      70     1  0.04264  0.08325
   255   0  1.059  5180740    0    125      72     2  0.04454  0.10567
   256   0  1.057  5180440    0    120      70     1  0.04353  0.10723
   257   0  1.059  5180571    0    117      75     1  0.04301  0.12013
   258   0  1.059  5181379    0    124      61     2  0.04297  0.07121
   259   0  1.062  5180836    0    118      63     3  0.03711  0.05324
   260   0  1.062  5180675    0    122      68     2  0.04109  0.06801

   484 OUT OF    485 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.55
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.08
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.429    32.122    34.496  88.686 108.464 111.867
 REC. CELL PARAMETERS   0.041649  0.033742  0.030742  83.803  70.600  67.349
 COORDINATES OF UNIT CELL A-AXIS   -14.243   -14.253    18.610
 COORDINATES OF UNIT CELL B-AXIS    14.443    25.574    13.008
 COORDINATES OF UNIT CELL C-AXIS   -21.611    21.446   -16.217
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.004    -0.012     0.004
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.035
 LAB COORDINATES OF ROTATION AXIS -0.999993  0.001823  0.003120
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000771  0.005093  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.35   1244.05
 DETECTOR ORIGIN (PIXELS) AT                     1265.34   1237.36
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.23
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     484 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   52     31     59     69     54     29     59     66     53
 SIGMAB (degree)       0.033  0.033  0.034  0.033  0.032  0.032  0.031  0.032  0.032
 SIGMAR (degree)       0.027  0.025  0.027  0.028  0.027  0.025  0.025  0.028  0.027


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    484 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  3  3  1  1  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  1  0  0  0     0  0  0  1  5 10 14 11  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  3  2  1  0  0  0     0  0  0  1  8 16 22 17  9  4  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  1  0  0  0     0  0  0  1  5  9 12  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  1  1  3  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  1  3  8  9  7  3  1  0  0  0     0  0  0  1  5 12 13  8  3  1  0  0  0     0  0  0  1  3  7  7  4  2  1  0  0  0
     0  0  0  3 17 38 46 34 16  7  3  1  0     0  0  0  6 27 58 63 43 20  8  3  1  0     0  0  0  4 16 32 33 22 10  4  1  0  0
     0  0  1  5 28 60 74 54 27 12  5  2  0     0  0  1  9 47 93100 68 33 14  5  2  0     0  0  1  6 29 53 52 35 16  7  2  1  0
     0  0  0  4 18 34 42 30 16  7  3  1  0     0  0  1  7 29 54 58 39 18  8  3  1  0     0  0  1  5 18 31 31 20  9  4  1  0  0
     0  0  0  2  5  9 11  6  3  1  0  0  0     0  0  0  3  8 14 14  8  4  2  0  0  0     0  0  0  2  5  9  8  5  2  1  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  1  1  2  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  7  6  4  2  1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 11 10  6  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  1
     0  0  0  2  4  7  6  4  2  1  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    636 REFLECTION INTENSITIES INTEGRATED


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
   261   0  1.064  5180433    0    138      70     1  0.04120  0.08344
   262   0  1.063  5179428    0    124      80     8  0.04041  0.06079
   263   0  1.063  5179721    0    136      73     3  0.04018  0.06473
   264   0  1.062  5180163    0    130      73     3  0.04096  0.06730
   265   0  1.063  5180650    0    120      66     3  0.04113  0.06474
   266   0  1.066  5180188    0    124      73     2  0.04519  0.08284
   267   0  1.067  5180363    0    134      64     0  0.04562  0.10505
   268   0  1.065  5179968    0    127      75     0  0.04792  0.10838
   269   0  1.068  5179827    0    125      77     2  0.04266  0.08848
   270   0  1.067  5180084    0    132      72     1  0.04203  0.08318
   271   0  1.064  5180431    0    128      69     0  0.04552  0.10357
   272   0  1.065  5179474    0    118      77     3  0.04405  0.07686
   273   0  1.063  5180514    0    117      69     7  0.04478  0.06844
   274   0  1.062  5180862    0    111      65     3  0.04204  0.07980
   275   0  1.064  5180524    0    127      71     1  0.04457  0.08146
   276   0  1.063  5179742    0    137      79     4  0.04242  0.07940
   277   0  1.063  5179328    0    134      86     5  0.04175  0.06638
   278   0  1.065  5179510    0    132      81     2  0.04310  0.08749
   279   0  1.065  5180170    0    123      78     3  0.04356  0.10391
   280   0  1.062  5180700    0    124      68     2  0.04359  0.08708

   478 OUT OF    480 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.50
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.07
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.436    32.132    34.506  88.671 108.479 111.880
 REC. CELL PARAMETERS   0.041645  0.033734  0.030735  83.809  70.589  67.339
 COORDINATES OF UNIT CELL A-AXIS   -14.243   -14.265    18.611
 COORDINATES OF UNIT CELL B-AXIS    14.445    25.584    13.011
 COORDINATES OF UNIT CELL C-AXIS   -21.615    21.458   -16.218
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.007    -0.017     0.009
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.035
 LAB COORDINATES OF ROTATION AXIS -0.999993  0.001823  0.003120
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000783  0.005032  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.30   1244.04
 DETECTOR ORIGIN (PIXELS) AT                     1265.27   1237.43
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.29
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     478 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   50     34     53     66     52     25     52     68     65
 SIGMAB (degree)       0.033  0.033  0.034  0.033  0.032  0.032  0.031  0.032  0.032
 SIGMAR (degree)       0.027  0.025  0.027  0.028  0.027  0.025  0.028  0.028  0.027


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    478 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  1  4 10 14 11  6  2  1  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  3  3  1  1  0  0  0     0  0  0  2  7 16 22 17 10  4  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  1  4  9 12  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0 -1  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  3  7  9  6  3  1  0  0  0     0  0  0  1  5 11 12  8  4  1  0  0  0     0  0  0  1  3  6  7  4  2  1  0  0  0
     0  0  0  4 17 36 48 32 16  7  3  1  0     0  0  1  7 30 56 64 39 20  8  3  1  0     0  0  0  5 18 34 35 20 10  4  1  0  0
     0  0  0  6 27 59 73 51 28 12  5  2  0     0  0  1 11 49 92100 64 33 14  5  2  0     0  0  1  8 31 54 54 35 17  7  3  1  0
     0  0  0  4 16 35 41 28 15  6  3  1  0     0  0  1  7 28 55 59 36 18  8  3  1  0     0  0  0  5 18 33 34 20 10  4  2  0  0
     0  0  0  1  4  8  8  6  3  1  0  0  0     0  0  0  2  6 13 13  8  4  1  0  0  0     0  0  0  1  4  8  9  5  2  1  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  2  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  8  8  4  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  8 13 12  8  4  2  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  8  8  5  2  1  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    628 REFLECTION INTENSITIES INTEGRATED


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
   281   0  1.064  5180284    0    140      72     7  0.04026  0.06664
   282   0  1.068  5180285    0    129      75     7  0.04096  0.07250
   283   0  1.066  5181065    0    132      63     1  0.04225  0.10668
   284   0  1.060  5179854    0    133      77     8  0.04399  0.06286
   285   0  1.064  5179956    0    138      75     7  0.04161  0.05391
   286   0  1.067  5180440    0    146      69     2  0.04241  0.07516
   287   0  1.061  5179772    0    133      76     1  0.03912  0.07410
   288   0  1.065  5180107    0    145      74     3  0.03796  0.05742
   289   0  1.064  5181215    0    134      67     3  0.03685  0.06650
   290   0  1.064  5181012    0    131      68     2  0.03684  0.08003
   291   0  1.064  5181147    0    122      64     7  0.04199  0.06182
   292   0  1.068  5181719    0    116      56     6  0.04344  0.06732
   293   0  1.066  5181113    0    118      65     4  0.04390  0.08213
   294   0  1.069  5180801    0    118      71     4  0.03971  0.07009
   295   0  1.065  5181154    0    129      62     6  0.04139  0.05772
   296   0  1.067  5180221    0    131      76     4  0.04118  0.06482
   297   0  1.069  5180445    0    120      67     2  0.04161  0.08605
   298   0  1.070  5181368    0    122      60     0  0.04351  0.09979
   299   0  1.068  5181530    0    119      56     4  0.04531  0.06913
   300   0  1.071  5180651    0    120      66     7  0.03748  0.04697

   466 OUT OF    469 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.57
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.07
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.431    32.134    34.508  88.668 108.462 111.860
 REC. CELL PARAMETERS   0.041641  0.033726  0.030729  83.827  70.612  67.363
 COORDINATES OF UNIT CELL A-AXIS   -14.244   -14.253    18.612
 COORDINATES OF UNIT CELL B-AXIS    14.443    25.586    13.013
 COORDINATES OF UNIT CELL C-AXIS   -21.618    21.459   -16.217
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.005    -0.010     0.001
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.035
 LAB COORDINATES OF ROTATION AXIS -0.999993  0.001823  0.003120
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000852  0.005200  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.39   1244.14
 DETECTOR ORIGIN (PIXELS) AT                     1265.27   1237.30
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.30
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     466 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   44     34     56     63     52     29     53     69     58
 SIGMAB (degree)       0.033  0.034  0.034  0.033  0.032  0.032  0.031  0.032  0.033
 SIGMAR (degree)       0.027  0.025  0.027  0.028  0.027  0.023  0.026  0.028  0.028


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    466 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  4  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  1  4 11 15 12  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  3  1  1  0  0  0     0  0  0  1  7 17 22 18 10  4  2  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  1  5 11 14 10  6  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  2  3  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  4  9 12  7  3  1  0  0  0     0  0  0  2  6 15 16  9  4  2  0  0  0     0  0  0  1  4  9  9  5  2  1  0  0  0
     0  0  0  3 16 40 50 36 18  8  3  1  0     0  0  1  6 29 59 67 44 21  9  3  1  0     0  0  1  4 18 34 36 24 10  5  2  0  0
     0  0  0  5 28 63 75 59 29 13  5  2  0     0  0  1 11 49 94100 71 34 15  5  2  0     0  0  1  8 30 55 57 37 18  8  3  1  0
     0  0  0  4 18 42 47 34 17  7  3  1  0     0  0  1  9 32 63 61 40 20  8  3  1  0     0  0  0  6 19 35 33 21 11  4  2  0  0
     0  0  0  1  6 11 11  8  3  2  0  0  0     0  0  0  2  9 16 14 10  4  2  0  0  0     0  0  0  1  6  9  8  5  2  1  0  0  0
     0  0  0  0  1  1  2  1  0  0  0 -1  0     0  0  0  1  1  2  2  1  0  0  0 -1  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  8  8  5  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 13 12  8  4  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  8  7  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0 -1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    633 REFLECTION INTENSITIES INTEGRATED


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
   301   0  1.071  5181197    0    118      65     6  0.04156  0.06824
   302   0  1.074  5181101    0    129      62     2  0.04057  0.07499
   303   0  1.071  5180230    0    115      67     7  0.03953  0.05038
   304   0  1.074  5181586    0    123      55     2  0.04430  0.08059
   305   0  1.074  5181830    0    119      57     0  0.04312  0.09121
   306   0  1.077  5181129    0    120      60     2  0.04044  0.06792
   307   0  1.076  5181158    0    131      68     2  0.04294  0.07888
   308   0  1.077  5180540    0    131      72     5  0.04185  0.07025
   309   0  1.078  5180988    0    131      68     5  0.03838  0.06255
   310   0  1.075  5180709    0    133      68     4  0.03929  0.06117
   311   0  1.074  5180678    0    127      71     0  0.03981  0.09283
   312   0  1.076  5180080    0    125      78     0  0.03840  0.08523
   313   0  1.076  5180562    0    122      72     1  0.04042  0.07714
   314   0  1.075  5181890    0    125      57     0  0.04193  0.07645
   315   0  1.075  5180810    0    117      72     1  0.04052  0.07466
   316   0  1.072  5179809    0    118      80     1  0.04029  0.07978
   317   0  1.070  5181059    0    132      68     2  0.04544  0.08266
   318   0  1.071  5180376    0    125      71     2  0.04043  0.06789
   319   0  1.075  5180814    0    121      72     2  0.03963  0.08246
   320   0  1.075  5181408    0    118      65     1  0.03952  0.08849

   474 OUT OF    474 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.51
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.07
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.438    32.137    34.507  88.659 108.484 111.890
 REC. CELL PARAMETERS   0.041646  0.033731  0.030735  83.816  70.587  67.332
 COORDINATES OF UNIT CELL A-AXIS   -14.245   -14.272    18.607
 COORDINATES OF UNIT CELL B-AXIS    14.445    25.587    13.017
 COORDINATES OF UNIT CELL C-AXIS   -21.613    21.464   -16.214
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.013    -0.018     0.018
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.035
 LAB COORDINATES OF ROTATION AXIS -0.999993  0.001823  0.003120
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000783  0.005079  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.28   1244.03
 DETECTOR ORIGIN (PIXELS) AT                     1265.25   1237.35
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.31
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     474 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   50     26     53     67     56     31     55     70     56
 SIGMAB (degree)       0.033  0.033  0.034  0.033  0.032  0.032  0.031  0.032  0.032
 SIGMAR (degree)       0.028  0.026  0.027  0.028  0.027  0.025  0.026  0.028  0.027


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    474 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  2  2  1  1  0  0  0     0  0  0  1  4  9 11 10  5  2  1  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  2  3  3  1  1  0  0  0     0  0  0  1  6 15 19 16  8  4  1  1  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  1  4 10 13 10  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  1  3  7  9  6  3  1  0  0  0     0  0  0  1  5 11 14  8  4  1  0  0  0     0  0  0  1  3  7  8  5  2  1  0  0  0
     0  0  0  4 15 34 43 32 16  6  2  1  0     0  0  0  7 27 56 63 42 20  8  3  1  0     0  0  0  5 18 34 37 23 11  4  2  0  0
     0  0  0  5 25 56 69 53 26 12  4  1  0     0  0  1 10 44 91100 69 33 15  5  2  1     0  0  1  8 31 56 56 37 18  7  3  1  0
     0  0  0  3 16 37 43 30 15  7  2  1  0     0  0  1  7 28 59 61 40 20  8  3  1  0     0  0  1  5 19 36 34 22 10  4  1  0  0
     0  0  0  1  4  7  8  6  3  1  0  0  0     0  0  0  2  7 13 12  8  4  2  0  0  0     0  0  0  1  5  8  7  5  2  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  2  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  8  5  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 13 12  8  4  2  1  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  9  8  5  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  1  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0 -1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    620 REFLECTION INTENSITIES INTEGRATED


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
   321   0  1.078  5181034    0    117      64     2  0.04330  0.08471
   322   0  1.078  5180812    0    125      67     2  0.03819  0.08544
   323   0  1.076  5181188    0    130      63     4  0.04269  0.06931
   324   0  1.079  5180446    0    127      73     4  0.04347  0.07834
   325   0  1.079  5180258    0    121      78     2  0.04212  0.09277
   326   0  1.077  5180803    0    117      65     2  0.04068  0.08447
   327   0  1.077  5181090    0    106      66     0  0.04247  0.11244
   328   0  1.078  5181428    0    113      60     2  0.04123  0.09753
   329   0  1.076  5181589    0    114      57     1  0.04013  0.08201
   330   0  1.078  5180485    0    112      70     0  0.03971  0.09499
   331   0  1.079  5180719    0    123      68     3  0.04219  0.08314
   332   0  1.081  5180437    0    117      75     3  0.04486  0.09612
   333   0  1.077  5180460    0    117      68     0  0.04378  0.09799
   334   0  1.079  5180767    0    127      66     7  0.04600  0.06849
   335   0  1.079  5180025    0    126      80     5  0.04734  0.08345
   336   0  1.078  5180636    0    119      69     7  0.04430  0.07240
   337   0  1.080  5180968    0    117      70     6  0.04481  0.07809
   338   0  1.084  5180374    0    118      70     4  0.04482  0.09295
   339   0  1.082  5180528    0    123      71     8  0.04288  0.06094
   340   0  1.082  5179946    0    133      76     1  0.04270  0.08130

   469 OUT OF    469 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.51
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.425    32.122    34.506  88.654 108.473 111.859
 REC. CELL PARAMETERS   0.041650  0.033738  0.030732  83.838  70.606  67.368
 COORDINATES OF UNIT CELL A-AXIS   -14.241   -14.250    18.608
 COORDINATES OF UNIT CELL B-AXIS    14.439    25.575    13.009
 COORDINATES OF UNIT CELL C-AXIS   -21.611    21.466   -16.213
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.013    -0.007     0.003
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.035
 LAB COORDINATES OF ROTATION AXIS -0.999993  0.001823  0.003120
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.001043  0.005098  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.39   1244.14
 DETECTOR ORIGIN (PIXELS) AT                     1265.02   1237.45
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.24
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     469 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   48     32     60     68     54     25     50     64     56
 SIGMAB (degree)       0.033  0.034  0.034  0.033  0.032  0.032  0.031  0.032  0.032
 SIGMAR (degree)       0.029  0.026  0.027  0.028  0.026  0.024  0.026  0.028  0.028


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    469 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  1  3  8 12 10  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  3  1  1  0  0  0     0  0  0  1  6 14 19 16  9  4  1  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  1  0  0  0     0  0  0  1  3  9 11  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  3  6  8  6  3  1  0  0  0     0  0  0  1  6 11 12  8  4  1  1  0  0     0  0  0  1  5  8  9  5  3  1  0  0  0
     0  0  0  3 16 34 43 32 15  7  3  1  0     0  0  1  6 29 57 64 45 21  9  3  1  0     0  0  1  5 22 40 41 27 13  5  2  0  0
     0  0  0  5 25 54 68 50 26 11  4  1  0     0  0  1 10 46 89100 69 34 15  6  2  0     0  0  1  8 32 58 62 41 20  9  3  1  0
     0  0  0  4 14 36 42 29 15  6  2  1  0     0  0  1  6 25 57 61 41 19  9  3  1  0     0  0  0  5 18 37 37 25 11  5  2  0  0
     0  0  0  1  4  9 11  7  3  1  0  0  0     0  0  0  2  7 13 16  9  4  2  0  0  0     0  0  0  1  5  9 10  6  2  1  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  3  3  2  1  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 12 11  7  3  1  0  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2 10 17 16 11  5  2  1  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5 10 10  6  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  2  3  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    621 REFLECTION INTENSITIES INTEGRATED


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
   341   0  1.085  5180016    0    134      75     2  0.03975  0.06654
   342   0  1.073  5180296    0    132      73     2  0.04310  0.07757
   343   0  1.077  5180965    0    122      66     2  0.03923  0.07512
   344   0  1.079  5180623    0    123      65     4  0.03591  0.04492
   345   0  1.077  5180787    0    117      67     3  0.03479  0.04601
   346   0  1.081  5181373    0    130      64     2  0.03729  0.06462
   347   0  1.081  5180822    0    129      69     3  0.04071  0.07315
   348   0  1.077  5181465    0    133      59     2  0.04091  0.07479
   349   0  1.078  5180981    0    126      66     3  0.04006  0.07224
   350   0  1.079  5180064    0    113      74     5  0.04279  0.06607
   351   0  1.079  5180821    0    115      67     2  0.04405  0.10474
   352   0  1.080  5180826    0    118      65     2  0.04369  0.09333
   353   0  1.083  5180769    0    111      67     1  0.04540  0.11904
   354   0  1.084  5181061    0    117      62     2  0.04486  0.09636
   355   0  1.083  5180610    0    121      64     1  0.03892  0.07514
   356   0  1.081  5179784    0    132      77     0  0.04058  0.08117
   357   0  1.080  5179807    0    129      74     2  0.04214  0.08546
   358   0  1.080  5180180    0    131      70     7  0.04082  0.06700
   359   0  1.082  5180257    0    127      73     2  0.03920  0.09778
   360   0  1.080  5180937    0    128      69     1  0.04089  0.10416

   475 OUT OF    476 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.52
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.05
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.439    32.146    34.505  88.654 108.483 111.893
 REC. CELL PARAMETERS   0.041644  0.033722  0.030736  83.820  70.590  67.330
 COORDINATES OF UNIT CELL A-AXIS   -14.246   -14.275    18.606
 COORDINATES OF UNIT CELL B-AXIS    14.441    25.599    13.020
 COORDINATES OF UNIT CELL C-AXIS   -21.613    21.458   -16.218
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.001    -0.016     0.021
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.035
 LAB COORDINATES OF ROTATION AXIS -0.999993  0.001823  0.003120
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000617  0.005141  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.25   1244.09
 DETECTOR ORIGIN (PIXELS) AT                     1265.44   1237.33
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.33
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     475 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     32     55     64     58     30     54     65     60
 SIGMAB (degree)       0.033  0.034  0.034  0.033  0.032  0.032  0.031  0.032  0.032
 SIGMAR (degree)       0.027  0.026  0.026  0.028  0.026  0.024  0.026  0.028  0.027


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    475 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  1  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  3  7 10  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0     0  0  0  1  6 13 17 13  7  3  1  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  1  4  8 11  8  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  3  9  9  7  3  1  0  0  0     0  0  0  2  6 14 14 10  4  2  0  0  0     0  0  0  1  5  8  8  6  2  1  0  0  0
     0  0  0  3 14 31 40 31 16  7  3  1  0     0  0  1  7 28 56 63 42 21  9  4  1  0     0  0  1  5 20 36 37 23 12  5  2  0  0
     0  0  0  4 22 53 64 49 24 11  4  1  0     0  0  1 10 45 94100 69 33 14  6  2  0     0  0  1  8 32 62 61 39 18  8  3  1  0
     0  0  0  3 15 33 41 28 14  6  2  0  0     0  0  1  7 29 57 64 41 19  9  3  1  0     0  0  1  5 20 37 39 23 11  5  1  0  0
     0  0  0  1  4  7  8  6  3  1  0  0  0     0  0  0  2  7 13 13  8  4  2  0  0  0     0  0  0  1  5  8  8  5  3  1  0  0  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6  9  9  5  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  9 16 15  9  4  2  1  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  9  9  5  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    638 REFLECTION INTENSITIES INTEGRATED


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
   361   0  1.083  5180274    0    122      75     8  0.03807  0.05458
   362   0  1.083  5180209    0    132      76     7  0.04136  0.05850
   363   0  1.077  5180012    0    127      77     8  0.04114  0.06573
   364   0  1.082  5180607    0    120      71     5  0.04289  0.07517
   365   0  1.083  5180836    0    126      65     4  0.04127  0.07762
   366   0  1.083  5180539    0    122      74     6  0.04046  0.07116
   367   0  1.081  5179669    0    142      80     4  0.04038  0.07373
   368   0  1.080  5179456    0    148      82     2  0.03732  0.07230
   369   0  1.080  5179668    0    148      80     3  0.03909  0.06923
   370   0  1.083  5180271    0    130      75     2  0.04212  0.07751
   371   0  1.080  5180748    0    127      71     3  0.04281  0.08215
   372   0  1.083  5180754    0    122      72     2  0.04173  0.08404
   373   0  1.087  5180941    0    116      67     1  0.04038  0.10412
   374   0  1.089  5181682    0    118      59     0  0.04220  0.11901
   375   0  1.091  5181621    0    115      59     6  0.04068  0.07265
   376   0  1.086  5181277    0    121      64     4  0.04395  0.07654
   377   0  1.089  5180427    0    126      71     1  0.04307  0.06482
   378   0  1.090  5180396    0    130      71     0  0.04099  0.07253
   379   0  1.094  5179928    0    130      75     2  0.03656  0.05689
   380   0  1.089  5181104    0    129      68     2  0.04330  0.09116

   481 OUT OF    481 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.54
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.428    32.121    34.493  88.692 108.469 111.873
 REC. CELL PARAMETERS   0.041654  0.033745  0.030746  83.792  70.592  67.340
 COORDINATES OF UNIT CELL A-AXIS   -14.244   -14.258    18.605
 COORDINATES OF UNIT CELL B-AXIS    14.440    25.573    13.011
 COORDINATES OF UNIT CELL C-AXIS   -21.608    21.442   -16.220
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.001    -0.010     0.014
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.035
 LAB COORDINATES OF ROTATION AXIS -0.999993  0.001823  0.003120
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000684  0.005270  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.27   1244.09
 DETECTOR ORIGIN (PIXELS) AT                     1265.37   1237.17
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.22
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     481 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   48     31     58     67     56     33     56     68     53
 SIGMAB (degree)       0.033  0.034  0.034  0.033  0.032  0.032  0.031  0.032  0.032
 SIGMAR (degree)       0.027  0.026  0.028  0.028  0.027  0.024  0.026  0.028  0.027


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    481 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  1  5 12 16 12  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  1  2  3  3  2  1  0  0  0     0  0  0  1  8 18 25 19 11  5  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  1  6 12 15 12  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  2  3  4  3  1  1  0  0  0
    -1  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0 -1  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  8 10  7  4  1  1  0  0     0  0  0  1  6 12 12  8  4  2  1  0  0     0  0  0  1  3  7  6  4  2  1  0  0  0
     0  0  0  3 17 41 51 37 18  8  3  1  0     0  0  0  5 27 59 67 44 21  9  3  1  0     0  0  0  4 16 31 34 21 10  4  1  0  0
     0  0  1  5 27 64 78 58 30 14  5  2  1     0  0  1 10 44 93100 68 34 15  6  2  0     0  0  1  7 26 50 49 31 16  7  2  1  0
     0  0  0  4 19 41 48 37 19  8  3  1  0     0  0  1  8 30 60 63 44 22  9  3  1  0     0  0  1  5 18 33 31 20 10  4  1  0  0
     0  0  0  1  5 11 12  8  4  2  0  0  0     0  0  0  2  8 16 16 10  4  2  0  0  0     0  0  0  1  4  8  8  5  2  1  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0 -1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  8  5  2  1  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  2  7 11 11  6  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  8  7  4  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    627 REFLECTION INTENSITIES INTEGRATED


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
   381   0  1.085  5180637    0    141      74     6  0.04013  0.06003
   382   0  1.086  5180075    0    141      75     2  0.03947  0.05922
   383   0  1.081  5180393    0    131      77     6  0.03660  0.05397
   384   0  1.079  5181126    0    131      66     1  0.04342  0.09508
   385   0  1.078  5179931    0    113      79     2  0.04335  0.09329
   386   0  1.077  5180604    0    123      68     3  0.03937  0.07454
   387   0  1.074  5180982    0    122      67     1  0.04265  0.09783
   388   0  1.082  5181362    0    115      61     1  0.04363  0.08690
   389   0  1.077  5181125    0    113      63     3  0.04261  0.07528
   390   0  1.080  5180637    0    118      70     4  0.04146  0.07000
   391   0  1.084  5180345    0    126      71     5  0.04193  0.06995
   392   0  1.080  5179807    0    127      78     5  0.03870  0.05430
   393   0  1.080  5180230    0    123      68     3  0.03778  0.06500
   394   0  1.079  5180962    0    118      65     3  0.03959  0.07938
   395   0  1.081  5181041    0    121      63     4  0.04254  0.06135
   396   0  1.079  5181374    0    125      58     5  0.03911  0.05695
   397   0  1.076  5180772    0    119      69     5  0.04065  0.07346
   398   0  1.077  5181169    0    125      62     0  0.04233  0.08224
   399   0  1.078  5180751    0    114      70     3  0.04012  0.07576
   400   0  1.076  5180778    0    129      62     4  0.03934  0.06717

   484 OUT OF    484 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.56
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.08
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.431    32.130    34.505  88.655 108.480 111.875
 REC. CELL PARAMETERS   0.041650  0.033734  0.030735  83.828  70.595  67.350
 COORDINATES OF UNIT CELL A-AXIS   -14.245   -14.262    18.605
 COORDINATES OF UNIT CELL B-AXIS    14.440    25.582    13.016
 COORDINATES OF UNIT CELL C-AXIS   -21.608    21.464   -16.216
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.011    -0.007     0.015
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.035
 LAB COORDINATES OF ROTATION AXIS -0.999993  0.001823  0.003120
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000814  0.005142  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.34   1244.10
 DETECTOR ORIGIN (PIXELS) AT                     1265.28   1237.35
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.27
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     484 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   48     27     55     72     61     30     56     67     56
 SIGMAB (degree)       0.033  0.033  0.034  0.033  0.032  0.032  0.031  0.032  0.032
 SIGMAR (degree)       0.028  0.028  0.028  0.028  0.026  0.024  0.026  0.027  0.028


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    484 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  1  4  9 12  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  1  8 15 20 15  9  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  1  4  9 12  9  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  1  4  7 10  5  3  1  0  0  0     0  0  0  2  7 12 15  7  4  1  0  0  0     0  0  0  1  4  8  9  4  2  1  0  0  0
     0  0  0  3 15 35 40 29 16  7  3  1  0     0  0  1  7 28 57 59 39 20  8  3  1  0     0  0  0  5 19 36 35 22 11  4  2  0  0
     0  0  0  5 27 55 69 49 26 12  4  1  0     0  0  1 10 49 92100 66 33 15  6  2  0     0  0  1  7 32 57 57 36 18  8  3  1  0
     0  0  0  4 15 30 42 28 13  6  2  1  0     0  0  1  7 28 50 61 39 17  8  3  1  0     0  0  1  5 18 32 34 20  9  4  1  0  0
     0  0  0  1  4  8  8  5  3  1  0  0  0     0  0  0  2  7 12 12  7  3  1  0  0  0     0  0  0  1  4  7  7  4  2  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  1  1  2  1  1  0  0  0 -1  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  1  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  1  0     0  0  0  0  0  0  0  0  0  0  0  1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  9  8  5  2  1  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  9 15 14  8  4  1  1  0  0     0  0  0  0  1  2  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  8  9  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    636 REFLECTION INTENSITIES INTEGRATED


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
   401   0  1.079  5180389    0    129      71     4  0.04358  0.07841
   402   0  1.077  5180208    0    121      74     0  0.04032  0.10405
   403   0  1.077  5180322    0    131      73     2  0.04067  0.09864
   404   0  1.078  5179918    0    121      76     2  0.03890  0.08474
   405   0  1.077  5180363    0    112      68     2  0.04174  0.09116
   406   0  1.076  5180218    0    113      70     1  0.04652  0.10441
   407   0  1.078  5181437    0    116      60     0  0.04340  0.12297
   408   0  1.081  5181655    0    122      56     2  0.03920  0.09072
   409   0  1.077  5180478    0    125      70     2  0.04208  0.08855
   410   0  1.077  5180225    0    127      70     1  0.03999  0.08208
   411   0  1.076  5180378    0    130      67     3  0.03931  0.06592
   412   0  1.083  5180165    0    123      69     2  0.04149  0.07025
   413   0  1.078  5180208    0    127      70     5  0.03892  0.06144
   414   0  1.076  5180280    0    132      67     4  0.04117  0.07261
   415   0  1.077  5179524    0    124      78     2  0.04353  0.08374
   416   0  1.074  5179301    0    129      81     7  0.04220  0.06127
   417   0  1.074  5179486    0    127      71     5  0.04343  0.06081
   418   0  1.079  5180836    0    120      61     3  0.04387  0.07673
   419   0  1.080  5180790    0    110      60     2  0.04508  0.08980
   420   0  1.078  5180872    0    103      58     0  0.04434  0.12679

   476 OUT OF    476 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.51
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.15
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.426    32.130    34.505  88.622 108.487 111.889
 REC. CELL PARAMETERS   0.041659  0.033735  0.030735  83.857  70.600  67.346
 COORDINATES OF UNIT CELL A-AXIS   -14.243   -14.266    18.597
 COORDINATES OF UNIT CELL B-AXIS    14.436    25.583    13.016
 COORDINATES OF UNIT CELL C-AXIS   -21.603    21.477   -16.207
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.021    -0.006     0.021
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.035
 LAB COORDINATES OF ROTATION AXIS -0.999993  0.001823  0.003120
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000961  0.005059  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.34   1244.02
 DETECTOR ORIGIN (PIXELS) AT                     1265.08   1237.38
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.21
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     476 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     30     61     68     55     29     56     66     56
 SIGMAB (degree)       0.032  0.033  0.033  0.033  0.032  0.032  0.031  0.032  0.032
 SIGMAR (degree)       0.027  0.027  0.028  0.028  0.026  0.024  0.026  0.028  0.028


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    476 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  2  6  9  6  4  1  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  4 10 14 11  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  3  6  8  6  3  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  2  5  7  5  3  1  0  0  0     0  0  0  1  4 12 12  8  4  1  0  0  0     0  0  0  1  4 10  9  6  3  1  0  0  0
     0  0  0  2 12 27 36 25 13  6  2  1  0     0  0  0  5 25 54 62 40 20  8  3  1  0     0  0  0  5 23 43 45 27 14  5  2  0  0
     0  0  0  4 22 47 58 43 23 10  4  1  0     0  0  1 10 44 90100 70 37 15  6  2  0     0  0  1  9 37 70 73 48 25 10  4  1  0
     0  0  0  3 13 29 35 24 12  5  2  0  0     0  0  1  6 27 55 61 40 20  8  3  1  0     0  0  1  5 21 44 45 27 13  5  2  0  0
     0  0  0  1  3  6  6  5  2  1  0  0  0     0  0  0  1  6 12 12  8  3  1  0  0  0     0  0  0  1  4 10  8  6  2  1  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  4  3  2  1  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  8 15 15  9  4  2  1  0  0     0  0  0  0  1  2  2  1  1  0  0  0 -1     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3 12 25 24 15  8  3  1  0  0     0  0  0  1  2  3  3  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 15 14  8  4  2  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  3  3  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    614 REFLECTION INTENSITIES INTEGRATED


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
   421   0  1.083  5180282    0    109      67     0  0.04765  0.14629
   422   0  1.077  5179665    0    117      75     1  0.04361  0.11641
   423   0  1.079  5179392    0    119      77     1  0.04553  0.09782
   424   0  1.075  5179624    0    127      76     2  0.04333  0.09400
   425   0  1.076  5178912    0    125      78     5  0.03989  0.07310
   426   0  1.079  5179536    0    124      76     2  0.03966  0.08142
   427   0  1.078  5180643    0    124      65     2  0.04614  0.09128
   428   0  1.076  5180325    0    116      63     4  0.04377  0.07550
   429   0  1.076  5180083    0    122      67     5  0.04179  0.06790
   430   0  1.076  5180558    0    122      63     3  0.03968  0.06091
   431   0  1.079  5180952    0    125      59     2  0.04012  0.06249
   432   0  1.079  5179583    0    123      73     3  0.04425  0.07475
   433   0  1.078  5180179    0    113      65     4  0.04449  0.07596
   434   0  1.079  5180592    0    123      67     2  0.04273  0.09223
   435   0  1.080  5179556    0    130      76     3  0.04523  0.08029
   436   0  1.080  5180001    0    136      72     1  0.04039  0.08681
   437   0  1.079  5179583    0    135      77     5  0.03968  0.06196
   438   0  1.078  5180160    0    127      72     8  0.04418  0.06823
   439   0  1.080  5180232    0    137      72     1  0.04395  0.10248
   440   0  1.077  5180035    0    133      77     4  0.04611  0.09457

   472 OUT OF    472 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.50
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.430    32.129    34.505  88.648 108.479 111.873
 REC. CELL PARAMETERS   0.041649  0.033734  0.030734  83.837  70.600  67.355
 COORDINATES OF UNIT CELL A-AXIS   -14.248   -14.262    18.601
 COORDINATES OF UNIT CELL B-AXIS    14.440    25.578    13.021
 COORDINATES OF UNIT CELL C-AXIS   -21.603    21.472   -16.212
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.023    -0.001     0.021
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.035
 LAB COORDINATES OF ROTATION AXIS -0.999993  0.001823  0.003120
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000917  0.005110  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.34   1244.13
 DETECTOR ORIGIN (PIXELS) AT                     1265.14   1237.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.26
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     472 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   49     29     55     62     54     31     57     70     56
 SIGMAB (degree)       0.033  0.033  0.034  0.034  0.032  0.032  0.031  0.032  0.032
 SIGMAR (degree)       0.027  0.026  0.027  0.028  0.026  0.023  0.025  0.027  0.027


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    472 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  3  8 10  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  1  0  0  0     0  0  0  1  6 14 18 13  7  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  4  8 10  8  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  3  7  9  7  3  1  0  0  0     0  0  0  2  6 13 14  9  4  1  0  0  0     0  0  0  1  4  9  9  6  3  1  0  0  0
     0  0  0  3 16 33 39 32 15  7  2  1  0     0  0  1  6 30 58 61 45 21  9  3  1  0     0  0  1  5 21 41 39 28 13  5  2  1  0
     0  0  0  4 25 55 66 49 25 12  4  1  0     0  0  1 10 47 94100 71 35 16  6  2  0     0  0  1  7 33 62 61 43 20  9  3  1  0
     0  0  0  3 15 33 40 29 14  6  2  1  0     0  0  1  6 29 57 64 43 20  9  3  1  0     0  0  1  5 21 37 40 26 12  5  2  1  0
     0  0  0  1  3  8  8  6  3  1  0  0  0     0  0  0  2  7 14 13 10  5  2  0  0  0     0  0  0  1  5  9  8  6  3  1  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  3  2  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 11 10  7  3  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  9 16 15 10  5  2  1  0  0     0  0  0  1  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 11 10  6  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  3  2  2  1  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    621 REFLECTION INTENSITIES INTEGRATED


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
   441   0  1.075  5179401    0    131      82     8  0.04031  0.06306
   442   0  1.077  5179274    0    129      79     2  0.04085  0.08506
   443   0  1.079  5179697    0    128      77     4  0.03857  0.07073
   444   0  1.077  5179476    0    140      78     2  0.04190  0.07686
   445   0  1.076  5179908    0    133      77     3  0.04248  0.09220
   446   0  1.076  5180189    0    120      74     4  0.04028  0.08499
   447   0  1.079  5180485    0    114      66     4  0.04142  0.06789
   448   0  1.077  5180910    0    116      62     4  0.03828  0.10803
   449   0  1.076  5180954    0    121      69     3  0.04445  0.13181
   450   0  1.080  5179620    0    126      74     2  0.03800  0.10275
   451   0  1.080  5179343    0    129      85     2  0.04044  0.08875
   452   0  1.081  5179787    0    128      77     2  0.03966  0.08177
   453   0  1.073  5180306    0    124      71     1  0.04024  0.09112
   454   0  1.078  5180539    0    116      74     2  0.04167  0.10006
   455   0  1.077  5180730    0    110      69     1  0.04103  0.10611
   456   0  1.075  5180885    0    115      69     2  0.04267  0.11403
   457   0  1.076  5180875    0    114      69     6  0.04183  0.07779
   458   0  1.080  5181077    0    116      63     3  0.03939  0.08816
   459   0  1.081  5181176    0    124      66     3  0.04110  0.09464
   460   0  1.080  5180690    0    126      70     7  0.04262  0.07433

   470 OUT OF    472 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.51
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.14
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.436    32.128    34.515  88.645 108.489 111.876
 REC. CELL PARAMETERS   0.041644  0.033736  0.030727  83.835  70.590  67.351
 COORDINATES OF UNIT CELL A-AXIS   -14.250   -14.268    18.603
 COORDINATES OF UNIT CELL B-AXIS    14.443    25.574    13.023
 COORDINATES OF UNIT CELL C-AXIS   -21.606    21.484   -16.214
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.032    -0.003     0.027
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.035
 LAB COORDINATES OF ROTATION AXIS -0.999993  0.001823  0.003120
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000996  0.005064  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.34   1244.09
 DETECTOR ORIGIN (PIXELS) AT                     1265.03   1237.43
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.33
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     470 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   49     32     57     67     53     27     52     65     54
 SIGMAB (degree)       0.033  0.033  0.034  0.033  0.032  0.032  0.031  0.032  0.032
 SIGMAR (degree)       0.028  0.026  0.027  0.028  0.027  0.024  0.026  0.029  0.028


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    470 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  1  1  0  1  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  1  5 11 15 12  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  4  3  2  1  0  0  0     0  0  0  2  9 20 26 19 11  4  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  1  6 12 15 11  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  1  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  4  7  9  7  3  1  0  0  0     0  0  0  1  6 10 12  8  4  1  0  0  0     0  0  0  1  3  6  7  5  2  0  0  0  0
     0  0  0  3 17 38 49 35 18  7  3  1  0     0  0  1  6 27 52 62 41 20  8  3  1  0     0  0  0  4 17 30 35 20 11  4  1  0  0
     0  0  1  7 30 67 82 57 32 13  5  2  0     0  0  1 12 48 94100 66 34 14  5  2  0     0  0  1  7 28 54 54 35 17  8  3  1  0
     0  0  0  4 19 42 48 33 16  7  3  1  0     0  0  1  7 30 60 61 38 19  8  3  1  0     0  0  1  4 18 34 33 22 10  4  2  0  0
     0  0  0  1  5  9  9  5  3  1  0  0  0     0  0  0  2  7 13 13  7  3  2  0  0  0     0  0  0  2  5  8  7  4  2  1  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  8  4  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 13 12  7  4  2  1  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  4  8  8  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    638 REFLECTION INTENSITIES INTEGRATED


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
   461   0  1.081  5179836    0    128      81     7  0.04059  0.07578
   462   0  1.082  5180216    0    127      75     4  0.04083  0.08065
   463   0  1.080  5180620    0    128      68     3  0.04253  0.07326
   464   0  1.078  5180904    0    126      66     7  0.04280  0.08474
   465   0  1.078  5179710    0    125      79     3  0.04350  0.10593
   466   0  1.077  5179439    0    127      78     1  0.04130  0.11218
   467   0  1.079  5179934    0    128      75     4  0.04395  0.08881
   468   0  1.083  5179963    0    131      74     8  0.04106  0.06540
   469   0  1.083  5179654    0    129      77     8  0.04071  0.05958
   470   0  1.082  5179054    0    130      88     4  0.04494  0.08901
   471   0  1.081  5180410    0    113      71     5  0.04318  0.07815
   472   0  1.078  5181154    0    111      62     6  0.04401  0.07087
   473   0  1.078  5180584    0    113      63     6  0.04325  0.06368
   474   0  1.079  5180123    0    120      71     6  0.04130  0.06369
   475   0  1.080  5179823    0    131      81     1  0.04374  0.09456
   476   0  1.076  5180545    0    133      69     1  0.04105  0.08983
   477   0  1.077  5180555    0    122      67     5  0.04295  0.06791
   478   0  1.078  5180929    0    120      66     4  0.04290  0.07498
   479   0  1.076  5180141    0    120      74     2  0.04305  0.09649
   480   0  1.078  5180339    0    128      72     4  0.04026  0.08151

   475 OUT OF    476 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.51
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.08
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.436    32.123    34.508  88.665 108.490 111.895
 REC. CELL PARAMETERS   0.041652  0.033747  0.030735  83.804  70.578  67.323
 COORDINATES OF UNIT CELL A-AXIS   -14.251   -14.274    18.598
 COORDINATES OF UNIT CELL B-AXIS    14.442    25.570    13.020
 COORDINATES OF UNIT CELL C-AXIS   -21.604    21.471   -16.220
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.024    -0.003     0.036
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.035
 LAB COORDINATES OF ROTATION AXIS -0.999993  0.001823  0.003120
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000670  0.005035  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.31   1244.07
 DETECTOR ORIGIN (PIXELS) AT                     1265.43   1237.45
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.29
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     475 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   49     31     56     69     54     27     53     66     59
 SIGMAB (degree)       0.033  0.034  0.034  0.033  0.032  0.032  0.031  0.032  0.033
 SIGMAR (degree)       0.028  0.026  0.027  0.028  0.027  0.024  0.027  0.029  0.028


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    475 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  3  4  3  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  1  6 14 18 14  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  3  4  3  2  1  1  0  0     0  0  0  2  9 22 29 22 12  5  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  1  6 14 18 14  8  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  2  3  4  3  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0 -1  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  4 10 13  9  4  2  0  0  0     0  0  0  2  5 13 15  9  4  2  0  0  0     0  0  0  1  3  8  8  5  2  1  0  0  0
     0  0  1  5 19 46 56 40 20  8  3  1  0     0  0  1  8 30 61 65 43 21  8  3  1  0     0  0  1  5 18 32 31 22 10  4  1  0  0
     0  0  0  6 30 72 86 63 33 14  6  2  0     0  0  1 11 47 96100 69 36 15  6  2  0     0  0  1  7 27 51 50 33 17  7  3  1  0
     0  0  0  4 20 45 53 40 20  9  3  1  0     0  0  1  7 29 58 62 45 21  9  3  1  0     0  0  1  4 16 29 31 21 10  4  1  0  0
     0  0  0  1  6 11 13  9  4  2  0  0  0     0  0  0  2  8 15 15 10  4  2  0  0  0     0  0  0  1  5  8  7  5  2  1  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  1  2  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0 -1  0  0  0  0  0  0  0     0  0 -1  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  1  0  1  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  1  4  7  6  4  2  1  0  0  0     0  0  0  0  0  1  1  1  1  1  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  1  2  7 10 10  6  3  1  1  0 -1     0  0  0  0  1  1  1  0  1  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  4  7  7  5  2  1  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  2  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    632 REFLECTION INTENSITIES INTEGRATED


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
   481   0  1.078  5179611    0    142      80     3  0.04142  0.07994
   482   0  1.077  5179740    0    140      81     3  0.03907  0.07566
   483   0  1.077  5180324    0    130      74     3  0.03972  0.07245
   484   0  1.077  5180222    0    135      75     1  0.04149  0.09417
   485   0  1.080  5180676    0    133      70     1  0.04150  0.10554
   486   0  1.077  5180131    0    132      76     3  0.04136  0.08409
   487   0  1.079  5180149    0    132      75     3  0.03851  0.07255
   488   0  1.074  5179538    0    127      86     5  0.04128  0.07753
   489   0  1.074  5179865    0    125      77     8  0.04173  0.05711
   490   0  1.078  5180367    0    136      71     1  0.04192  0.10183
   491   0  1.078  5178914    0    134      84     2  0.04417  0.08225
   492   0  1.081  5179070    0    136      91     3  0.04292  0.09338
   493   0  1.079  5179467    0    138      80     0  0.04423  0.11167
   494   0  1.084  5179795    0    127      74     1  0.04521  0.11354
   495   0  1.086  5179345    0    126      79     3  0.04755  0.10855
   496   0  1.081  5179333    0    119      79     4  0.04281  0.11911
   497   0  1.080  5180653    0    114      65     2  0.04440  0.12126
   498   0  1.080  5180393    0    119      72     0  0.04742  0.13327
   499   0  1.079  5179836    0    122      76     0  0.04693  0.12819
   500   0  1.081  5179566    0    121      79     1  0.04651  0.10891

   487 OUT OF    487 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.52
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.432    32.137    34.507  88.648 108.487 111.885
 REC. CELL PARAMETERS   0.041653  0.033729  0.030735  83.830  70.590  67.341
 COORDINATES OF UNIT CELL A-AXIS   -14.251   -14.269    18.596
 COORDINATES OF UNIT CELL B-AXIS    14.439    25.584    13.027
 COORDINATES OF UNIT CELL C-AXIS   -21.601    21.473   -16.218
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.019     0.003     0.032
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.035
 LAB COORDINATES OF ROTATION AXIS -0.999993  0.001823  0.003120
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000727  0.005171  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.28   1244.06
 DETECTOR ORIGIN (PIXELS) AT                     1265.33   1237.26
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.29
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     487 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   54     34     63     66     55     31     52     59     60
 SIGMAB (degree)       0.033  0.033  0.034  0.033  0.032  0.032  0.032  0.032  0.032
 SIGMAR (degree)       0.029  0.027  0.028  0.028  0.026  0.023  0.025  0.027  0.028


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    487 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  1  4  9 11  8  4  2  1  0  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  1  6 13 16 13  7  3  1  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  1  3  8 10  7  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  1  3  7  8  6  3  1  0  0  0     0  0  0  2  6 13 13  9  4  1  1  0  0     0  0  0  1  4  9  9  6  2  1  0  0  0
     0  0  0  3 14 32 41 29 14  6  2  1  0     0  0  1  6 29 58 66 43 21  8  3  1  0     0  0  1  5 22 38 41 26 13  5  2  0  0
     0  0  1  5 22 50 60 46 24 10  4  1  0     0  0  1 10 47 93100 70 35 15  6  2  0     0  0  1  9 35 62 63 42 21  9  3  1  0
     0  0  1  3 14 32 35 26 13  6  2  1  0     0  0  1  6 29 59 59 42 20  8  3  1  0     0  0  1  6 22 40 37 25 12  5  2  0  0
     0  0  0  1  3  7  7  6  3  1  0  0  0     0  0  0  2  7 14 13 10  4  2  0  0  0     0  0  0  2  5  9  9  6  3  1  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  1  2  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  3  3  2  1  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 11 11  7  3  1  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  1  1  0  0  0  1  0  0  0
     0  0  1  3 11 17 16 11  5  2  1  0  0     0  0  1  1  2  2  2  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  2  6 11 10  7  3  1  0  0  0     0  0  1  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  2  3  2  1  1  0  0  0  0     0  1  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    640 REFLECTION INTENSITIES INTEGRATED


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
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  0  5 10 13 10  5  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  2  7 16 20 15  8  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  2  0  0  0  0  0     0  0  0  0  5 10 13 10  5  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  4  9 11  7  3  0  0  0  0     0  0  0  2  7 16 16  9  4  0  0  0  0     0  0  0  0  5 10 10  5  3  0  0  0  0
     0  0  0  5 18 38 46 33 16  6  0  0  0     0  0  0  9 33 63 67 45 20  8  0  0  0     0  0  0  6 22 39 38 24 11  4  0  0  0
     0  0  0  7 29 58 68 51 26 10  0  0  0     0  0  2 13 51 94100 70 33 13  0  0  0     0  0  2  9 34 58 57 38 17  6  0  0  0
     0  0  0  4 18 37 44 32 15  6  0  0  0     0  0  0  8 32 62 67 45 20  7  0  0  0     0  0  0  6 21 38 38 24 11  4  0  0  0
     0  0  0  2  4  9 11  7  3  0  0  0  0     0  0  0  3  8 17 17 11  5  2  0  0  0     0  0  0  3  5 10 10  6  3  0  0  0  0
     0  0  0  0  0  0  3  0  0  0  0  0  0     0  0  0  0  0  4  4  3  0  0  0  0  0     0  0  0  0  0  0  3  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 10  9  5  3  0  0  0  0     0  0  0  0  0  0  5  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3 10 15 14  9  4  2  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 10  9  6  2  0  0  0  0     0  0  0  0  0  4  4  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  2 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  4  6  7  5  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  4  9 10  7  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  4  6  7  5  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  2  0  0  0  0  0  0     0  0  0  0  0  2  2  0  0  0  0  0  0
     0  0  0  3  5  9  9  5  3  0  0  0  0     0  0  0  3  9 19 20 12  6  3  0  0  0     0  0  0  3  8 17 17 10  5  3  0  0  0
     0  0  0  4 15 30 32 22 11  4  0  0  0     0  0  0  8 34 68 71 47 22  9  0  0  0     0  0  0  7 30 61 62 38 18  8  0  0  0
     0  0  0  4 20 41 46 32 16  7  0  0  0     0  0  0 10 47 95100 68 33 14  0  0  0     0  0  0  9 42 84 84 55 26 11  0  0  0
     0  0  0  3 14 27 30 20 10  4  0  0  0     0  0  0  7 30 60 63 41 19  8  0  0  0     0  0  0  6 25 50 50 31 14  6  0  0  0
     0  0  0  3  6  9  9  6  3  0  0  0  0     0  0  0  3 11 18 17 11  5  2  0  0  0     0  0  0  2  7 13 13  8  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  4  3  0  0  0  0  0  0     0  0  0  0  0  2  2  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  4  4  4  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  9 16 15  9  4  3  0  0  0     0  0  0  0  3  4  4  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3 12 22 20 13  6  4  0  0  0     0  0  0  0  3  4  4  4  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 12 12  7  4  3  0  0  0     0  0  0  0  0  4  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  4  3  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  3 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  6  9  7  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  2  0  0  0  0  0     0  0  0  3  6 12 15 12  6  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  2  0  0  0  0  0  0     0  0  0  2  6 11 12  9  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  4  4  3  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  5  7  5  3  0  0  0  0     0  0  0  0  3 10 13  9  4  2  0  0  0     0  0  0  0  3  8 10  6  3  0  0  0  0
     0  0  0  3 10 25 33 25 13  6  0  0  0     0  0  0  4 20 47 58 42 22  9  0  0  0     0  0  0  3 16 36 42 28 14  6  0  0  0
     0  0  0  5 22 47 58 44 23 10  0  0  0     0  0  2  9 41 86100 72 37 16  6  0  0     0  0  0  7 31 62 69 48 24 10  0  0  0
     0  0  0  6 23 42 46 32 16  7  0  0  0     0  0  2 11 39 72 77 53 25 11  0  0  0     0  0  2  8 27 50 51 34 16  7  0  0  0
     0  0  0  4 11 17 16 10  5  0  0  0  0     0  0  0  6 18 27 25 16  7  3  0  0  0     0  0  0  4 11 17 16 10  4  0  0  0  0
     0  0  0  0  3  3  3  0  0  0  0  0  0     0  0  0  0  3  4  4  2  0  0  0  0  0     0  0  0  0  2  3  3  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  5 10 11  7  4  2  0  0  0     0  0  0  0  0  3  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  9 18 18 12  6  3  0  0  0     0  0  0  0  2  3  3  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  8 14 13  9  4  3  0  0  0     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  5  4  3  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  4 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  4  4  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  0  0  0  0  0     0  0  0  0  5 12 18 16  9  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  4  3  3  0  0  0  0     0  0  0  2  8 18 25 20 11  5  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  2  0  0  0  0  0  0     0  0  0  2  5 11 13  9  5  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  2  3  3  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  9 14 11  6  2  0  0  0     0  0  0  0  5 13 17 12  6  2  0  0  0     0  0  0  0  3  8  9  6  3  0  0  0  0
     0  0  0  4 17 41 56 44 24 11  0  0  0     0  0  0  7 28 60 71 51 27 12  5  0  0     0  0  0  5 18 34 37 25 13  6  0  0  0
     0  0  0  7 29 62 78 59 31 14  6  0  0     0  0  0 12 46 90100 70 36 16  6  0  0     0  0  0  8 29 52 53 36 18  8  0  0  0
     0  0  0  5 18 36 41 29 15  6  0  0  0     0  0  0  8 29 53 55 36 18  8  0  0  0     0  0  0  6 18 31 31 19  9  4  0  0  0
     0  0  0  3  5  8  8  5  3  0  0  0  0     0  0  0  3  8 12 11  7  4  0  0  0  0     0  0  0  0  5  7  6  4  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  5  9  9  6  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  8 13 13  8  4  2  0  0  0     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  5  8  7  4  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
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
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  2  0  0  0  0  0     0  0  0  0  3  9 14 12  7  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  0  0  0  0  0     0  0  0  0  5 15 22 19 11  5  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  3  8 11  9  5  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  2  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  6  9  7  3  0  0  0  0     0  0  0  0  4 10 11  8  3  2  0  0  0     0  0  0  0  3  5  5  3  0  0  0  0  0
     0  0  0  3 14 35 48 37 20  8  0  0  0     0  0  0  6 27 56 63 42 21  9  0  0  0     0  0  0  4 17 30 29 18  9  4  0  0  0
     0  0  0  5 24 59 77 58 31 14  7  0  0     0  0  0 11 47 94100 67 34 15  6  0  0     0  0  0  8 29 49 46 30 14  6  0  0  0
     0  0  0  3 12 31 38 27 14  6  0  0  0     0  0  0  5 24 49 50 32 15  7  0  0  0     0  0  0  4 14 25 24 14  7  3  0  0  0
     0  0  0  0  3  4  5  3  2  0  0  0  0     0  0  0  0  4  8  7  4  2  0  0  0  0     0  0  0  0  3  4  3  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  4  6  6  3  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  7 10  9  5  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  4  5  4  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
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
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  5  7  6  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  8 12  9  5  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  2  5  7  5  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  5  6  4  2  0  0  0  0     0  0  0  3  5 11 10  6  3  0  0  0  0     0  0  0  0  3  5  4  3  0  0  0  0  0
     0  0  0  3 10 28 35 24 12  5  0  0  0     0  0  0  6 30 61 59 34 16  7  0  0  0     0  0  0  5 19 30 24 13  6  3  0  0  0
     0  0  0  3 18 47 56 40 20  9  0  0  0     0  0  3 12 54100 94 58 28 12  0  0  0     0  0  2 10 32 49 41 24 11  5  0  0  0
     0  0  0  2 10 26 31 22 10  5  0  0  0     0  0  0  7 31 57 54 32 15  6  0  0  0     0  0  0  6 19 29 24 14  6  3  0  0  0
     0  0  0  0  2  5  5  3  2  0  0  0  0     0  0  0  2  6 12 10  6  3  0  0  0  0     0  0  0  0  4  6  5  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  4  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  7  6  3  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  4  5  4  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
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
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  2  6  9  7  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  0  0  0  0  0     0  0  0  0  5 13 19 16  9  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  2  0  0  0  0  0     0  0  0  0  3  9 13 10  5  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  4  3  2  0  0  0  0     0  0  0  0  3  6  6  4  2  0  0  0  0     0  0  0  0  2  3  3  2  0  0  0  0  0
     0  0  0  0  8 25 34 25 12  5  0  0  0     0  0  0  3 19 45 49 31 15  6  0  0  0     0  0  0  3 13 25 24 14  7  3  0  0  0
     0  0  0  3 21 55 69 53 27 12  0  0  0     0  0  0  8 43 94100 68 33 14  6  0  0     0  0  0  7 29 53 51 33 16  7  0  0  0
     0  0  0  3 14 36 47 33 17  7  0  0  0     0  0  0  6 28 62 67 44 21  9  0  0  0     0  0  0  5 19 36 35 22 10  4  0  0  0
     0  0  0  0  3  6  8  5  3  0  0  0  0     0  0  0  3  5 11 12  7  3  0  0  0  0     0  0  0  0  4  7  6  4  2  0  0  0  0
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
     0  0  0  3  7 12 10  6  3  2  0  0  0     0  0  0  0  0  3  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  5  9  8  5  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
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
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  2  0  0  0  0  0     0  0  0  0  4 10 13 10  5  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  4  4  3  0  0  0  0     0  0  0  2  8 19 25 20 11  5  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  2  0  0  0  0     0  0  0  0  5 13 17 14  7  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  4  3  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  7  7  4  2  0  0  0  0     0  0  0  0  5 10 10  6  3  0  0  0  0     0  0  0  0  3  6  6  3  0  0  0  0  0
     0  0  0  4 16 35 42 31 15  7  0  0  0     0  0  0  6 25 51 56 38 18  8  0  0  0     0  0  0  4 16 30 31 20 10  4  0  0  0
     0  0  0  6 29 64 78 59 32 14  6  0  0     0  0  0 10 46 92100 71 37 16  6  0  0     0  0  0  7 29 54 55 38 19  8  0  0  0
     0  0  0  4 18 42 53 39 20  9  0  0  0     0  0  0  7 28 60 68 47 23 10  0  0  0     0  0  0  5 18 35 37 24 12  5  0  0  0
     0  0  0  0  3  8 10  7  4  0  0  0  0     0  0  0  2  5 12 13  8  4  2  0  0  0     0  0  0  0  3  7  7  5  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  4  8  8  5  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  8 14 14  9  5  2  0  0  0     0  0  0  0  0  3  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  5  9  9  6  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  2  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  9 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  4  4  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  2  0  0  0  0  0     0  0  0  3  6 12 14 10  5  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  2  7 15 18 14  7  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  3  7  9  7  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  8 14 13  7  3  0  0  0  0     0  0  0  4 12 21 19 11  5  2  0  0  0     0  0  0  3  7 13 12  7  3  0  0  0  0
     0  0  0  6 24 47 51 34 16  7  0  0  0     0  0  2 10 40 75 78 50 24 10  0  0  0     0  0  0  7 27 49 49 31 15  6  0  0  0
     0  0  0  5 26 55 65 47 24 11  0  0  0     0  0  0  9 44 91100 70 35 15  6  0  0     0  0  0  7 31 62 65 44 21  9  0  0  0
     0  0  0  3  9 25 31 23 12  5  0  0  0     0  0  0  3 17 41 48 33 17  7  0  0  0     0  0  0  3 12 28 32 21 10  4  0  0  0
     0  0  0  0  0  4  5  4  2  0  0  0  0     0  0  0  0  3  6  8  6  3  0  0  0  0     0  0  0  0  2  4  5  4  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  4  4  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  7 13 13  8  4  3  0  0  0     0  0  0  0  3  3  3  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  9 17 16 11  5  3  0  0  0     0  0  0  0  3  3  3  4  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  7  8  5  3  2  0  0  0     0  0  0  0  0  2  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  2  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 STANDARD DEVIATION (DEGREES) OF BEAM DIVERGENCE
 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1    0.023  0.024  0.024  0.024  0.023  0.023  0.023  0.024  0.024
     2    0.024  0.024  0.025  0.024  0.022  0.023  0.023  0.023  0.024
     3    0.024  0.024  0.024  0.024  0.023  0.023  0.023  0.023  0.023
     4    0.024  0.024  0.025  0.024  0.023  0.023  0.023  0.024  0.023
     5    0.023  0.023  0.024  0.023  0.023  0.023  0.023  0.023  0.023
     6    0.023  0.024  0.024  0.024  0.023  0.023  0.022  0.023  0.023
     7    0.023  0.023  0.024  0.024  0.022  0.022  0.023  0.023  0.024
     8    0.024  0.024  0.024  0.024  0.023  0.023  0.022  0.023  0.023
     9    0.023  0.024  0.024  0.024  0.022  0.023  0.022  0.023  0.024
    10    0.023  0.023  0.024  0.023  0.022  0.023  0.022  0.023  0.023
    11    0.023  0.024  0.025  0.023  0.022  0.022  0.022  0.023  0.023
    12    0.023  0.024  0.024  0.024  0.023  0.024  0.022  0.023  0.023
    13    0.024  0.024  0.025  0.024  0.023  0.022  0.022  0.023  0.023
    14    0.023  0.023  0.024  0.024  0.023  0.023  0.022  0.023  0.023
    15    0.024  0.024  0.025  0.024  0.023  0.023  0.022  0.023  0.023
    16    0.023  0.023  0.024  0.024  0.023  0.023  0.022  0.023  0.023
    17    0.023  0.024  0.025  0.024  0.023  0.023  0.022  0.023  0.023
    18    0.023  0.024  0.024  0.024  0.022  0.023  0.022  0.023  0.023
    19    0.023  0.024  0.024  0.024  0.023  0.023  0.022  0.023  0.023
    20    0.023  0.024  0.025  0.024  0.023  0.023  0.022  0.023  0.023
    21    0.023  0.024  0.024  0.024  0.023  0.023  0.022  0.023  0.023
    22    0.024  0.024  0.024  0.024  0.023  0.023  0.022  0.023  0.023
    23    0.023  0.024  0.024  0.023  0.023  0.023  0.022  0.023  0.023
    24    0.023  0.024  0.024  0.024  0.023  0.023  0.022  0.023  0.023
    25    0.023  0.024  0.024  0.024  0.023  0.023  0.022  0.023  0.023


 STANDARD DEVIATION (DEGREES) OF REFLECTING RANGE
 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1    0.027  0.024  0.026  0.029  0.026  0.020  0.025  0.028  0.026
     2    0.027  0.031  0.028  0.027  0.025  0.021  0.026  0.027  0.027
     3    0.026  0.025  0.026  0.027  0.026  0.022  0.025  0.028  0.026
     4    0.027  0.024  0.026  0.028  0.026  0.022  0.027  0.031  0.028
     5    0.027  0.026  0.027  0.028  0.026  0.022  0.026  0.027  0.027
     6    0.027  0.024  0.026  0.028  0.025  0.022  0.025  0.028  0.026
     7    0.026  0.026  0.026  0.027  0.025  0.022  0.025  0.027  0.027
     8    0.026  0.024  0.026  0.028  0.025  0.021  0.025  0.028  0.026
     9    0.026  0.027  0.027  0.027  0.025  0.022  0.025  0.027  0.026
    10    0.026  0.025  0.026  0.027  0.026  0.021  0.025  0.027  0.026
    11    0.026  0.026  0.027  0.027  0.026  0.020  0.025  0.028  0.026
    12    0.026  0.027  0.026  0.027  0.026  0.022  0.025  0.027  0.027
    13    0.027  0.025  0.027  0.030  0.027  0.023  0.025  0.028  0.027
    14    0.027  0.026  0.027  0.029  0.028  0.024  0.027  0.028  0.027
    15    0.027  0.027  0.027  0.028  0.027  0.022  0.027  0.029  0.027
    16    0.027  0.026  0.027  0.027  0.026  0.023  0.026  0.029  0.028
    17    0.029  0.030  0.028  0.029  0.027  0.024  0.025  0.027  0.029
    18    0.027  0.027  0.026  0.027  0.025  0.022  0.025  0.027  0.028
    19    0.027  0.026  0.027  0.028  0.027  0.022  0.026  0.030  0.028
    20    0.027  0.030  0.028  0.028  0.027  0.023  0.025  0.027  0.027
    21    0.027  0.037  0.030  0.027  0.025  0.022  0.025  0.028  0.032
    22    0.027  0.031  0.028  0.027  0.025  0.021  0.025  0.027  0.028
    23    0.028  0.029  0.027  0.028  0.028  0.024  0.028  0.031  0.028
    24    0.028  0.025  0.026  0.030  0.028  0.023  0.028  0.032  0.029
    25    0.027  0.031  0.028  0.027  0.025  0.021  0.024  0.027  0.028


 ***** SUGGESTED VALUES FOR INPUT PARAMETERS *****
 BEAM_DIVERGENCE=   0.23246  BEAM_DIVERGENCE_E.S.D.=   0.02325
 REFLECTING_RANGE=  0.18509  REFLECTING_RANGE_E.S.D.=  0.02644



 NUMBER OF STRONG REFLECTIONS USED FOR LEARNING REFERENCE PROFILES

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1       51     36     63     70     59     34     58     69     63
     2       47     32     57     67     57     27     48     63     54
     3       47     32     60     63     50     27     53     69     60
     4       49     31     55     62     58     31     54     65     57
     5       47     29     55     63     51     30     60     67     52
     6       47     34     60     65     54     28     53     70     65
     7       45     31     56     69     53     29     51     63     49
     8       49     32     59     71     56     29     57     70     64
     9       46     35     62     70     55     29     57     66     54
    10       48     25     53     71     52     27     51     67     57
    11       45     28     58     69     54     29     57     63     51
    12       48     38     57     69     50     25     49     65     63
    13       52     31     59     69     54     29     59     66     53
    14       50     34     53     66     52     25     52     68     65
    15       44     34     56     63     52     29     53     69     58
    16       50     26     53     67     56     31     55     70     56
    17       48     32     60     68     54     25     50     64     56
    18       47     32     55     64     58     30     54     65     60
    19       48     31     58     67     56     33     56     68     53
    20       48     27     55     72     61     30     56     67     56
    21       47     30     61     68     55     29     56     66     56
    22       49     29     55     62     54     31     57     70     56
    23       49     32     57     67     53     27     52     65     54
    24       49     31     56     69     54     27     53     66     59
    25       54     34     63     66     55     31     52     59     60
 margin    1203    789   1436   1677   1363    721   1353   1660   1431



 NUMBER OF GRID POINTS IN REFERENCE PROFILES WITH  SIGNAL >   2.% OF MAXIMUM

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1      142    131    144    149    121     98    114    138    128
     2      136    140    152    137    115    107    115    136    137
     3      136    133    139    132    119    107    114    134    125
     4      141    131    146    143    119    109    122    143    136
     5      137    129    136    136    118    106    117    131    127
     6      136    131    142    139    119    107    114    132    127
     7      128    132    143    140    112    107    118    130    128
     8      135    127    140    140    118    107    116    134    124
     9      137    129    149    134    113    109    113    129    131
    10      131    124    137    136    118    109    115    123    123
    11      134    135    151    135    119    100    107    136    129
    12      133    132    140    139    125    111    118    128    127
    13      139    130    149    144    123    109    111    134    127
    14      133    131    142    139    126    111    115    130    124
    15      139    140    146    142    121    104    118    138    133
    16      141    131    143    144    120    112    115    133    128
    17      139    141    151    144    122    114    114    126    131
    18      134    141    141    140    114    104    113    137    131
    19      136    138    147    144    123    108    117    145    128
    20      132    141    154    143    118    110    110    125    129
    21      131    142    149    132    114    108    111    126    140
    22      142    139    152    144    113    111    112    127    133
    23      138    136    142    135    124    108    119    140    132
    24      138    135    145    149    123    107    121    143    142
    25      141    149    156    141    120    105    113    123    131
 margin     168    170    178    160    137    126    137    157    168

 MINIMUM INTEGRATION VOLUME IS   98 GRID POINTS



 CORRELATIONS BETWEEN REFERENCE PROFILES WITH CENTRAL PROFILE #1

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1     1.00   0.96   0.95   0.98   0.98   0.95   0.94   0.97   0.97
     2     1.00   0.88   0.95   0.98   0.97   0.97   0.96   0.98   0.95
     3     1.00   0.96   0.95   0.99   0.98   0.97   0.98   0.98   0.98
     4     1.00   0.96   0.97   0.98   0.98   0.98   0.97   0.96   0.94
     5     1.00   0.87   0.96   0.97   0.94   0.95   0.97   0.98   0.94
     6     1.00   0.97   0.97   0.97   0.96   0.97   0.95   0.96   0.96
     7     1.00   0.93   0.96   0.98   0.99   0.97   0.98   0.99   0.95
     8     1.00   0.95   0.96   0.97   0.97   0.97   0.97   0.97   0.94
     9     1.00   0.91   0.96   0.98   0.98   0.97   0.98   0.98   0.96
    10     1.00   0.97   0.97   0.98   0.97   0.97   0.98   0.98   0.97
    11     1.00   0.93   0.95   0.99   0.98   0.97   0.95   0.97   0.97
    12     1.00   0.93   0.97   0.97   0.97   0.95   0.98   0.99   0.95
    13     1.00   0.97   0.97   0.98   0.98   0.98   0.98   0.98   0.96
    14     1.00   0.96   0.97   0.96   0.96   0.98   0.98   0.99   0.97
    15     1.00   0.88   0.93   0.98   0.97   0.97   0.96   0.98   0.94
    16     1.00   0.96   0.96   0.99   0.98   0.98   0.98   0.96   0.96
    17     1.00   0.85   0.95   0.97   0.93   0.95   0.97   0.99   0.94
    18     1.00   0.95   0.95   0.98   0.97   0.98   0.97   0.96   0.98
    19     1.00   0.95   0.97   0.97   0.96   0.98   0.96   0.97   0.93
    20     1.00   0.91   0.94   0.98   0.97   0.98   0.97   0.98   0.96
    21     1.00   0.82   0.95   0.98   0.96   0.97   0.99   0.99   0.92
    22     1.00   0.88   0.92   0.98   0.97   0.97   0.98   0.99   0.95
    23     1.00   0.87   0.93   0.98   0.96   0.97   0.96   0.97   0.97
    24     1.00   0.96   0.96   0.97   0.95   0.96   0.96   0.97   0.96
    25     1.00   0.89   0.94   0.98   0.98   0.97   0.97   0.98   0.95
 margin    1.00   0.94   0.97   0.98   0.98   0.98   0.98   0.99   0.97



 CORRELATIONS BETWEEN REFERENCE PROFILES WITH RUN-AVERAGED PROFILES

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1     0.99   0.96   0.97   0.98   0.99   0.99   1.00   1.00   0.99
     2     0.99   0.96   0.97   0.98   0.99   0.99   1.00   0.99   0.98
     3     1.00   1.00   1.00   1.00   0.99   0.98   0.99   1.00   1.00
     4     0.99   0.97   0.99   1.00   1.00   0.99   0.98   0.98   0.97
     5     1.00   0.99   1.00   0.99   0.99   0.98   0.99   1.00   0.99
     6     1.00   0.97   0.99   1.00   0.99   0.99   0.99   0.99   0.99
     7     0.99   0.99   0.99   0.99   0.98   0.98   0.99   0.99   0.99
     8     0.99   0.98   1.00   1.00   1.00   1.00   1.00   0.99   0.99
     9     0.99   0.98   0.99   0.99   0.98   0.97   0.99   0.99   0.99
    10     0.99   1.00   1.00   1.00   1.00   0.99   0.98   0.99   0.99
    11     1.00   1.00   1.00   1.00   1.00   0.99   0.99   0.99   0.99
    12     0.99   0.99   1.00   1.00   1.00   0.96   0.97   0.98   0.99
    13     1.00   0.97   0.98   0.99   0.99   1.00   1.00   0.99   0.99
    14     1.00   0.99   0.99   0.99   0.99   0.99   0.99   1.00   1.00
    15     0.99   0.99   1.00   1.00   1.00   0.98   0.98   0.99   0.99
    16     1.00   1.00   1.00   0.99   0.99   0.99   1.00   0.99   0.99
    17     0.99   0.97   0.99   0.99   0.98   0.97   1.00   0.99   0.97
    18     1.00   1.00   1.00   0.99   0.99   0.99   1.00   0.99   0.99
    19     0.99   0.98   0.99   1.00   1.00   0.99   0.97   0.98   0.97
    20     0.99   0.99   1.00   1.00   1.00   1.00   1.00   1.00   0.99
    21     0.98   0.90   0.97   0.99   0.99   0.99   0.97   0.96   0.91
    22     1.00   0.98   0.98   0.99   1.00   1.00   0.99   1.00   0.99
    23     0.99   0.99   1.00   0.99   0.98   0.95   0.96   0.98   0.99
    24     0.99   0.97   0.97   0.98   0.98   0.98   0.98   0.98   0.98
    25     0.99   0.98   0.98   0.99   0.99   0.99   1.00   1.00   0.99
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

       -1... 0      -0.141      0.061        -2.          6
        0... 1       0.259      0.111         5.         56
        1... 2       0.486      0.104        14.         80
        2... 3       0.611      0.087        33.        161
        3... 6       0.719      0.085        83.        599
        6... 9       0.754      0.085       197.        862
        9...12       0.766      0.083       361.       1018
       12...15       0.771      0.083       572.       1120
       15...18       0.777      0.080       845.       1059
       18...21       0.786      0.079      1197.        971
       21...24       0.788      0.083      1595.        892
       24...27       0.792      0.080      2083.        808
       27...30       0.797      0.081      2643.        710
       30...33       0.797      0.086      3367.        616
       33...36       0.793      0.077      4364.        572
       36...39       0.792      0.089      5678.        577
       39...42       0.798      0.083      7980.        613
       42...45       0.791      0.084     11796.        708
       45...48       0.800      0.091     21225.        969
       48...51       0.802      0.108     43161.        728

   13125 REFLECTIONS SAVED ON FILE "INTEGRATE.HKL"


 total elapsed wall-clock time for INTEGRATE      16.0 sec

 TASK     cpu time (sec)    elapsed wall-clock time (sec)
   1         255.3                    15.9
