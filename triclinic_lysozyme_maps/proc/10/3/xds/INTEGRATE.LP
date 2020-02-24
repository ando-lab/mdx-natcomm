
 ***** INTEGRATE ***** (VERSION Mar 15, 2019  BUILT=20190806)   6-Jan-2020

 INPUT PARAMETER VALUES
 ----------------------
 MAXIMUM_NUMBER_OF_PROCESSORS=   0  (0: automatic choice of OpenMP threads)
 number of OpenMP threads used  20
 MAXIMUM_NUMBER_OF_JOBS=   0  (0:automatic choice of forked main tasks)
 number of forked INTEGRATE tasks   1
 CLUSTER_NODES=
 SECONDS=     0
 NAME_TEMPLATE_OF_DATA_FRAMES=~/data/lysozyme/nitrate/ubatch4/lys_nitr_10_3_????.cbf CBF
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
 UNIT_CELL_CONSTANTS=    27.43    32.15    34.53  88.660 108.471 111.869
 ROTATION_AXIS= -0.99999 -0.00410  0.00137
 OSCILLATION_RANGE=  0.100000 DEGREES
 STARTING_ANGLE= 340.000 STARTING_FRAME=       1
 X-RAY_WAVELENGTH=   0.97700 ANGSTROM
 INCIDENT_BEAM_DIRECTION=  0.000567  0.005055  1.023529
 NUMBER OF DETECTOR SEGMENTS   1
 NX=  2463 NY=  2527 QX= 0.17200 QY= 0.17200
 ORGX=   1265.61  ORGY=   1237.41
 DETECTOR_DISTANCE=   231.207
 DIRECTION_OF_DETECTOR_X-AXIS=   1.00000   0.00000   0.00000
 DIRECTION_OF_DETECTOR_Y-AXIS=   0.00000   1.00000   0.00000
 SEGMENT=    1 2463    1 2527  ! x1,x2,y1,y2 (pixel units)
    SEGMENT_ORGX=      0.00   SEGMENT_ORGY=      0.00
    SEGMENT_DISTANCE=     0.00
    DIRECTION_OF_SEGMENT_X-AXIS=  1.000000  0.000000  0.000000
    DIRECTION_OF_SEGMENT_Y-AXIS=  0.000000  1.000000  0.000000
 RESOLUTION RANGE RECORDED BY DETECTOR (ANGSTROM)    44.448     1.211
 NUMBER OF TRUSTED DETECTOR PIXELS    5212633
 MEAN CONTENTS OF TRUSTED PIXELS IN BACKGROUND TABLE    142.541
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
     1   0  1.060  5204266    0    149      86     3  0.04578  0.10733
     2   0  1.021  5203976    0    145      94     6  0.04312  0.08412
     3   0  1.021  5204264    0    139      89     2  0.04594  0.11701
     4   0  1.021  5204135    0    136      91     1  0.04350  0.11375
     5   0  1.018  5204311    0    134      87     6  0.04486  0.07462
     6   0  1.022  5204551    0    140      87     3  0.04413  0.08565
     7   0  1.020  5204556    0    150      89     4  0.04389  0.08833
     8   0  1.018  5203320    0    142     100     5  0.04373  0.08561
     9   0  1.017  5204112    0    143      88     4  0.04420  0.10508
    10   0  1.013  5204848    0    133      81     1  0.04502  0.14233
    11   0  1.013  5205580    0    118      75     3  0.04353  0.11202
    12   0  1.016  5206062    0    118      69     4  0.04445  0.10635
    13   0  1.010  5206305    0    118      63     6  0.04317  0.07352
    14   0  1.010  5206072    0    122      62     6  0.04204  0.06124
    15   0  1.009  5205520    0    120      69     7  0.04252  0.07023
    16   0  1.007  5205589    0    123      72     1  0.04341  0.13083
    17   0  1.008  5205694    0    122      73     4  0.04260  0.10880
    18   0  1.009  5205714    0    125      71     2  0.04625  0.15401
    19   0  1.008  5205407    0    120      78     1  0.04346  0.11966
    20   0  1.006  5205710    0    114      76     2  0.04772  0.12011

   493 OUT OF    495 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.61
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.14
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.426    32.131    34.510  88.636 108.464 111.899
 REC. CELL PARAMETERS   0.041658  0.033737  0.030727  83.847  70.618  67.332
 COORDINATES OF UNIT CELL A-AXIS   -14.363   -14.235    18.527
 COORDINATES OF UNIT CELL B-AXIS    14.403    25.573    13.074
 COORDINATES OF UNIT CELL C-AXIS   -21.500    21.492   -16.334
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.012    -0.012     0.020
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.038
 LAB COORDINATES OF ROTATION AXIS -0.999991 -0.004096  0.001370
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000958  0.005114  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.34   1243.99
 DETECTOR ORIGIN (PIXELS) AT                     1265.08   1237.28
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.17
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     493 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   46     33     56     71     56     30     57     69     55
 SIGMAB (degree)       0.035  0.035  0.035  0.034  0.033  0.034  0.034  0.035  0.035
 SIGMAR (degree)       0.031  0.030  0.029  0.030  0.028  0.025  0.029  0.031  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    493 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  1  2  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  1  4 10 11  8  4  2  1  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  3  3  2  1  0  0  1     0  0  0  2  6 14 17 14  8  3  1  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  1  3  7 10  8  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  7  5  2  1  0  0  0     0  0  0  2  5 11 11  6  3  1  0  0  0     0  0  0  1  3  7  7  4  2  0  0  0  0
     0  0  0  4 15 35 42 30 13  6  2  1  0     0  0  1  6 26 56 63 40 17  7  2  1  0     0  0  0  5 19 37 39 22 10  4  1  0  0
     0  0  1  5 24 54 65 51 26 11  4  1  0     0  0  1  8 44 89100 70 34 14  5  1  0     0  0  1  7 32 59 64 41 19  8  3  1  0
     0  0  0  2 12 28 38 28 15  6  2  1  0     0  0  0  4 22 48 58 40 20  8  3  1  0     0  0  0  4 16 32 36 23 11  5  1  0  0
     0  0  0  1  2  4  6  4  2  1  0  0  0     0  0  0  1  4  8 10  7  3  1  0  0  0     0  0  0  1  3  5  7  4  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  9  8  4  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  2  8 13 13  8  4  2  1  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  4  8  8  5  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  1  1  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    625 REFLECTION INTENSITIES INTEGRATED


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
    21   0  1.006  5205840    0    119      73     1  0.04716  0.12647
    22   0  1.008  5205827    0    114      69     1  0.04631  0.12244
    23   0  1.007  5205539    0    125      73     1  0.04748  0.14094
    24   0  1.006  5204566    0    126      78     3  0.04377  0.08094
    25   0  1.006  5204598    0    130      85     1  0.04309  0.10502
    26   0  1.006  5204654    0    134      84     1  0.04267  0.10915
    27   0  1.005  5203450    0    141      93     4  0.04063  0.09527
    28   0  1.003  5203197    0    144      98     1  0.04202  0.12238
    29   0  1.003  5203873    0    153      95     4  0.04328  0.11076
    30   0  1.002  5204232    0    142      90     3  0.04860  0.11654
    31   0  1.000  5204020    0    123      92     1  0.04643  0.13113
    32   0  1.003  5204282    0    134      86     1  0.04647  0.13777
    33   0  1.004  5204056    0    134      90     2  0.04903  0.12924
    34   0  1.003  5203380    0    142      97     2  0.04827  0.13870
    35   0  1.001  5204011    0    147      89     2  0.04332  0.12484
    36   0  0.998  5204120    0    137      84     4  0.04368  0.08973
    37   0  0.995  5204829    0    124      74     6  0.04353  0.07284
    38   0  0.998  5206209    0    123      65     1  0.04822  0.12975
    39   0  0.995  5205488    0    124      70     2  0.04666  0.13364
    40   0  0.998  5205699    0    129      76     4  0.04320  0.12070

   487 OUT OF    490 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.50
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.426    32.127    34.512  88.643 108.474 111.904
 REC. CELL PARAMETERS   0.041664  0.033744  0.030727  83.833  70.604  67.324
 COORDINATES OF UNIT CELL A-AXIS   -14.363   -14.230    18.531
 COORDINATES OF UNIT CELL B-AXIS    14.403    25.573    13.064
 COORDINATES OF UNIT CELL C-AXIS   -21.497    21.486   -16.349
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.000    -0.001     0.010
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.038
 LAB COORDINATES OF ROTATION AXIS -0.999991 -0.004096  0.001370
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000756  0.005204  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.35   1244.08
 DETECTOR ORIGIN (PIXELS) AT                     1265.35   1237.25
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.14
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     487 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   53     33     65     66     53     27     56     68     56
 SIGMAB (degree)       0.034  0.035  0.035  0.034  0.033  0.034  0.034  0.034  0.034
 SIGMAR (degree)       0.029  0.028  0.028  0.030  0.028  0.025  0.028  0.029  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    487 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  1  4  9 11  8  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  1  5 14 19 15  8  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  3  7 10  8  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  7  4  2  1  0  0  0     0  0  0  2  6 10 10  6  2  1  0  0  0     0  0  0  1  4  6  5  3  1  0  0  0  0
     0  0  0  3 15 33 39 28 13  5  2  0  0     0  0  1  6 26 54 58 38 17  7  2  1  0     0  0  0  4 17 32 32 20  9  4  1  0  0
     0  0  0  5 22 55 67 50 25 10  4  1  0     0  0  1  9 42 88100 71 34 13  5  1  0     0  0  1  8 28 52 56 39 18  7  2  0  0
     0  0  0  2 12 28 36 27 14  6  2  0  0     0  0  0  4 24 48 55 39 20  8  3  1  0     0  0  0  3 17 29 32 21 11  4  1  0  0
     0  0  0  0  2  4  5  4  2  1  0  0  0     0  0  0  1  4  8 10  7  3  1  0  0  0     0  0  0  1  3  5  6  4  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  6  4  2  1  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 11 11  7  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  7  4  2  1  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    643 REFLECTION INTENSITIES INTEGRATED


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
    41   0  0.998  5204688    0    126      80     3  0.03978  0.10979
    42   0  0.995  5204906    0    125      82     8  0.04318  0.07758
    43   0  0.995  5205788    0    134      73     6  0.04347  0.09991
    44   0  0.996  5205225    0    125      74     5  0.04378  0.07691
    45   0  0.994  5204959    0    134      85     3  0.04115  0.10037
    46   0  0.997  5205333    0    140      79     2  0.04360  0.10244
    47   0  0.999  5204341    0    140      89     5  0.04196  0.08068
    48   0  1.000  5205284    0    146      76     4  0.04413  0.07876
    49   0  0.995  5204464    0    153      83     6  0.04193  0.06275
    50   0  0.994  5203659    0    152      93     9  0.04075  0.06544
    51   0  0.993  5203580    0    153      94     7  0.04438  0.08243
    52   0  0.994  5203631    0    155      92     6  0.04401  0.07512
    53   0  0.993  5204264    0    136      89     3  0.04304  0.09106
    54   0  0.994  5204927    0    134      86     3  0.04470  0.10698
    55   0  0.989  5206004    0    120      71     1  0.04720  0.14415
    56   0  0.985  5205067    0    124      76     2  0.04554  0.12430
    57   0  0.984  5204730    0    128      84     1  0.04717  0.14719
    58   0  0.986  5204988    0    131      81     1  0.04808  0.12720
    59   0  0.983  5204592    0    134      86     0  0.04616  0.12915
    60   0  0.983  5204374    0    137      84     3  0.04435  0.08499

   495 OUT OF    495 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.59
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.418    32.130    34.516  88.644 108.466 111.935
 REC. CELL PARAMETERS   0.041683  0.033748  0.030723  83.824  70.609  67.292
 COORDINATES OF UNIT CELL A-AXIS   -14.363   -14.236    18.515
 COORDINATES OF UNIT CELL B-AXIS    14.407    25.574    13.066
 COORDINATES OF UNIT CELL C-AXIS   -21.502    21.490   -16.345
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.005    -0.003     0.019
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.038
 LAB COORDINATES OF ROTATION AXIS -0.999991 -0.004096  0.001370
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000533  0.005063  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.35   1244.08
 DETECTOR ORIGIN (PIXELS) AT                     1265.65   1237.43
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.15
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     495 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   51     33     53     66     58     32     53     68     63
 SIGMAB (degree)       0.034  0.035  0.035  0.034  0.033  0.034  0.034  0.035  0.034
 SIGMAR (degree)       0.030  0.027  0.029  0.030  0.029  0.025  0.027  0.031  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    495 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  1  2  3  3  2  1  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  2  7 14 15 12  6  2  1  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  3  4  3  2  1  0  0  0     0  0  0  2  9 20 26 21 11  5  2  1  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0     0  0  0  1  4 11 14 12  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  7  9  6  3  1  0  0  0     0  0  0  2  6  9  9  6  2  1  0  0  0     0  0  0  1  3  4  4  3  1  0  0  0  0
     0  0  1  4 18 42 45 31 16  6  2  1  0     0  0  1  6 27 57 55 35 16  7  2  0  0     0  0  1  4 15 27 24 14  6  2  1  0  0
     0  0  1  7 30 66 80 61 31 13  5  1  0     0  0  1 10 47 93100 69 33 13  5  1  0     0  0  1  6 24 45 43 28 13  5  2  0  0
     0  0  0  3 15 36 47 36 18  8  3  1  0     0  0  0  5 24 51 60 41 19  8  3  1  0     0  0  0  3 12 25 28 17  8  3  1  0  0
     0  0  0  0  2  7  9  6  3  1  0  0  0     0  0  0  1  4 11 12  8  4  1  0  0  0     0  0  0  1  2  6  5  4  2  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  5  4  3  1  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  9  8  5  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  5  5  3  1  1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    644 REFLECTION INTENSITIES INTEGRATED


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
    61   0  0.983  5205012    0    129      80     1  0.04417  0.10152
    62   0  0.983  5204988    0    129      79     4  0.04251  0.09123
    63   0  0.985  5205669    0    128      75     3  0.04264  0.09863
    64   0  0.985  5204181    0    129      91     0  0.04252  0.11767
    65   0  0.984  5204505    0    126      87     1  0.04476  0.12406
    66   0  0.983  5204875    0    128      86     2  0.04684  0.12743
    67   0  0.980  5205211    0    124      77     0  0.04629  0.14370
    68   0  0.980  5204739    0    147      85     2  0.04678  0.14225
    69   0  0.980  5203901    0    138      95     6  0.04595  0.09603
    70   0  0.978  5203304    0    152      98     8  0.04392  0.09122
    71   0  0.977  5203428    0    144     103     8  0.04247  0.08293
    72   0  0.977  5204116    0    139      92     4  0.04316  0.10787
    73   0  0.976  5205349    0    136      79     3  0.04440  0.09505
    74   0  0.976  5205306    0    131      82     3  0.04434  0.10435
    75   0  0.977  5206049    0    124      75     4  0.04479  0.09560
    76   0  0.976  5206254    0    126      71     3  0.04296  0.10949
    77   0  0.974  5204979    0    137      82     3  0.04206  0.12259
    78   0  0.974  5203947    0    139      86     7  0.04150  0.08503
    79   0  0.972  5203723    0    142      87     7  0.04444  0.08150
    80   0  0.971  5204704    0    142      82     5  0.04529  0.09324

   499 OUT OF    502 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.55
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.411    32.129    34.509  88.641 108.464 111.933
 REC. CELL PARAMETERS   0.041692  0.033749  0.030729  83.830  70.612  67.296
 COORDINATES OF UNIT CELL A-AXIS   -14.360   -14.234    18.508
 COORDINATES OF UNIT CELL B-AXIS    14.406    25.571    13.070
 COORDINATES OF UNIT CELL C-AXIS   -21.497    21.490   -16.338
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.011    -0.004     0.022
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.038
 LAB COORDINATES OF ROTATION AXIS -0.999991 -0.004096  0.001370
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000600  0.005040  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.29   1244.00
 DETECTOR ORIGIN (PIXELS) AT                     1265.50   1237.38
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.10
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     499 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   52     33     59     68     54     28     60     70     62
 SIGMAB (degree)       0.034  0.035  0.034  0.034  0.033  0.034  0.034  0.035  0.035
 SIGMAR (degree)       0.030  0.027  0.029  0.031  0.029  0.026  0.028  0.031  0.029


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    499 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  1  2  3  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  2  6 11 14 10  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  3  2  1  0  0  0     0  0  0  2  8 17 22 18 10  4  1  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  1  0  0  0     0  0  0  1  4  8 12 10  6  2  1  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4  7  8  4  2  0  0  0  0     0  0  0  1  5  9  9  5  2  1  0  0  0     0  0  0  1  3  4  4  3  1  0  0  0  0
     0  0  0  4 17 34 43 31 14  6  2  0  0     0  0  0  6 25 51 57 37 17  7  2  0  0     0  0  0  4 15 28 28 17  8  3  1  0  0
     0  0  1  5 26 59 73 58 28 12  4  1  0     0  0  1  9 45 92100 70 33 14  5  1  0     0  0  1  7 28 51 49 31 15  6  2  0  0
     0  0  1  3 14 30 43 33 17  7  3  1  0     0  0  1  6 24 51 59 41 20  8  3  1  0     0  0  0  4 16 30 30 19  9  4  1  0  0
     0  0  0  0  2  5  7  5  3  1  0  0  0     0  0  0  1  3  8 10  8  4  1  0  0  0     0  0  0  1  2  5  5  4  2  1  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0 -1  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  6  5  3  1  0  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  6 10  9  5  2  1  0  0  1     0  0  0  0  1  1  1  0  0  0  0  0  2     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  6  6  3  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    643 REFLECTION INTENSITIES INTEGRATED


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
    81   0  0.969  5205003    0    134      80     2  0.04532  0.14167
    82   0  0.973  5205351    0    125      72     0  0.04448  0.14322
    83   0  0.971  5204733    0    126      78     1  0.04271  0.10523
    84   0  0.970  5203496    0    138      92     1  0.04261  0.10924
    85   0  0.968  5204325    0    138      86     1  0.04619  0.12687
    86   0  0.966  5204376    0    142      81     6  0.04782  0.09546
    87   0  0.966  5204333    0    131      84     2  0.04623  0.12296
    88   0  0.962  5204930    0    134      73     0  0.04563  0.11972
    89   0  0.965  5204324    0    131      85     1  0.04702  0.12997
    90   0  0.963  5203638    0    131      93     4  0.04411  0.10369
    91   0  0.964  5203629    0    134      94     6  0.04392  0.09031
    92   0  0.961  5205475    0    132      72     4  0.04814  0.09183
    93   0  0.962  5204888    0    127      79     8  0.04453  0.08813
    94   0  0.961  5204867    0    127      75     1  0.04555  0.11542
    95   0  0.961  5204336    0    126      81     1  0.04805  0.12310
    96   0  0.963  5204283    0    132      87     2  0.04783  0.11390
    97   0  0.959  5204319    0    136      83     4  0.04594  0.10258
    98   0  0.961  5203935    0    135      86     4  0.04418  0.09048
    99   0  0.958  5203653    0    136      93     2  0.04523  0.09871
   100   0  0.960  5203227    0    147     101     1  0.04494  0.12021

   464 OUT OF    465 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.54
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.11
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.420    32.140    34.517  88.659 108.471 111.904
 REC. CELL PARAMETERS   0.041673  0.033731  0.030723  83.817  70.600  67.319
 COORDINATES OF UNIT CELL A-AXIS   -14.364   -14.230    18.521
 COORDINATES OF UNIT CELL B-AXIS    14.410    25.579    13.079
 COORDINATES OF UNIT CELL C-AXIS   -21.500    21.488   -16.354
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.004     0.004     0.015
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.038
 LAB COORDINATES OF ROTATION AXIS -0.999991 -0.004096  0.001370
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000613  0.005112  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.24   1243.99
 DETECTOR ORIGIN (PIXELS) AT                     1265.44   1237.28
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.17
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     464 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     31     57     69     58     29     48     62     52
 SIGMAB (degree)       0.034  0.035  0.035  0.034  0.033  0.034  0.035  0.035  0.034
 SIGMAR (degree)       0.030  0.027  0.028  0.031  0.029  0.026  0.028  0.030  0.029


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    464 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  1  2  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  1  4  9 12 10  5  2  1  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  1  6 14 19 15  9  4  1  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  0  3  7 11  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4  6  7  4  2  1  0  0  0     0  0  0  1  6 10  9  6  2  1  0  0  0     0  0  0  1  4  6  5  3  1  1  0  0  0
     0  0  0  4 15 33 41 31 15  6  2  0  0     0  0  1  6 25 54 60 40 18  7  2  0  0     0  0  0  5 17 33 32 21  9  4  1  0  0
     0  0  0  5 22 53 67 51 28 11  4  1  0     0  0  1 10 44 91100 69 35 14  5  1  0     0  0  1  8 29 55 54 35 17  7  2  1  0
     0  0  0  2 11 28 39 31 15  6  2  1  0     0  0  1  5 23 51 59 42 19  8  3  1  0     0  0  1  4 15 30 33 21  9  4  1  0  0
     0  0  0  0  2  5  7  6  3  1  0  0  0     0  0  0  1  4  9 11  8  4  2  0  0  0     0  0  0  1  3  5  6  4  2  1  0  0  0
     0  0  0  0  0  0  1  1  1  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  7  4  2  1  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 12 11  7  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  6  7  4  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    643 REFLECTION INTENSITIES INTEGRATED


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
   101   0  0.961  5203376    0    143      90     2  0.04686  0.11724
   102   0  0.959  5203680    0    144      90     2  0.04894  0.11940
   103   0  0.959  5203983    0    131      86     0  0.04463  0.12633
   104   0  0.958  5205157    0    116      72     0  0.04624  0.13816
   105   0  0.959  5205713    0    105      65     7  0.04370  0.07289
   106   0  0.954  5205311    0    115      74     1  0.04754  0.11645
   107   0  0.958  5204377    0    127      84     4  0.04517  0.10032
   108   0  0.959  5204299    0    135      89     0  0.04362  0.11634
   109   0  0.957  5204457    0    140      88     1  0.04534  0.11908
   110   0  0.958  5204001    0    126      93     5  0.04813  0.10919
   111   0  0.962  5204644    0    137      81     6  0.04515  0.08907
   112   0  0.958  5203978    0    135      91     5  0.04545  0.09888
   113   0  0.955  5203044    0    140      95     3  0.04399  0.09275
   114   0  0.955  5204166    0    126      84     6  0.04813  0.08495
   115   0  0.954  5205065    0    125      79     1  0.04719  0.13552
   116   0  0.954  5205052    0    117      79     2  0.05154  0.13452
   117   0  0.955  5204493    0    126      86     3  0.04739  0.11383
   118   0  0.954  5204121    0    139      91     1  0.04608  0.12680
   119   0  0.954  5203348    0    148      94     6  0.04187  0.07672
   120   0  0.951  5202959    0    151     106     5  0.04374  0.08927

   475 OUT OF    477 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.54
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.414    32.125    34.502  88.640 108.470 111.898
 REC. CELL PARAMETERS   0.041679  0.033744  0.030735  83.841  70.610  67.332
 COORDINATES OF UNIT CELL A-AXIS   -14.357   -14.225    18.521
 COORDINATES OF UNIT CELL B-AXIS    14.403    25.569    13.069
 COORDINATES OF UNIT CELL C-AXIS   -21.492    21.486   -16.337
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.010    -0.004     0.011
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.038
 LAB COORDINATES OF ROTATION AXIS -0.999991 -0.004096  0.001370
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000804  0.005029  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.30   1244.08
 DETECTOR ORIGIN (PIXELS) AT                     1265.24   1237.48
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.06
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     475 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   48     32     57     65     50     26     54     70     55
 SIGMAB (degree)       0.034  0.035  0.034  0.034  0.033  0.035  0.035  0.035  0.035
 SIGMAR (degree)       0.030  0.028  0.029  0.031  0.029  0.027  0.028  0.029  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    475 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  1  1  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  3  7  9  7  3  2  1  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  1  5 11 16 13  7  3  1  0  0
     1  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  2  5  8  8  4  2  1  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0 -1  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  2  4  6  7  3  2  1  0  0  0     0  0  0  2  6 11 10  5  2  1  0  0  0     0  0  0  1  4  7  6  3  1  0  0  0  0
     0  0  1  4 13 28 34 26 11  5  2  0  0     0  0  1  7 28 52 55 38 16  6  2  0  0     0  0  1  5 21 37 36 23 10  4  2  0  0
     0  0  1  5 21 47 60 48 24  9  3  1  0     0  0  1 10 44 90100 72 35 13  4  1  0     0  0  1  9 34 64 65 45 21  8  3  1  0
     0  0  0  2 10 23 33 26 13  5  2  1  0     0  0  1  5 23 46 56 39 19  8  3  1  0     0  0  1  5 18 35 38 25 12  5  1  0  0
     0  0  0  0  1  3  5  4  3  1  0  0  0     0  0  0  1  3  7  9  7  4  1  0  0  0     0  0  0  1  3  6  7  5  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  6  9  9  6  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  9 17 16 10  5  2  1  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  9  9  6  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    645 REFLECTION INTENSITIES INTEGRATED


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
   121   0  0.953  5202826    0    145     104     3  0.04477  0.10851
   122   0  0.954  5203092    0    142     101     0  0.04518  0.11802
   123   0  0.951  5204730    0    139      85     0  0.04749  0.12329
   124   0  0.952  5203997    0    136      92     0  0.04387  0.12487
   125   0  0.952  5204287    0    133      89     2  0.04519  0.10493
   126   0  0.951  5205031    0    125      85     2  0.04653  0.11306
   127   0  0.952  5205272    0    122      84     1  0.04740  0.12084
   128   0  0.950  5204506    0    133      84     0  0.04377  0.13111
   129   0  0.948  5204558    0    133      84     0  0.04524  0.12447
   130   0  0.949  5204856    0    147      83     5  0.04478  0.09853
   131   0  0.947  5203496    0    149      97     7  0.04030  0.07981
   132   0  0.946  5203516    0    147      92     3  0.04074  0.08198
   133   0  0.945  5204113    0    144      87     3  0.04608  0.08248
   134   0  0.946  5204156    0    145      85     4  0.04533  0.08010
   135   0  0.946  5203598    0    144      88     8  0.04342  0.06836
   136   0  0.950  5203729    0    152      95     3  0.04665  0.08880
   137   0  0.948  5204550    0    141      81     3  0.04548  0.08095
   138   0  0.949  5205138    0    126      75     2  0.04358  0.08786
   139   0  0.949  5205013    0    122      75     3  0.04623  0.10097
   140   0  0.949  5204028    0    128      93     1  0.04800  0.13941

   483 OUT OF    486 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.56
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.431    32.125    34.511  88.631 108.477 111.881
 REC. CELL PARAMETERS   0.041648  0.033739  0.030728  83.855  70.609  67.352
 COORDINATES OF UNIT CELL A-AXIS   -14.360   -14.223    18.546
 COORDINATES OF UNIT CELL B-AXIS    14.406    25.574    13.058
 COORDINATES OF UNIT CELL C-AXIS   -21.497    21.489   -16.342
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.006    -0.008    -0.007
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.038
 LAB COORDINATES OF ROTATION AXIS -0.999991 -0.004096  0.001370
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.001051  0.005183  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.40   1244.12
 DETECTOR ORIGIN (PIXELS) AT                     1265.02   1237.31
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.14
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     483 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   46     36     57     66     57     32     53     68     59
 SIGMAB (degree)       0.035  0.035  0.035  0.034  0.033  0.034  0.035  0.035  0.035
 SIGMAR (degree)       0.030  0.028  0.030  0.030  0.028  0.025  0.027  0.030  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    483 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  1  1  2  2  1  1  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  4  7 10  7  4  2  0  0  0
     0  0  0  0  0  0  1  1  1  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  1  5 11 16 15  8  3  1  0  0
     0  0  0  0  0  1  1  0  0  0  0  0 -1     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  3  6  9  8  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4  6  6  4  2  1  0  0  0     0  0  1  2  6  9  9  5  2  1  0  0  0     0  0  0  1  4  6  6  3  1  0  0  0  0
     0  0  1  4 14 29 36 25 12  5  2  0  0     0  0  1  7 27 51 56 36 16  7  2  0  0     0  0  1  5 19 35 36 22 10  4  1  0  0
     0  0  1  5 21 48 61 50 25 10  4  1  0     0  0  1 10 44 89100 73 36 14  5  2  0     0  0  1  8 34 63 66 45 21  8  3  1  0
     0  0  0  2 12 25 35 29 15  6  2  1  0     0  0  1  5 25 48 58 43 21  9  3  1  0     0  0  1  5 19 36 39 27 13  5  2  0  0
     0  0  0  0  2  4  7  5  2  1  0  0  0     0  0  0  1  4  9 12  8  4  1  1  0  0     0  0  0  1  4  7  8  5  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  6 10  9  6  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  2  9 16 16 11  5  2  1  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5 10 10  6  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    641 REFLECTION INTENSITIES INTEGRATED


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
   141   0  0.949  5204291    0    134      91     0  0.04780  0.17368
   142   0  0.948  5203906    0    136      88     3  0.04503  0.13574
   143   0  0.947  5203406    0    143      89     0  0.04649  0.14593
   144   0  0.945  5202380    0    137     104     8  0.04440  0.09151
   145   0  0.945  5203399    0    139      89     7  0.04666  0.09530
   146   0  0.945  5204487    0    126      77     0  0.05002  0.16253
   147   0  0.944  5204041    0    122      88     1  0.05109  0.16095
   148   0  0.944  5203892    0    127      88     0  0.04849  0.15457
   149   0  0.943  5202576    0    134      98     1  0.04638  0.13685
   150   0  0.943  5203042    0    139      91     6  0.04649  0.09190
   151   0  0.942  5202962    0    131      93     3  0.04711  0.09322
   152   0  0.943  5203980    0    131      79     1  0.04882  0.10282
   153   0  0.944  5203537    0    128      86     1  0.04832  0.11987
   154   0  0.941  5202881    0    125      94     3  0.04737  0.10033
   155   0  0.942  5202185    0    135     104     7  0.04718  0.10186
   156   0  0.945  5202852    0    138      96     9  0.04601  0.08411
   157   0  0.944  5203726    0    143      85     5  0.04690  0.10424
   158   0  0.941  5203341    0    139      82     2  0.04711  0.10350
   159   0  0.938  5203400    0    139      89     0  0.04893  0.13208
   160   0  0.942  5203411    0    133      89     1  0.04645  0.11545

   500 OUT OF    502 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.54
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.11
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.421    32.128    34.508  88.632 108.471 111.903
 REC. CELL PARAMETERS   0.041669  0.033742  0.030730  83.847  70.612  67.329
 COORDINATES OF UNIT CELL A-AXIS   -14.361   -14.225    18.529
 COORDINATES OF UNIT CELL B-AXIS    14.406    25.574    13.063
 COORDINATES OF UNIT CELL C-AXIS   -21.494    21.489   -16.341
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.008    -0.002     0.003
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.038
 LAB COORDINATES OF ROTATION AXIS -0.999991 -0.004096  0.001370
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000845  0.005125  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.34   1244.06
 DETECTOR ORIGIN (PIXELS) AT                     1265.23   1237.33
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.11
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     500 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   52     33     62     69     56     29     57     69     56
 SIGMAB (degree)       0.035  0.036  0.035  0.034  0.033  0.035  0.035  0.035  0.035
 SIGMAR (degree)       0.030  0.027  0.029  0.031  0.030  0.027  0.028  0.029  0.029


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    500 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  1  2  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  1  1  1  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  2  5  9 11  9  4  2  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0     0  0  0  2  6 14 19 15  8  4  1  0  0
     0  0  0  0  0  0  0  0  1  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  3  7 11 10  6  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  1  2  4  7  7  4  2  0  0  0  0     0  0  1  2  6 10 10  5  2  1  0  0  0     0  0  0  1  4  6  6  3  1  0  0  0  0
     0  0  1  5 17 33 39 28 13  5  2  0  0     0  0  1  8 29 56 58 38 16  7  2  0  0     0  0  1  6 20 35 34 20  9  4  1  0  0
     0  0  1  6 24 54 69 52 26 11  4  1  0     0  0  2 11 46 91100 71 33 14  5  1  0     0  0  1  9 32 56 57 38 17  7  2  0  0
     0  0  0  3 12 29 39 32 17  7  2  1  0     0  0  1  6 24 50 58 42 21  8  3  1  0     0  0  1  5 16 30 33 22 11  4  1  0  0
     0  0  0  1  2  5  6  6  3  1  0  0  0     0  0  0  1  4  9 11  8  4  2  0  0  0     0  0  0  1  3  6  7  5  2  1  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  1  1  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  8  8  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  3  8 13 12  8  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  7  4  2  1  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    640 REFLECTION INTENSITIES INTEGRATED


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
   161   0  0.940  5204935    0    126      73     1  0.04715  0.13751
   162   0  0.941  5204396    0    134      78     1  0.04549  0.12998
   163   0  0.940  5203972    0    130      78     4  0.04526  0.08698
   164   0  0.941  5203845    0    131      81     4  0.04628  0.09576
   165   0  0.940  5204064    0    125      81     4  0.04622  0.10898
   166   0  0.940  5203620    0    130      92     1  0.04937  0.14803
   167   0  0.940  5203931    0    132      88     2  0.04284  0.12962
   168   0  0.939  5204103    0    134      92     0  0.04550  0.14627
   169   0  0.936  5203600    0    139      96     1  0.04575  0.13567
   170   0  0.935  5203452    0    140      94     2  0.04761  0.12080
   171   0  0.934  5202099    0    142     105     2  0.04892  0.12529
   172   0  0.937  5202530    0    148     103     1  0.04640  0.12656
   173   0  0.936  5202600    0    144     104     3  0.04765  0.11616
   174   0  0.937  5203092    0    129      97     1  0.04296  0.11477
   175   0  0.935  5205612    0    127      66     0  0.04607  0.11895
   176   0  0.934  5205150    0    119      74     1  0.04824  0.11525
   177   0  0.933  5204542    0    118      77     3  0.04698  0.09787
   178   0  0.935  5204003    0    127      90     3  0.04762  0.11869
   179   0  0.937  5204670    0    130      82     0  0.04861  0.14978
   180   0  0.936  5204196    0    129      80     1  0.04808  0.12825

   474 OUT OF    475 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.56
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.16
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.403    32.133    34.513  88.644 108.471 111.921
 REC. CELL PARAMETERS   0.041703  0.033742  0.030726  83.827  70.604  67.307
 COORDINATES OF UNIT CELL A-AXIS   -14.357   -14.227    18.504
 COORDINATES OF UNIT CELL B-AXIS    14.406    25.574    13.075
 COORDINATES OF UNIT CELL C-AXIS   -21.495    21.492   -16.346
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.010     0.005     0.017
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.038
 LAB COORDINATES OF ROTATION AXIS -0.999991 -0.004096  0.001370
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000390  0.005122  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.26   1244.05
 DETECTOR ORIGIN (PIXELS) AT                     1265.75   1237.32
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.08
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     474 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   51     33     55     64     56     28     50     65     58
 SIGMAB (degree)       0.035  0.035  0.035  0.035  0.034  0.035  0.035  0.035  0.035
 SIGMAR (degree)       0.031  0.028  0.028  0.030  0.029  0.027  0.029  0.031  0.029


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    474 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  1  1  1  1  0  0  0  0  0  0     0  0  1  2  3  4  3  2  1  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  1  2  2  3  2  1  0  0  0  0     0  0  1  3  8 13 16 11  5  2  1  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  1  2  3  3  3  2  1  0  0  0     0  0  0  2 10 18 24 20 10  5  2  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  1  4  9 13 12  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  1  2  2  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  8  8  5  3  1  0  0  0     0  0  0  1  5  9  9  6  3  1  0  0  0     0  0  0  1  3  5  4  3  1  0  0  0  0
     0  0  1  5 18 35 44 30 15  6  2  0  0     0  0  1  5 24 48 59 37 18  7  2  1  0     0  0  1  4 16 27 27 15  7  3  1  0  0
     0  0  1  6 26 57 74 59 30 13  5  1  0     0  0  1  9 40 86100 73 35 14  5  1  0     0  0  1  7 27 49 48 31 14  6  2  0  0
     0  0  1  3 13 29 42 35 18  8  3  1  0     0  0  1  5 23 47 57 42 21  9  3  1  0     0  0  0  4 15 27 28 19  9  4  1  0  0
     0  0  0  1  2  5  7  7  4  1  1  0  0     0  0  0  1  4  9 11  9  5  2  1  0  0     0  0  0  1  3  6  7  5  2  1  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0 -1  0  0  0  0  0 -1 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  7  6  3  1  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  8 12 10  6  2  1  0  0  0     0  0  0  0  1  2  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  1
     0  0  0  1  4  7  5  4  1  1  0  0  0     0  0  0  0  0  1  0  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    628 REFLECTION INTENSITIES INTEGRATED


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
   181   0  0.933  5203755    0    132      86     3  0.04856  0.10827
   182   0  0.935  5203145    0    144      96     3  0.04619  0.09260
   183   0  0.936  5204062    0    145      86     2  0.04252  0.09094
   184   0  0.937  5204581    0    140      78     4  0.04278  0.07936
   185   0  0.936  5204109    0    129      85     2  0.04377  0.09791
   186   0  0.934  5204721    0    133      76     2  0.04694  0.09102
   187   0  0.936  5204978    0    133      76     5  0.04490  0.07780
   188   0  0.934  5204139    0    125      93     4  0.04367  0.09315
   189   0  0.935  5204386    0    127      91     3  0.04413  0.10734
   190   0  0.935  5205194    0    125      85     1  0.04136  0.12947
   191   0  0.933  5204527    0    128      84     0  0.04271  0.12496
   192   0  0.930  5203931    0    135      90     6  0.04554  0.09217
   193   0  0.930  5203790    0    131      92     6  0.04458  0.08348
   194   0  0.928  5204045    0    129      82     5  0.04543  0.08632
   195   0  0.929  5204370    0    127      84     0  0.04785  0.13140
   196   0  0.929  5203726    0    132      84     0  0.04656  0.12089
   197   0  0.928  5203524    0    145      83     0  0.04398  0.09085
   198   0  0.930  5202629    0    148      96     2  0.04785  0.08864
   199   0  0.930  5202899    0    153      93     1  0.04336  0.08159
   200   0  0.928  5203678    0    141      87     4  0.04236  0.08497

   490 OUT OF    495 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.52
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.04
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.425    32.137    34.520  88.650 108.486 111.886
 REC. CELL PARAMETERS   0.041664  0.033729  0.030723  83.826  70.589  67.339
 COORDINATES OF UNIT CELL A-AXIS   -14.365   -14.223    18.533
 COORDINATES OF UNIT CELL B-AXIS    14.407    25.580    13.072
 COORDINATES OF UNIT CELL C-AXIS   -21.495    21.490   -16.365
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.000     0.013     0.004
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.038
 LAB COORDINATES OF ROTATION AXIS -0.999991 -0.004096  0.001370
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000552  0.005170  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.27   1244.10
 DETECTOR ORIGIN (PIXELS) AT                     1265.54   1237.31
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.15
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     490 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   46     30     58     68     57     29     60     72     57
 SIGMAB (degree)       0.035  0.035  0.035  0.035  0.034  0.035  0.035  0.036  0.036
 SIGMAR (degree)       0.029  0.027  0.029  0.030  0.029  0.026  0.028  0.030  0.029


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    490 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  1  2  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  1  5  8 10  7  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  2  6 12 16 14  8  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  3  6  9  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  1  1  0  0  0  0  0  0     0  0  0  1  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  2  4  6  7  4  2  1  0  0  0     0  0  0  2  6 10 11  6  3  1  0  0  0     0  0  0  1  4  6  7  4  2  0  0  0  0
     0  0  1  4 16 29 37 26 14  5  2  0  0     0  0  1  7 28 51 60 39 19  7  3  1  0     0  0  1  6 20 34 36 22 10  4  1  0  0
     0  0  1  5 23 47 59 49 26 11  4  1  0     0  0  1 10 44 89100 74 37 15  5  1  0     0  0  1  8 32 59 61 42 21  8  3  1  0
     0  0  0  3 12 25 33 29 16  7  2  1  0     0  0  1  5 23 48 57 44 22 10  3  1  0     0  0  1  4 17 32 36 25 12  5  2  0  0
     0  0  0  0  2  4  7  6  4  1  0  0  0     0  0  0  1  4  9 12  9  5  2  1  0  0     0  0  0  1  3  6  7  5  3  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  2  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  8  8  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  3  8 14 14  9  4  2  1  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  8  9  6  3  1  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    634 REFLECTION INTENSITIES INTEGRATED


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
   201   0  0.928  5204024    0    136      80     2  0.04436  0.09668
   202   0  0.928  5204389    0    124      79     0  0.04647  0.11628
   203   0  0.926  5204583    0    122      73     1  0.04600  0.10148
   204   0  0.929  5205011    0    124      69     6  0.05047  0.09213
   205   0  0.929  5203826    0    132      90     6  0.04965  0.10647
   206   0  0.928  5203365    0    136      94     6  0.04543  0.09810
   207   0  0.924  5203530    0    127      89     6  0.04744  0.10021
   208   0  0.925  5204116    0    136      89     2  0.04924  0.14704
   209   0  0.927  5204528    0    119      82     0  0.05063  0.15423
   210   0  0.928  5204012    0    119      82     0  0.04942  0.14004
   211   0  0.927  5203933    0    130      84     2  0.05182  0.13612
   212   0  0.924  5202887    0    131      96     4  0.04693  0.10562
   213   0  0.924  5202565    0    130      91     3  0.04703  0.09874
   214   0  0.926  5203075    0    129      94     4  0.05182  0.11705
   215   0  0.923  5203690    0    132      86     0  0.04927  0.15954
   216   0  0.925  5204106    0    118      82     0  0.04999  0.15098
   217   0  0.925  5204236    0    118      82     1  0.05295  0.14842
   218   0  0.923  5204652    0    124      78     1  0.04984  0.12890
   219   0  0.919  5202964    0    133      96     4  0.04531  0.09276
   220   0  0.918  5203205    0    147      88     3  0.04628  0.09058

   471 OUT OF    472 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.64
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.16
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.418    32.142    34.513  88.635 108.469 111.896
 REC. CELL PARAMETERS   0.041670  0.033725  0.030725  83.848  70.614  67.336
 COORDINATES OF UNIT CELL A-AXIS   -14.362   -14.223    18.526
 COORDINATES OF UNIT CELL B-AXIS    14.410    25.585    13.073
 COORDINATES OF UNIT CELL C-AXIS   -21.496    21.492   -16.346
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.007     0.004     0.000
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.038
 LAB COORDINATES OF ROTATION AXIS -0.999991 -0.004096  0.001370
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000745  0.005048  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.29   1244.03
 DETECTOR ORIGIN (PIXELS) AT                     1265.31   1237.40
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.14
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     471 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   50     30     59     67     54     28     48     65     56
 SIGMAB (degree)       0.035  0.036  0.035  0.034  0.034  0.035  0.035  0.036  0.036
 SIGMAR (degree)       0.032  0.028  0.029  0.031  0.029  0.026  0.028  0.030  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    471 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0 -1  0
    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  1  2  2  3  2  1  0  0  0  0
     0  0  0  0  1  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  2  6 10 13 10  5  2  1  0  0
    -1  0  0  0  0  0  0  1  1  0  0  0  0     0  0  0  0  1  2  3  3  2  1  0  0  0     0  0  0  2  7 16 21 20 11  5  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  1  0  0  0     0  0  0  1  4  8 12 10  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  1  2  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  1  1  0  0  0  0  0  0     0  0  0  1  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  3  5  7  7  5  2  1  0  0  0     0  0  0  3  6  9  9  6  2  1  0  0  0     0  0  0  1  3  5  5  3  1  0  0  0  0
     0  0  1  7 21 37 42 31 14  6  2  0  0     0  0  2  9 32 56 56 37 17  7  2  0  0     0  0  1  6 20 31 30 19  8  3  1  0  0
     0  0  1  7 28 58 72 58 30 13  4  1  0     0  0  2 12 48 91100 75 37 15  5  1  0     0  0  1  8 29 52 52 37 18  7  2  1  0
     0  0  0  3 14 31 41 31 17  8  3  1  0     0  0  1  5 25 51 60 44 22  9  3  1  0     0  0  1  4 15 28 32 22 11  4  1  0  0
     0  0  0  1  2  5  6  5  3  1  0  0  0     0  0  0  1  4  9 10  8  4  2  0  0  0     0  0  0  1  3  5  5  4  2  1  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  1  1  2  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0 -1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  2  5  7  7  5  2  1  0  0  0     0  0  0  1  1  1  1  1  0  0  0  0  0     0  0  0  0  1  0  1  0  0  0  0  0  0
     0  0  1  2  7 12 11  8  3  2  0  0  0     0  0  0  1  1  1  1  1  0  0  0  0 -1     0  0  0  0  1  1  0  0  0  0  0  0 -1
     0  0  0  1  4  6  7  4  2  1  0  0  0     0  0  0  1  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    640 REFLECTION INTENSITIES INTEGRATED


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
   221   0  0.920  5202553    0    153      99     0  0.04516  0.09973
   222   0  0.922  5202503    0    152     104     1  0.04478  0.09618
   223   0  0.920  5203285    0    149      92     2  0.04788  0.11181
   224   0  0.924  5203938    0    143      88     2  0.04883  0.11308
   225   0  0.921  5203215    0    131      91     3  0.04839  0.11792
   226   0  0.922  5202979    0    141      97     2  0.04928  0.12338
   227   0  0.923  5203085    0    143      99     1  0.04608  0.13208
   228   0  0.922  5202821    0    147      98     4  0.04585  0.09518
   229   0  0.923  5202984    0    152     100     7  0.04741  0.08629
   230   0  0.921  5202726    0    144      97     5  0.04464  0.09006
   231   0  0.925  5202782    0    146      97     1  0.04673  0.09335
   232   0  0.924  5203360    0    145      90     2  0.04624  0.09459
   233   0  0.921  5203662    0    155      88     2  0.04679  0.09656
   234   0  0.921  5202819    0    143      94     3  0.04762  0.09780
   235   0  0.922  5203467    0    137      92     4  0.04538  0.09222
   236   0  0.922  5203958    0    128      86     4  0.04665  0.10482
   237   0  0.922  5203572    0    130      94     5  0.04672  0.10991
   238   0  0.921  5204440    0    134      88     4  0.04420  0.11162
   239   0  0.920  5203749    0    129      90     4  0.04573  0.10720
   240   0  0.920  5203650    0    127      95     5  0.04764  0.11905

   493 OUT OF    495 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.56
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.11
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.420    32.135    34.509  88.625 108.464 111.904
 REC. CELL PARAMETERS   0.041668  0.033733  0.030727  83.859  70.622  67.332
 COORDINATES OF UNIT CELL A-AXIS   -14.362   -14.220    18.530
 COORDINATES OF UNIT CELL B-AXIS    14.409    25.583    13.061
 COORDINATES OF UNIT CELL C-AXIS   -21.496    21.490   -16.341
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.005     0.002    -0.007
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.038
 LAB COORDINATES OF ROTATION AXIS -0.999991 -0.004096  0.001370
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000745  0.005145  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.39   1244.05
 DETECTOR ORIGIN (PIXELS) AT                     1265.41   1237.30
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.12
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     493 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   52     35     56     65     58     34     56     64     59
 SIGMAB (degree)       0.035  0.036  0.036  0.035  0.034  0.035  0.035  0.036  0.036
 SIGMAR (degree)       0.030  0.027  0.029  0.031  0.029  0.025  0.027  0.031  0.029


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    493 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  1  2  3  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  1     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  2  6 11 12  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  3  2  1  0  0  0     0  0  0  2  7 15 21 17 10  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  1  4  8 13 12  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  1  1  1  0  0  0  0  0     0  0  0  1  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  1  3  5  8  7  5  3  1  0  0  0     0  0  1  3  6 11  9  6  3  1  0  0  0     0  0  0  2  4  6  5  3  2  0  0  0  0
     0  0  1  5 18 36 40 30 14  5  2  0  0     0  0  1  8 29 55 56 40 18  7  2  1  0     0  0  1  6 18 31 29 20  9  3  1  0  0
     0  0  1  7 27 55 70 53 29 12  5  1  0     0  0  1 11 46 89100 70 36 15  5  2  0     0  0  1  8 28 50 53 36 18  7  2  1  0
     0  0  1  3 14 29 43 36 20  8  3  1  0     0  0  1  6 25 50 63 47 24  9  3  1  0     0  0  1  4 14 28 33 23 11  4  1  0  0
     0  0  0  1  2  5  7  6  4  2  1  0  0     0  0  0  1  4 10 12  9  5  2  1  0  0     0  0  0  1  3  6  6  4  2  1  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
    -1  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  1  1  1  0  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  2  4  6  6  4  1  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0    -1  0  0  1  0  1  0  0  0  0  0  0  0
     0  0  0  2  6 10  9  6  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  1  1  3  6  6  4  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  1  1  1  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0 -1  0  0  0  0  0  0  0  0  0  0  0

    632 REFLECTION INTENSITIES INTEGRATED


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
   241   0  0.920  5204486    0    128      83     3  0.05362  0.14119
   242   0  0.924  5204102    0    125      90     3  0.04860  0.13012
   243   0  0.922  5204546    0    136      84     2  0.04662  0.13916
   244   0  0.922  5203765    0    137      91     2  0.05123  0.13775
   245   0  0.921  5202174    0    146     105     2  0.05064  0.13716
   246   0  0.919  5201799    0    170     106     8  0.04776  0.08754
   247   0  0.918  5201238    0    172     110     6  0.04492  0.08063
   248   0  0.918  5201828    0    175     105     4  0.04514  0.08750
   249   0  0.917  5201244    0    167     108     5  0.04672  0.09709
   250   0  0.919  5201759    0    164     106     2  0.04584  0.10632
   251   0  0.918  5202530    0    156     101     2  0.04342  0.11792
   252   0  0.919  5202726    0    137     100     2  0.04580  0.12208
   253   0  0.917  5203066    0    130      95     8  0.04380  0.11010
   254   0  0.916  5204027    0    131      89     3  0.04340  0.11965
   255   0  0.917  5203877    0    139      94     5  0.04490  0.11764
   256   0  0.919  5203030    0    141      98     4  0.04897  0.11688
   257   0  0.918  5204072    0    130      89     3  0.04919  0.13299
   258   0  0.917  5205035    0    121      81     1  0.05184  0.20162
   259   0  0.916  5204754    0    107      84     1  0.05481  0.19332
   260   0  0.919  5203877    0    115      91     3  0.04988  0.12719

   503 OUT OF    505 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.57
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.14
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.408    32.136    34.507  88.624 108.457 111.912
 REC. CELL PARAMETERS   0.041686  0.033734  0.030728  83.859  70.629  67.324
 COORDINATES OF UNIT CELL A-AXIS   -14.358   -14.222    18.514
 COORDINATES OF UNIT CELL B-AXIS    14.407    25.580    13.069
 COORDINATES OF UNIT CELL C-AXIS   -21.495    21.492   -16.333
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.011     0.001     0.004
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.038
 LAB COORDINATES OF ROTATION AXIS -0.999991 -0.004096  0.001370
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000678  0.005031  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.31   1244.07
 DETECTOR ORIGIN (PIXELS) AT                     1265.42   1237.47
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.08
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     503 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   50     35     66     71     57     28     57     68     55
 SIGMAB (degree)       0.035  0.036  0.035  0.035  0.034  0.034  0.035  0.036  0.035
 SIGMAR (degree)       0.031  0.027  0.029  0.031  0.030  0.027  0.028  0.031  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    503 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  1  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  1  2  3  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  1  1  0  0  0  0     0  0  0  1  1  2  2  1  1  1  0  0  0     0  0  1  3  8 12 15 10  6  3  1  0  0
     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  1  3  3  3  2  1  0  0  0     0  0  1  3  9 19 26 21 13  5  2  1  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  1  2  2  1  1  0  0  0     0  0  0  1  4  9 13 13  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  1  2  2  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  7  8  5  2  1  0  0  0     0  0  0  2  5  8 10  5  2  1  0  0  0     0  0  0  1  3  4  4  2  1  0  0  0  0
     0  0  1  6 19 34 43 31 16  7  2  1  0     0  0  1  8 27 47 53 35 17  7  2  0  0     0  0  1  5 15 24 25 16  7  3  1  0  0
     0  0  1  7 28 59 79 62 34 14  5  1  0     0  0  2 13 47 87100 73 36 14  5  1  0     0  0  1  8 29 47 48 34 16  6  2  1  1
     0  0  0  3 13 30 42 36 19  9  3  1  0     0  0  1  6 24 48 58 42 21  9  3  1  0     0  0  1  5 16 27 29 20  9  4  1  0  0
     0  0  0  0  2  4  6  6  4  1  0  0  0     0  0  0  1  4  8  9  7  4  1  0  0  0     0  0  0  1  3  5  5  4  2  1  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0 -1  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  5  5  3  1  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  2  7 10  9  6  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  1  4  6  6  4  2  1  0  0  0     0  0  0  0  1  0  1  0  0  0  0  0  0     0  0  0  0  1  0  1  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    649 REFLECTION INTENSITIES INTEGRATED


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
   261   0  0.918  5203421    0    126      96     4  0.04923  0.11240
   262   0  0.917  5203911    0    138      93     3  0.04677  0.11294
   263   0  0.917  5203519    0    134      89     2  0.04573  0.09487
   264   0  0.916  5203809    0    138      90     2  0.04685  0.11375
   265   0  0.917  5203472    0    134      90     8  0.04422  0.08733
   266   0  0.918  5203897    0    132      86     9  0.04417  0.09057
   267   0  0.919  5203684    0    135      97    10  0.04719  0.12843
   268   0  0.917  5204025    0    136      91     7  0.04705  0.13409
   269   0  0.918  5204259    0    127      84     7  0.04485  0.10397
   270   0  0.914  5205734    0    122      72     5  0.04655  0.12445
   271   0  0.915  5205505    0    121      74     4  0.04514  0.10668
   272   0  0.913  5204603    0    123      81     5  0.04295  0.08635
   273   0  0.912  5203721    0    140      93     7  0.04514  0.08653
   274   0  0.913  5203417    0    133      97     4  0.04844  0.10243
   275   0  0.912  5204806    0    136      82     4  0.04743  0.11343
   276   0  0.912  5204745    0    123      80     7  0.04722  0.10461
   277   0  0.914  5205215    0    120      73     6  0.04861  0.07946
   278   0  0.913  5205752    0    120      73     4  0.04705  0.13560
   279   0  0.914  5205424    0    123      74     8  0.04677  0.09537
   280   0  0.912  5203878    0    124      88     4  0.04334  0.11243

   473 OUT OF    474 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.57
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.15
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.414    32.136    34.510  88.628 108.463 111.903
 REC. CELL PARAMETERS   0.041676  0.033732  0.030726  83.856  70.622  67.332
 COORDINATES OF UNIT CELL A-AXIS   -14.360   -14.221    18.523
 COORDINATES OF UNIT CELL B-AXIS    14.409    25.581    13.067
 COORDINATES OF UNIT CELL C-AXIS   -21.495    21.492   -16.339
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.009     0.002    -0.001
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.038
 LAB COORDINATES OF ROTATION AXIS -0.999991 -0.004096  0.001370
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000668  0.005112  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.37   1244.05
 DETECTOR ORIGIN (PIXELS) AT                     1265.50   1237.34
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.09
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     473 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   48     31     60     64     57     30     46     62     55
 SIGMAB (degree)       0.035  0.036  0.036  0.035  0.034  0.035  0.035  0.036  0.035
 SIGMAR (degree)       0.031  0.026  0.027  0.030  0.028  0.027  0.029  0.030  0.029


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    473 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  1  2  2  2  2  1  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  1  2  7 11 12  9  5  2  1  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  2  1  0  0  0     0  0  0  2  8 17 21 18 11  4  2  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  1  3  9 13 12  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4  7  7  5  2  1  0  0  0     0  0  0  2  5 10  9  6  3  1  0  0  0     0  0  0  1  3  5  5  3  1  0  0  0  0
     0  0  1  6 19 35 41 30 15  6  2  0  0     0  0  1  8 27 53 58 41 18  8  3  1  0     0  0  1  5 17 30 29 20  8  3  1  0  0
     0  0  1  7 27 55 70 57 31 13  5  1  0     0  0  1 11 45 87100 74 38 16  6  2  0     0  0  1  8 29 51 51 36 18  7  2  1  0
     0  0  0  3 12 29 41 36 19  8  3  1  0     0  0  1  5 22 48 60 46 23  9  3  1  0     0  0  1  4 15 27 31 22 10  4  1  0  0
     0  0  0  0  2  4  7  6  3  2  0  0  0     0  0  0  1  4  7  9  8  4  2  0  0  0     0  0  0  1  3  5  5  4  2  1  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  2  5  7  6  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  3  7 11  9  6  3  1  0  0  0     0  0  0  1  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  4  6  6  4  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  1  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    610 REFLECTION INTENSITIES INTEGRATED


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
   281   0  0.911  5203989    0    137      90     2  0.04284  0.12556
   282   0  0.911  5204014    0    150      91    10  0.04282  0.08386
   283   0  0.910  5204400    0    140      88     8  0.04324  0.07689
   284   0  0.910  5204562    0    133      85     9  0.04516  0.10618
   285   0  0.911  5203854    0    133      90     8  0.04806  0.10639
   286   0  0.914  5204559    0    133      86     6  0.04660  0.12967
   287   0  0.915  5205239    0    129      80     7  0.04529  0.11160
   288   0  0.916  5205230    0    136      83     8  0.04246  0.12232
   289   0  0.915  5204854    0    145      81     5  0.04472  0.10280
   290   0  0.911  5203480    0    150      93     5  0.04393  0.08774
   291   0  0.914  5203529    0    157      94     5  0.04396  0.08705
   292   0  0.915  5202688    0    163      96     4  0.04524  0.07387
   293   0  0.916  5202709    0    156      98     4  0.04657  0.08350
   294   0  0.916  5203588    0    155      88     3  0.04458  0.08728
   295   0  0.917  5203378    0    142      88     2  0.04701  0.09523
   296   0  0.919  5203069    0    134      93     2  0.04919  0.12062
   297   0  0.917  5204496    0    136      79     1  0.04616  0.15459
   298   0  0.917  5204746    0    145      79     2  0.05368  0.15313
   299   0  0.918  5203489    0    145      92     2  0.04944  0.16179
   300   0  0.917  5202400    0    134      98     3  0.04536  0.16331

   515 OUT OF    516 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.56
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.11
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.419    32.126    34.500  88.628 108.466 111.892
 REC. CELL PARAMETERS   0.041666  0.033740  0.030736  83.858  70.620  67.343
 COORDINATES OF UNIT CELL A-AXIS   -14.358   -14.216    18.535
 COORDINATES OF UNIT CELL B-AXIS    14.406    25.576    13.056
 COORDINATES OF UNIT CELL C-AXIS   -21.491    21.482   -16.336
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.004    -0.004    -0.011
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.038
 LAB COORDINATES OF ROTATION AXIS -0.999991 -0.004096  0.001370
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000924  0.005189  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.37   1244.06
 DETECTOR ORIGIN (PIXELS) AT                     1265.15   1237.25
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.02
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     515 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   50     29     59     76     56     34     71     70     53
 SIGMAB (degree)       0.036  0.036  0.036  0.035  0.034  0.035  0.035  0.036  0.036
 SIGMAR (degree)       0.030  0.028  0.029  0.030  0.029  0.026  0.027  0.030  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    515 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  1  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  1  3  7  9  7  3  1  1  0  0
     0  0  0  0  0  0  1  0  1  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  1  5 11 15 13  8  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  2  6  9  8  5  2  1  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0
     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  1  1  1  1  0  0  0  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0
     0  0  0  2  4  7  7  4  2  1  0  0  0     0  0  1  3  7 12 12  7  3  1  0  0  0     0  0  0  2  5  8  8  4  2  1  0  0  0
     0  0  1  5 15 30 36 27 13  5  2  0  0     0  0  1  9 30 54 59 41 19  8  3  1  0     0  0  1  7 23 37 38 25 11  4  1  0  0
     0  0  1  6 24 48 62 50 28 11  4  1  0     0  0  2 12 46 86100 73 40 15  6  2  0     0  0  1 11 36 61 63 42 23  9  3  1  0
     0  0  1  3 11 25 34 29 16  7  2  1  1     0  0  1  5 21 45 56 42 22  9  3  1  0     0  0  1  5 18 33 36 25 12  5  2  0  0
     0  0  0  0  2  4  6  5  3  1  0  0  0     0  0  0  1  4  8 10  8  4  2  0  0  0     0  0  0  1  3  6  6  5  2  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  3  2  1  1  0  0  0  0     0  0  0  1  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 10 10  6  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  3 11 17 17 10  5  2  1  0  0     0  0  0  0  1  1  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6  9 10  7  3  1  0  0  0     0  0  0  1  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    667 REFLECTION INTENSITIES INTEGRATED


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
   301   0  0.915  5203071    0    143      93     6  0.05379  0.12952
   302   0  0.915  5202831    0    128     101     6  0.05055  0.13205
   303   0  0.915  5203461    0    132      89     9  0.04681  0.09346
   304   0  0.916  5203536    0    129      93     9  0.04709  0.07299
   305   0  0.918  5203795    0    113      94     7  0.04708  0.10296
   306   0  0.917  5205136    0    121      77     4  0.04928  0.12950
   307   0  0.919  5204809    0    120      82     2  0.04962  0.11968
   308   0  0.920  5204396    0    125      86     4  0.04901  0.11554
   309   0  0.920  5205114    0    134      73     6  0.04696  0.09624
   310   0  0.916  5203903    0    128      89     6  0.04286  0.11810
   311   0  0.917  5204146    0    122      92     4  0.04387  0.11592
   312   0  0.920  5204753    0    125      83     4  0.05050  0.10707
   313   0  0.920  5204618    0    114      82     3  0.04808  0.13793
   314   0  0.918  5205110    0    114      72     4  0.04676  0.15949
   315   0  0.918  5204685    0    117      76     4  0.04874  0.12271
   316   0  0.918  5204242    0    120      84     9  0.04651  0.09474
   317   0  0.919  5204202    0    125      85     4  0.04774  0.12200
   318   0  0.921  5203806    0    129      86     2  0.04866  0.12002
   319   0  0.919  5203565    0    127      91     2  0.04694  0.12596
   320   0  0.917  5202856    0    134      97     1  0.05274  0.14436

   456 OUT OF    456 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.57
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.13
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.415    32.132    34.509  88.640 108.468 111.887
 REC. CELL PARAMETERS   0.041673  0.033733  0.030728  83.846  70.614  67.344
 COORDINATES OF UNIT CELL A-AXIS   -14.359   -14.218    18.527
 COORDINATES OF UNIT CELL B-AXIS    14.407    25.576    13.069
 COORDINATES OF UNIT CELL C-AXIS   -21.495    21.489   -16.343
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.008     0.002    -0.001
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.038
 LAB COORDINATES OF ROTATION AXIS -0.999991 -0.004096  0.001370
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000613  0.005072  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.35   1244.10
 DETECTOR ORIGIN (PIXELS) AT                     1265.55   1237.45
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.07
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     456 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   49     30     53     60     50     28     42     64     58
 SIGMAB (degree)       0.035  0.036  0.036  0.035  0.034  0.035  0.035  0.035  0.036
 SIGMAR (degree)       0.031  0.026  0.028  0.030  0.028  0.025  0.029  0.031  0.029


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    456 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  1  1  2  3  1  1  0  0  0  0
     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  2  5 10 10  8  4  2  1  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  2  3  3  2  1  0  0  0     0  0  0  2  7 14 18 15  8  4  1  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  1  4  7 10  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  2  5  8  9  5  2  1  0  0  0     0  0  0  3  6 12 13  8  3  1  0  0  0     0  0  0  2  4  7  7  4  1  1  0  0  0
     0  0  1  6 17 33 38 28 14  6  2  0  0     0  0  1  8 27 54 62 41 21  8  3  1  0     0  0  1  6 20 36 36 22 11  4  2  0  0
     0  0  1  7 24 49 61 50 27 11  4  1  0     0  0  1 12 43 83100 74 38 15  6  2  0     0  0  1 11 34 56 61 42 21  8  3  1  0
     0  0  0  3 13 26 33 28 15  7  2  1  0     0  0  1  6 23 46 56 44 21  9  3  1  0     0  0  1  5 19 31 34 25 12  5  2  0  0
     0  0  0  0  2  4  6  5  3  1  0  0  0     0  0  0  1  4  9 11  7  4  2  1  0  0     0  0  0  1  3  6  7  5  3  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  9  8  5  2  1  0  0  0     0  0  0  1  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  9 14 14  9  4  2  1  0  0     0  0  0  1  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  8  8  6  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    600 REFLECTION INTENSITIES INTEGRATED


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
   321   0  0.917  5203353    0    139      93     0  0.05092  0.16846
   322   0  0.918  5203229    0    133      94     2  0.04736  0.15562
   323   0  0.917  5203447    0    138      94     6  0.04883  0.12518
   324   0  0.917  5203876    0    146      82     5  0.04952  0.09233
   325   0  0.920  5202843    0    137      94     8  0.04567  0.07033
   326   0  0.921  5202978    0    131      95     3  0.04823  0.08827
   327   0  0.919  5203830    0    139      86     2  0.05088  0.11257
   328   0  0.920  5205180    0    130      69     4  0.04822  0.10346
   329   0  0.919  5204985    0    126      74     8  0.04803  0.08794
   330   0  0.920  5205927    0    121      67     1  0.05284  0.14874
   331   0  0.922  5205162    0    117      73     2  0.05116  0.11143
   332   0  0.923  5204635    0    121      80     3  0.04384  0.10379
   333   0  0.920  5204334    0    131      80     3  0.04598  0.11310
   334   0  0.920  5204174    0    147      85     1  0.04613  0.12148
   335   0  0.922  5203294    0    147      97     2  0.04536  0.11142
   336   0  0.923  5203418    0    150      89     3  0.04330  0.08224
   337   0  0.923  5204036    0    143      89     1  0.04818  0.09738
   338   0  0.923  5204827    0    140      80     2  0.04394  0.09614
   339   0  0.924  5204512    0    133      82     2  0.04698  0.09866
   340   0  0.923  5205374    0    135      73     2  0.04541  0.10152

   506 OUT OF    507 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.50
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.11
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.408    32.133    34.503  88.638 108.452 111.890
 REC. CELL PARAMETERS   0.041679  0.033733  0.030730  83.855  70.632  67.344
 COORDINATES OF UNIT CELL A-AXIS   -14.354   -14.214    18.523
 COORDINATES OF UNIT CELL B-AXIS    14.405    25.579    13.066
 COORDINATES OF UNIT CELL C-AXIS   -21.499    21.482   -16.334
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.002    -0.005    -0.007
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.038
 LAB COORDINATES OF ROTATION AXIS -0.999991 -0.004096  0.001370
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000747  0.005041  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.32   1244.05
 DETECTOR ORIGIN (PIXELS) AT                     1265.34   1237.43
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.07
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     506 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   50     33     61     73     56     29     60     71     56
 SIGMAB (degree)       0.035  0.036  0.036  0.035  0.034  0.035  0.035  0.036  0.036
 SIGMAR (degree)       0.031  0.027  0.028  0.030  0.029  0.026  0.028  0.030  0.029


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    506 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  2  5 10 12  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  3  2  1  0  0  0     0  0  0  2  8 15 19 17 10  4  2  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  1  4  8 12 10  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  2  4  7  8  5  2  1  0  0  0     0  0  0  2  6 11 10  6  3  1  0  0  0     0  0  0  1  4  6  5  3  1  0  0  0  0
     0  0  1  6 17 35 41 32 15  7  2  0  0     0  0  1  9 29 56 59 41 19  8  3  1  0     0  0  1  6 17 31 30 19  9  3  1  0  0
     0  0  1  8 28 56 68 56 29 12  5  1  0     0  0  2 14 50 92100 73 36 15  5  1  1     0  0  1 10 32 52 52 35 17  7  2  0  0
     0  0  0  4 15 31 40 34 17  8  3  1  0     0  0  1  8 28 52 59 45 21  9  3  1  0     0  0  1  6 18 30 31 22 10  4  1  0  0
     0  0  0  1  2  5  7  6  4  2  0  0  0     0  0  0  1  4 10 11  8  5  2  0  0  0     0  0  0  1  3  6  6  4  2  1  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  6  4  2  1  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  3  7 11 10  6  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4  7  6  4  2  1  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    654 REFLECTION INTENSITIES INTEGRATED


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
   341   0  0.922  5205647    0    130      70     5  0.04644  0.07088
   342   0  0.924  5205812    0    122      66     5  0.04346  0.06556
   343   0  0.926  5205851    0    127      66     4  0.04540  0.06998
   344   0  0.927  5205770    0    133      69     2  0.04533  0.07820
   345   0  0.926  5205121    0    131      77     2  0.04557  0.08962
   346   0  0.927  5205795    0    136      74     0  0.04611  0.11798
   347   0  0.931  5204916    0    126      79     0  0.04796  0.11524
   348   0  0.931  5204065    0    134      84     0  0.04547  0.09786
   349   0  0.932  5204243    0    140      81     1  0.04498  0.08343
   350   0  0.933  5203907    0    156      89     2  0.04759  0.09366
   351   0  0.934  5203597    0    154      88     7  0.04605  0.06849
   352   0  0.934  5202188    0    157      98     2  0.04329  0.07362
   353   0  0.939  5203766    0    157      85     2  0.04507  0.06795
   354   0  0.941  5204642    0    148      81     4  0.04966  0.07663
   355   0  0.943  5204097    0    135      86     3  0.04609  0.10129
   356   0  0.945  5204448    0    122      83     1  0.04804  0.11305
   357   0  0.945  5205041    0    126      80     0  0.04920  0.12276
   358   0  0.948  5204973    0    137      73     5  0.04753  0.07783
   359   0  0.951  5204305    0    139      79     2  0.04528  0.09316
   360   0  0.950  5203459    0    138      98     4  0.04551  0.10324

   488 OUT OF    488 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.51
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.05
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.416    32.135    34.512  88.633 108.465 111.887
 REC. CELL PARAMETERS   0.041670  0.033730  0.030725  83.855  70.620  67.347
 COORDINATES OF UNIT CELL A-AXIS   -14.355   -14.213    18.536
 COORDINATES OF UNIT CELL B-AXIS    14.405    25.587    13.057
 COORDINATES OF UNIT CELL C-AXIS   -21.502    21.482   -16.349
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.011    -0.001    -0.019
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.038
 LAB COORDINATES OF ROTATION AXIS -0.999991 -0.004096  0.001370
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000817  0.005191  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.30   1244.06
 DETECTOR ORIGIN (PIXELS) AT                     1265.23   1237.24
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.12
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     488 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   46     29     56     69     60     34     59     63     57
 SIGMAB (degree)       0.035  0.036  0.036  0.035  0.034  0.035  0.035  0.036  0.035
 SIGMAR (degree)       0.028  0.026  0.029  0.030  0.027  0.024  0.027  0.030  0.029


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    488 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  4  8 10  8  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  1  0  0  0     0  0  0  2  6 13 17 15  8  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  3  7 10  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  3  6  7  4  2  1  0  0  0     0  0  0  2  5  9  9  6  3  1  0  0  0     0  0  0  1  3  5  5  3  2  1  0  0  0
     0  0  1  5 16 31 37 29 14  6  2  0  0     0  0  1  8 28 53 57 40 19  7  3  0  0     0  0  1  5 17 29 31 21 10  4  1  0  0
     0  0  1  6 26 51 64 51 26 11  4  1  0     0  0  1 13 47 88100 71 34 14  5  1  0     0  0  1  9 29 51 54 36 18  7  2  1  0
     0  0  0  3 13 27 38 29 15  6  2  0  0     0  0  1  7 25 48 57 40 19  8  3  0  0     0  0  1  5 16 28 30 21 10  4  1  0  0
     0  0  0  1  2  5  6  6  3  1  0  0  0     0  0  0  1  5  9  9  9  4  2  0  0  0     0  0  0  1  3  5  5  5  2  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0 -1 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  6  5  4  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 10 10  6  3  1  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  6  6  4  2  1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
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
   361   0  0.953  5204865    0    135      78     2  0.04961  0.12073
   362   0  0.955  5205444    0    131      72     0  0.05042  0.14174
   363   0  0.954  5205125    0    122      74     0  0.05085  0.12008
   364   0  0.955  5205239    0    125      72     1  0.05029  0.10293
   365   0  0.956  5204544    0    135      79     0  0.04759  0.11714
   366   0  0.958  5203589    0    127      86     3  0.04583  0.08638
   367   0  0.958  5202981    0    129      92    10  0.04672  0.06959
   368   0  0.958  5204688    0    128      79     5  0.05025  0.09343
   369   0  0.958  5206492    0    120      62     0  0.05433  0.14285
   370   0  0.960  5205519    0    107      68     1  0.04612  0.11540
   371   0  0.961  5205312    0    116      75     2  0.04886  0.11487
   372   0  0.962  5205645    0    122      66     1  0.04742  0.10933
   373   0  0.963  5206081    0    109      65     1  0.04552  0.09560
   374   0  0.966  5206201    0    108      63     1  0.04765  0.09961
   375   0  0.966  5206678    0    116      62     2  0.04289  0.09616
   376   0  0.968  5206248    0    117      68     7  0.04379  0.07188
   377   0  0.971  5205355    0    134      75     3  0.04376  0.07846
   378   0  0.970  5204208    0    139      82     5  0.04299  0.06915
   379   0  0.969  5204954    0    149      74     1  0.04760  0.10013
   380   0  0.971  5204533    0    145      80     2  0.04744  0.09467

   447 OUT OF    448 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.55
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.413    32.130    34.499  88.638 108.462 111.874
 REC. CELL PARAMETERS   0.041670  0.033732  0.030736  83.856  70.622  67.360
 COORDINATES OF UNIT CELL A-AXIS   -14.352   -14.211    18.534
 COORDINATES OF UNIT CELL B-AXIS    14.403    25.579    13.061
 COORDINATES OF UNIT CELL C-AXIS   -21.496    21.476   -16.337
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.003    -0.007    -0.014
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.038
 LAB COORDINATES OF ROTATION AXIS -0.999991 -0.004096  0.001370
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000830  0.004992  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.36   1244.09
 DETECTOR ORIGIN (PIXELS) AT                     1265.27   1237.54
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.04
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     447 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   45     34     55     55     56     23     43     66     57
 SIGMAB (degree)       0.035  0.036  0.036  0.035  0.034  0.035  0.035  0.036  0.036
 SIGMAR (degree)       0.029  0.026  0.028  0.029  0.028  0.026  0.029  0.030  0.029


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    447 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  3  7 10  9  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  2  5 12 16 14  8  3  1  0  0
     0  0  0  0  0  0  0  0  1  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  3  7  9  8  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  4  5  6  4  2  1  0  0  0     0  0  1  3  7 10 10  6  2  1  0  0  0     0  0  0  2  5  7  7  4  2  1  0  0  0
     0  0  1  4 14 30 37 31 15  6  2  0  0     0  0  2  8 28 56 60 44 20  8  3  1  0     0  0  1  6 21 38 37 26 12  5  1  0  0
     0  0  1  6 23 49 60 49 27 12  4  1  0     0  0  2 13 45 91100 73 37 16  6  2  0     0  0  2 10 33 62 62 42 20  9  3  1  0
     0  0  0  3 11 26 33 27 15  6  2  1  0     0  0  1  7 24 48 56 41 22  8  3  1  0     0  0  1  6 18 31 33 23 12  4  1  0  0
     0  0  0  1  2  5  6  6  3  1  0  0  0     0  0  0  1  5  8 11  9  4  2  1  0  0     0  0  0  1  4  6  7  6  2  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  9  8  6  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  8 14 13  8  4  2  1  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  7  7  5  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    617 REFLECTION INTENSITIES INTEGRATED


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
   381   0  0.970  5203759    0    141      83     2  0.04598  0.09162
   382   0  0.969  5204656    0    147      77     3  0.04773  0.09142
   383   0  0.966  5204111    0    131      81     2  0.04432  0.08442
   384   0  0.965  5204327    0    127      75     2  0.04737  0.07599
   385   0  0.968  5205650    0    125      65     1  0.04938  0.09578
   386   0  0.969  5204884    0    120      77     1  0.04574  0.09534
   387   0  0.973  5204873    0    134      75     0  0.04534  0.10123
   388   0  0.972  5203900    0    132      87     5  0.04651  0.07977
   389   0  0.973  5204360    0    137      82     7  0.04783  0.08264
   390   0  0.972  5204312    0    146      88     6  0.04878  0.08486
   391   0  0.974  5203632    0    132      91     6  0.04811  0.08765
   392   0  0.975  5204638    0    130      78     2  0.04455  0.10897
   393   0  0.974  5205071    0    141      78     6  0.04982  0.09315
   394   0  0.971  5204061    0    135      87     6  0.04677  0.08494
   395   0  0.975  5204646    0    131      75     0  0.04647  0.11812
   396   0  0.976  5205549    0    131      68     0  0.05042  0.10485
   397   0  0.978  5204637    0    125      80     0  0.04945  0.11572
   398   0  0.977  5203784    0    132      84     2  0.04476  0.08985
   399   0  0.979  5204261    0    132      80     2  0.04922  0.09219
   400   0  0.981  5204472    0    128      84     1  0.04640  0.11319

   476 OUT OF    477 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.52
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.05
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.411    32.138    34.512  88.641 108.459 111.884
 REC. CELL PARAMETERS   0.041676  0.033726  0.030724  83.849  70.623  67.348
 COORDINATES OF UNIT CELL A-AXIS   -14.356   -14.208    18.530
 COORDINATES OF UNIT CELL B-AXIS    14.405    25.589    13.061
 COORDINATES OF UNIT CELL C-AXIS   -21.502    21.479   -16.354
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.015     0.007    -0.018
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.038
 LAB COORDINATES OF ROTATION AXIS -0.999991 -0.004096  0.001370
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000587  0.005196  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.30   1244.07
 DETECTOR ORIGIN (PIXELS) AT                     1265.53   1237.25
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.12
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     476 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   50     30     58     66     54     28     56     68     60
 SIGMAB (degree)       0.035  0.036  0.036  0.035  0.034  0.035  0.035  0.035  0.035
 SIGMAR (degree)       0.030  0.026  0.029  0.030  0.028  0.026  0.028  0.030  0.028


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    476 STRONG REFLECTIONS *****


     0  1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  1  4  9 12 10  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  2  7 15 20 18 10  4  2  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  1  4  9 12 11  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  6  8  5  2  1  0  0  0     0  0  0  2  6 10 11  6  3  1  0  0  0     0  0  0  1  3  5  6  3  1  0  0  0  0
     0  0  1  4 16 34 43 32 15  7  2  0  0     0  0  1  7 28 53 60 41 18  8  3  1  0     0  0  1  4 16 27 28 18  8  3  1  0  0
     0  0  1  5 25 55 71 60 30 13  5  1  0     0  0  2 11 45 87100 77 36 15  5  2  0     0  0  1  7 26 46 48 35 15  6  2  1  0
     0  0  0  4 15 31 42 33 18  8  3  1  0     0  0  1  7 26 50 60 43 23  9  3  1  0     0  0  1  5 16 27 30 20 10  4  1  0  0
     0  0  0  1  3  6  8  6  4  1  0  0  0     0  0  0  2  5 10 11  8  4  2  1  0  0     0  0  0  1  3  6  6  4  2  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  1  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  5  5  4  1  1  0  0 -1     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 10  9  6  2  1  0  0 -1     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  6  6  4  2  0  0  0 -1     0  0  0  0  1  1  0  0  0  0  0  0 -1     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0 -1  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    641 REFLECTION INTENSITIES INTEGRATED


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
   401   0  0.978  5204531    0    130      84     1  0.05047  0.12072
   402   0  0.981  5204780    0    126      82     4  0.04475  0.08732
   403   0  0.977  5205849    0    132      70     2  0.04636  0.09466
   404   0  0.980  5205524    0    140      70     1  0.04913  0.10819
   405   0  0.980  5205505    0    136      77     3  0.04461  0.10009
   406   0  0.977  5203636    0    155      88     4  0.04557  0.08811
   407   0  0.979  5203092    0    154      93     5  0.04780  0.08447
   408   0  0.978  5203156    0    157      95     6  0.04457  0.07544
   409   0  0.980  5203779    0    151      89     3  0.04834  0.10046
   410   0  0.979  5203159    0    142      90     2  0.05107  0.09739
   411   0  0.977  5203363    0    147      88     2  0.04652  0.10239
   412   0  0.981  5203470    0    150      86     3  0.04425  0.07914
   413   0  0.978  5203814    0    155      84     1  0.04851  0.09147
   414   0  0.984  5203632    0    147      86     0  0.04701  0.09067
   415   0  0.985  5203871    0    143      87     2  0.04713  0.09034
   416   0  0.986  5204084    0    146      91     5  0.05107  0.09460
   417   0  0.987  5204216    0    140      86     4  0.05044  0.09909
   418   0  0.986  5204158    0    138      86     5  0.04653  0.08240
   419   0  0.986  5203855    0    135      84     7  0.04698  0.06862
   420   0  0.988  5203618    0    133      89     8  0.04758  0.08450

   510 OUT OF    510 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.55
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.418    32.134    34.514  88.639 108.477 111.872
 REC. CELL PARAMETERS   0.041666  0.033728  0.030725  83.849  70.606  67.359
 COORDINATES OF UNIT CELL A-AXIS   -14.353   -14.210    18.543
 COORDINATES OF UNIT CELL B-AXIS    14.400    25.589    13.055
 COORDINATES OF UNIT CELL C-AXIS   -21.502    21.476   -16.360
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.022     0.001    -0.023
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.038
 LAB COORDINATES OF ROTATION AXIS -0.999991 -0.004096  0.001370
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000565  0.005209  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.38   1244.04
 DETECTOR ORIGIN (PIXELS) AT                     1265.63   1237.20
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.11
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     510 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   52     32     60     75     56     35     61     70     58
 SIGMAB (degree)       0.035  0.036  0.036  0.035  0.034  0.035  0.035  0.036  0.036
 SIGMAR (degree)       0.029  0.026  0.029  0.030  0.028  0.025  0.028  0.029  0.028


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    510 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  2  6  9  7  4  2  1  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  1  5 11 14 12  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  3  6  8  8  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0 -1  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  3  6  7  5  2  1  0  0  0     0  0  0  2  6 12 12  7  3  1  0  0  0     0  0  0  1  5  8  7  4  2  0  0  0  0
     0  0  1  3 12 27 35 27 14  6  2  0  0     0  0  1  7 26 53 61 44 21  8  3  1  0     0  0  1  6 20 36 37 25 11  4  1  0  0
     0  0  0  5 20 44 58 47 25 10  4  1  0     0  0  1 11 44 86100 74 38 15  5  1  0     0  0  1 10 35 59 62 41 20  8  3  1  0
     0  0  0  3 12 25 33 27 14  6  2  1  0     0  0  1  6 25 49 58 41 21  9  3  1  0     0  0  1  6 20 34 34 23 12  5  1  0  0
     0  0  0  1  2  5  7  6  4  1  0  0  0     0  0  0  1  5  9 13  9  5  2  1  0  0     0  0  0  1  5  7  8  5  3  1  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  8  8  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  9 14 14  8  4  1  1  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  8  8  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    660 REFLECTION INTENSITIES INTEGRATED


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
   421   0  0.987  5204123    0    133      85     7  0.04727  0.07229
   422   0  0.985  5204781    0    130      78     3  0.04642  0.09616
   423   0  0.989  5204359    0    132      80     0  0.04669  0.12501
   424   0  0.989  5204598    0    145      83     0  0.04860  0.12119
   425   0  0.990  5204854    0    140      76     0  0.04851  0.11422
   426   0  0.991  5204592    0    130      78     3  0.04545  0.07768
   427   0  0.994  5204917    0    125      73     2  0.04624  0.08747
   428   0  0.994  5206048    0    127      64     3  0.04861  0.10091
   429   0  0.988  5205445    0    124      75     1  0.04784  0.11131
   430   0  0.989  5205615    0    123      73     0  0.05112  0.11241
   431   0  0.990  5205292    0    116      71     0  0.04730  0.11110
   432   0  0.990  5205042    0    112      70     1  0.04696  0.10544
   433   0  0.991  5204703    0    120      82     8  0.04552  0.08354
   434   0  0.996  5204959    0    122      78     4  0.04385  0.09176
   435   0  0.994  5205294    0    136      74     1  0.04720  0.10896
   436   0  0.999  5204565    0    134      83     3  0.04560  0.08179
   437   0  0.999  5203695    0    146      91     6  0.04475  0.08940
   438   0  0.997  5203912    0    143      89     7  0.04789  0.07889
   439   0  0.997  5203811    0    136      88     4  0.04620  0.09209
   440   0  0.998  5205041    0    138      77     2  0.04646  0.11800

   483 OUT OF    483 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.56
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.08
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.409    32.135    34.499  88.635 108.446 111.881
 REC. CELL PARAMETERS   0.041673  0.033728  0.030733  83.864  70.640  67.355
 COORDINATES OF UNIT CELL A-AXIS   -14.351   -14.208    18.533
 COORDINATES OF UNIT CELL B-AXIS    14.400    25.589    13.056
 COORDINATES OF UNIT CELL C-AXIS   -21.503    21.468   -16.337
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.020    -0.009    -0.021
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.038
 LAB COORDINATES OF ROTATION AXIS -0.999991 -0.004096  0.001370
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000775  0.005103  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.35   1244.03
 DETECTOR ORIGIN (PIXELS) AT                     1265.34   1237.33
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.09
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     483 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   49     35     62     70     51     25     51     68     60
 SIGMAB (degree)       0.035  0.036  0.035  0.035  0.034  0.035  0.035  0.035  0.035
 SIGMAR (degree)       0.031  0.026  0.028  0.030  0.028  0.026  0.028  0.030  0.029


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    483 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  1  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  1  4  9 12 10  5  2  1  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  2  3  3  2  1  0  0  0     0  0  0  2  7 16 20 17 10  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  1  4  8 13 10  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  3  6  7  5  2  1  0  0  0     0  0  0  2  6 10 11  7  3  1  0  0  0     0  0  0  1  3  5  6  3  1  1  0  0  0
     0  0  1  4 15 33 43 32 16  6  2  1  0     0  0  1  8 29 55 62 43 19  8  3  1  0     0  0  1  5 17 29 31 19  9  3  1  0  0
     0  0  1  6 26 57 69 55 29 12  4  1  0     0  0  1 12 45 89100 72 35 15  5  1  0     0  0  1  8 27 48 49 33 16  6  2  0  0
     0  0  1  4 15 31 42 31 16  7  2  1  0     0  0  1  8 26 48 58 40 20  9  3  1  0     0  0  1  5 16 26 28 18  9  4  1  0  0
     0  0  0  1  3  5  8  6  3  1  0  0  0     0  0  0  2  5  9 11  8  4  2  1  0  0     0  0  0  1  4  6  6  4  2  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0
     0  0  0  1  4  6  6  4  1  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  1  1  0  0  1  0  0  0  0
     0  0  0  2  6 10  9  6  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0
     0  0  0  1  4  5  5  3  1  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  1  0  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0 -1 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    627 REFLECTION INTENSITIES INTEGRATED


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
   441   0  1.000  5204646    0    137      79     2  0.04760  0.10859
   442   0  0.997  5203511    0    133      88     4  0.04751  0.08515
   443   0  0.996  5204286    0    138      83     3  0.04700  0.09163
   444   0  0.998  5204388    0    142      81     5  0.04807  0.07870
   445   0  1.000  5204059    0    123      81     5  0.04546  0.09479
   446   0  1.000  5205160    0    118      75     0  0.05200  0.15792
   447   0  1.001  5206018    0    108      60     4  0.05119  0.09627
   448   0  1.001  5206105    0     96      61     0  0.05216  0.14354
   449   0  1.000  5206894    0    107      54     0  0.05326  0.13856
   450   0  1.002  5206062    0    112      66     3  0.05034  0.09408
   451   0  0.998  5205207    0    117      76     3  0.04245  0.08218
   452   0  1.000  5204822    0    132      74     1  0.04618  0.09218
   453   0  1.003  5204371    0    127      87     2  0.04799  0.09859
   454   0  1.003  5205287    0    133      72     6  0.05012  0.08264
   455   0  1.002  5205905    0    129      69     1  0.04850  0.13240
   456   0  1.005  5204986    0    126      81     2  0.04439  0.11108
   457   0  1.008  5204562    0    135      82     1  0.04699  0.09508
   458   0  1.004  5203946    0    124      89     1  0.04669  0.09002
   459   0  1.006  5204801    0    121      80     1  0.04505  0.10810
   460   0  1.007  5205018    0    112      76     1  0.04812  0.10627

   450 OUT OF    453 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.56
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.13
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.415    32.132    34.506  88.633 108.468 111.880
 REC. CELL PARAMETERS   0.041669  0.033731  0.030731  83.857  70.618  67.354
 COORDINATES OF UNIT CELL A-AXIS   -14.352   -14.209    18.539
 COORDINATES OF UNIT CELL B-AXIS    14.397    25.590    13.052
 COORDINATES OF UNIT CELL C-AXIS   -21.500    21.470   -16.352
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.026    -0.001    -0.022
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.038
 LAB COORDINATES OF ROTATION AXIS -0.999991 -0.004096  0.001370
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000713  0.005141  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.34   1244.11
 DETECTOR ORIGIN (PIXELS) AT                     1265.40   1237.36
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.09
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     450 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   48     25     55     64     55     27     50     60     51
 SIGMAB (degree)       0.035  0.036  0.036  0.035  0.034  0.035  0.035  0.035  0.035
 SIGMAR (degree)       0.030  0.027  0.029  0.030  0.029  0.026  0.028  0.029  0.028


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    450 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  3  7 10  8  4  2  1  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  1  2  3  1  1  0  0  0     0  0  0  1  5 11 17 15  8  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  1  3  7 10  9  5  2  1  0  0
     1  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0
     0  0  0  1  2  5  6  4  2  1  0  0  0     0  0  0  2  5 10 11  7  3  1  0  0  0     0  0  0  1  4  6  6  4  2  0  0  0  0
     0  0  0  3 12 30 39 29 14  6  2  1  0     0  0  1  6 26 56 62 42 19  8  3  1  0     0  0  1  5 19 36 36 23 10  4  1  0  0
     0  0  1  4 21 47 63 51 27 11  4  1  0     0  0  1 10 43 88100 71 35 14  5  1  0     0  0  1  8 31 56 57 37 17  7  2  1  0
     0  0  0  4 13 26 35 28 15  6  2  1  0     0  0  1  7 26 48 55 38 20  8  3  1  0     0  0  1  6 18 29 30 19 10  4  1  0  0
     0  0  0  1  3  6  7  5  3  1  0  0  0     0  0  0  2  6 10 11  8  4  2  0  0  0     0  0  0  1  4  6  6  4  2  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  8  8  5  2  1  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  8 13 12  7  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  7  7  4  2  1  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    594 REFLECTION INTENSITIES INTEGRATED


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
   461   0  1.006  5204550    0    130      83     2  0.04994  0.09485
   462   0  1.007  5204545    0    136      81     1  0.04789  0.10008
   463   0  1.007  5203648    0    146      87     1  0.04595  0.08603
   464   0  1.006  5203434    0    151      96     1  0.04577  0.09102
   465   0  1.007  5203610    0    141      90     3  0.04544  0.07308
   466   0  1.010  5205265    0    133      74     2  0.04383  0.07622
   467   0  1.011  5206001    0    138      69     2  0.05067  0.09802
   468   0  1.013  5204926    0    135      79     3  0.04811  0.08181
   469   0  1.013  5204265    0    138      85     3  0.04661  0.09013
   470   0  1.012  5204884    0    136      75     3  0.04934  0.10672
   471   0  1.013  5204745    0    138      77     6  0.04708  0.07678
   472   0  1.015  5203983    0    134      84     3  0.04903  0.09198
   473   0  1.014  5204797    0    142      73     6  0.04828  0.07543
   474   0  1.014  5205172    0    138      74     6  0.04987  0.06778
   475   0  1.012  5204412    0    139      80     6  0.04381  0.06485
   476   0  1.012  5203642    0    139      90     2  0.04884  0.09332
   477   0  1.011  5204018    0    134      85     1  0.04862  0.09595
   478   0  1.015  5204899    0    139      75     8  0.04644  0.07006
   479   0  1.018  5204057    0    136      85     4  0.04364  0.08264
   480   0  1.017  5204832    0    131      76     1  0.04739  0.10360

   504 OUT OF    505 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.55
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.408    32.137    34.516  88.639 108.457 111.875
 REC. CELL PARAMETERS   0.041677  0.033724  0.030720  83.856  70.627  67.358
 COORDINATES OF UNIT CELL A-AXIS   -14.354   -14.203    18.531
 COORDINATES OF UNIT CELL B-AXIS    14.396    25.594    13.058
 COORDINATES OF UNIT CELL C-AXIS   -21.506    21.474   -16.362
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.033     0.010    -0.021
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.038
 LAB COORDINATES OF ROTATION AXIS -0.999991 -0.004096  0.001370
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000526  0.005206  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.22   1244.07
 DETECTOR ORIGIN (PIXELS) AT                     1265.53   1237.24
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.14
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     504 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   54     35     59     67     53     26     58     76     63
 SIGMAB (degree)       0.035  0.036  0.036  0.035  0.034  0.035  0.034  0.035  0.035
 SIGMAR (degree)       0.029  0.025  0.028  0.030  0.027  0.027  0.028  0.030  0.028


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    504 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  1  4  9 13 11  5  2  1  0  0
     0  0  0  0  0  0  1  0  0  0  1  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  1  7 16 21 18 10  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  1  4 10 13 11  6  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  1  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  7  8  5  3  1  0  0  0     0  0  0  2  5 11 13  7  3  1  0  0  0     0  0  0  1  3  5  6  3  2  0  0  0  0
     0  0  1  3 15 35 45 34 16  7  2  0  0     0  0  1  6 26 54 63 43 19  8  3  1  0     0  0  1  4 15 28 30 20  8  4  1  0  0
     0  0  1  6 27 59 74 59 30 12  4  1  0     0  0  1 11 47 88100 75 36 14  5  1  0     0  0  1  7 27 45 47 33 16  6  2  1  0
     0  0  0  4 15 33 43 33 16  7  2  1  0     0  0  1  7 26 50 57 40 19  8  3  1  0     0  0  1  4 14 25 26 18  9  3  1  0  0
     0  0  0  1  3  6  9  7  4  1  0  0  0     0  0  0  2  5  9 12  8  5  1  0  0  0     0  0  0  1  3  5  6  4  2  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  6  6  3  1  0  0  0  0     1  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6  9  8  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  5  5  3  1  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    642 REFLECTION INTENSITIES INTEGRATED


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
   481   0  1.016  5205418    0    136      74     0  0.04802  0.10571
   482   0  1.016  5205417    0    136      77     0  0.04333  0.09270
   483   0  1.012  5205247    0    139      76     3  0.04730  0.09172
   484   0  1.015  5205315    0    134      70     1  0.04811  0.10038
   485   0  1.019  5206253    0    130      65     0  0.05267  0.10641
   486   0  1.019  5205751    0    138      72     0  0.04729  0.09990
   487   0  1.020  5204144    0    132      81     1  0.04678  0.09839
   488   0  1.021  5204404    0    147      80     0  0.04958  0.09685
   489   0  1.022  5203844    0    147      85     2  0.04844  0.08162
   490   0  1.024  5204418    0    146      81     3  0.04589  0.08266
   491   0  1.025  5204093    0    136      79     8  0.04582  0.05873
   492   0  1.029  5203797    0    133      83     2  0.04651  0.09824
   493   0  1.027  5205202    0    125      67     3  0.04734  0.10871
   494   0  1.027  5205087    0    127      67     2  0.04903  0.10261
   495   0  1.028  5205056    0    130      76     5  0.04459  0.08459
   496   0  1.025  5204012    0    130      84     1  0.04603  0.10400
   497   0  1.027  5203815    0    128      83     3  0.04543  0.08283
   498   0  1.026  5204404    0    132      76     2  0.04547  0.08871
   499   0  1.025  5204216    0    132      78     3  0.04610  0.09491
   500   0  1.031  5203609    0    129      85     4  0.04641  0.08219

   470 OUT OF    472 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.52
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.07
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.413    32.133    34.505  88.636 108.461 111.871
 REC. CELL PARAMETERS   0.041667  0.033727  0.030730  83.860  70.625  67.364
 COORDINATES OF UNIT CELL A-AXIS   -14.349   -14.209    18.539
 COORDINATES OF UNIT CELL B-AXIS    14.391    25.593    13.055
 COORDINATES OF UNIT CELL C-AXIS   -21.507    21.465   -16.350
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.036    -0.008    -0.020
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.038
 LAB COORDINATES OF ROTATION AXIS -0.999991 -0.004096  0.001370
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000703  0.005088  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.36   1244.10
 DETECTOR ORIGIN (PIXELS) AT                     1265.44   1237.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.11
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     470 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   46     29     54     67     59     31     57     63     56
 SIGMAB (degree)       0.035  0.036  0.035  0.035  0.034  0.035  0.034  0.034  0.035
 SIGMAR (degree)       0.029  0.025  0.028  0.030  0.028  0.025  0.027  0.030  0.028


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    470 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  1  3  7 10  8  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0     0  0  0  1  5 13 17 14  8  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  3  7 10  8  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  1  2  5  6  5  2  1  0  0  0     0  0  0  1  5 10 11  7  3  1  0  0  0     0  0  0  1  3  6  6  3  1  0  0  0  0
     0  0  0  3 12 28 38 29 14  6  2  0  0     0  0  1  6 26 53 59 40 18  7  3  0  0     0  0  1  5 18 31 31 20  9  3  1  0  0
     0  0  1  5 21 51 65 49 26 11  4  1  0     0  0  1 10 44 89100 69 34 14  5  1  0     0  0  1  8 30 51 52 35 17  6  2  0  0
     0  0  0  3 12 29 36 28 15  6  2  1  0     0  0  1  7 24 50 57 40 20  8  3  0  0     0  0  0  5 16 29 31 20 10  4  1  0  0
     0  0  0  1  3  5  6  5  3  1  0  0  0     0  0  0  2  6  9 10  8  4  2  0  0  0     0  0  0  1  4  6  5  4  2  1  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  5  6  3  1  1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 10  9  6  3  1  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  6  6  4  2  1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
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
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  2  5 11 13  9  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  4  3  0  0  0  0  0     0  0  0  2  8 16 19 16  8  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  0  0  0  0  0     0  0  0  0  4  9 12 10  5  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4  8  9  5  3  0  0  0  0     0  0  0  2  6 12 13  7  3  0  0  0  0     0  0  0  2  4  7  7  4  2  0  0  0  0
     0  0  0  5 18 37 44 31 14  5  0  0  0     0  0  0  8 31 60 64 42 17  6  0  0  0     0  0  0  6 20 35 35 21  8  3  0  0  0
     0  0  0  6 27 57 67 52 25  9  0  0  0     0  0  2 11 48 92100 72 33 12  0  0  0     0  0  0  9 32 55 55 37 16  6  0  0  0
     0  0  0  3 15 32 41 31 14  5  0  0  0     0  0  0  6 27 55 63 44 19  7  0  0  0     0  0  0  5 18 33 35 22 10  3  0  0  0
     0  0  0  0  3  7  9  6  3  0  0  0  0     0  0  0  0  6 13 15 10  4  0  0  0  0     0  0  0  0  4  8  9  5  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  2  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6  9  8  4  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  9 13 12  7  3  0  0  0  0     0  0  0  0  0  2  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  5  8  8  5  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
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
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  6  5  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  4  5  7  7  5  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  2  5  9  9  7  4  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  5  6  4  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4  7  9  8  4  3  0  0  0  0     0  0  0  4  9 14 14  8  4  2  0  0  0     0  0  0  2  5  9  9  6  3  0  0  0  0
     0  0  0  6 17 30 30 20 10  4  0  0  0     0  0  2  9 32 60 63 42 20  8  0  0  0     0  0  0  6 22 43 46 29 14  6  0  0  0
     0  0  0  5 20 39 44 33 17  7  0  0  0     0  0  2 12 46 89100 73 36 14  5  0  0     0  0  0  9 35 68 73 50 24  9  0  0  0
     0  0  0  3 10 21 25 19 10  4  0  0  0     0  0  0  6 27 51 59 42 21  9  0  0  0     0  0  0  5 21 38 41 29 14  6  0  0  0
     0  0  0  0  3  4  5  4  2  0  0  0  0     0  0  0  3  5 11 13  9  4  2  0  0  0     0  0  0  0  4  9  9  6  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4  8  8  6  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  7 13 13  9  4  3  0  0  0     0  0  0  0  0  2  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  7  8  5  3  2  0  0  0     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  3  3  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
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
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  0  0  0  0  0  0     0  0  0  0  3  6  8  7  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  0  0  0  0  0     0  0  0  0  5 11 15 13  7  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  3  0  0  0  0  0     0  0  0  0  3  7 11  9  5  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  4  5  4  2  0  0  0  0     0  0  0  0  3  7  8  5  3  0  0  0  0     0  0  0  0  3  5  5  3  2  0  0  0  0
     0  0  0  3 10 25 30 22 11  5  0  0  0     0  0  0  5 22 45 48 32 16  6  0  0  0     0  0  0  4 16 30 29 18  8  4  0  0  0
     0  0  0  4 21 47 59 47 25 11  0  0  0     0  0  0 10 45 89100 73 37 15  6  0  0     0  0  0  9 35 61 61 41 20  8  4  0  0
     0  0  0  3 15 32 42 34 19  8  0  0  0     0  0  0  8 33 62 72 54 28 12  5  0  0     0  0  0  7 26 44 46 32 16  7  0  0  0
     0  0  0  0  4  8 10  8  5  3  0  0  0     0  0  0  3  8 15 17 13  7  3  0  0  0     0  0  0  3  6 11 11  8  4  2  0  0  0
     0  0  0  0  0  2  2  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4  7  6  4  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3 10 15 14  8  4  2  0  0  0     0  0  0  0  0  2  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  8 12 11  7  4  3  0  0  0     0  0  0  0  0  2  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  4  3  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  4 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  4  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  0  0  0  0  0     0  0  0  0  4 11 17 14  7  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  4  3  0  0  0  0     0  0  0  0  5 16 24 22 12  5  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  3  0  0  0  0  0     0  0  0  0  3  7 12 11  6  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  2  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  7 10  7  3  0  0  0  0     0  0  0  0  5 11 13  8  4  0  0  0  0     0  0  0  0  4  7  6  4  2  0  0  0  0
     0  0  0  3 15 38 52 40 20  8  0  0  0     0  0  0  6 28 58 66 44 21  9  0  0  0     0  0  0  6 20 35 34 21  9  4  0  0  0
     0  0  0  4 22 56 78 63 33 14  6  0  0     0  0  0 10 43 88100 71 35 14  5  0  0     0  0  0  9 31 54 53 34 17  6  0  0  0
     0  0  0  3 11 27 38 31 16  7  0  0  0     0  0  0  5 22 44 51 36 18  7  0  0  0     0  0  0  5 16 28 28 18  9  3  0  0  0
     0  0  0  0  2  4  6  5  3  0  0  0  0     0  0  0  0  4  8  9  6  3  0  0  0  0     0  0  0  0  4  5  5  4  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  4  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  6 10  8  4  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4 10 14 12  7  3  2  0  0  0     0  0  0  0  0  3  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  5  7  7  4  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  2  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  5 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  2  7 12 11  5  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  3  0  0  0  0  0     0  0  0  0  3 11 20 19 11  5  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  2  5 10  9  5  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  4  6  5  2  0  0  0  0     0  0  0  0  4  7  9  5  3  0  0  0  0     0  0  0  0  3  5  4  3  0  0  0  0  0
     0  0  0  2 10 29 43 35 17  7  0  0  0     0  0  0  5 24 52 59 41 18  8  0  0  0     0  0  0  5 19 32 29 18  8  3  0  0  0
     0  0  0  3 17 50 74 61 32 13  0  0  0     0  0  0  9 42 90100 69 33 13  5  0  0     0  0  0  9 32 53 49 30 14  5  0  0  0
     0  0  0  2  8 23 35 27 14  6  0  0  0     0  0  0  5 20 42 47 30 14  6  0  0  0     0  0  0  4 14 24 22 13  6  3  0  0  0
     0  0  0  0  0  3  4  3  2  0  0  0  0     0  0  0  0  3  6  6  4  2  0  0  0  0     0  0  0  0  5  3  3  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  6  8  6  3  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  5  8 12  9  5  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
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
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  4  6  6  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  7 11 10  6  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  2  4  6  6  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  4  5  4  2  0  0  0  0     0  0  0  2  6 11 10  6  3  0  0  0  0     0  0  0  3  5  6  4  3  0  0  0  0  0
     0  0  0  2  9 24 33 25 12  5  0  0  0     0  0  2  8 32 61 60 37 16  6  0  0  0     0  0  0  8 24 33 27 14  6  3  0  0  0
     0  0  0  3 16 41 56 44 23  9  0  0  0     0  0  3 13 55100100 65 30 12  5  0  0     0  0  3 13 39 54 44 26 11  4  0  0  0
     0  0  0  2  8 22 31 24 12  5  0  0  0     0  0  0  7 28 54 56 35 16  7  0  0  0     0  0  0  6 20 29 24 14  6  3  0  0  0
     0  0  0  0  2  4  5  5  3  0  0  0  0     0  0  0  2  5 10 10  7  3  0  0  0  0     0  0  0  0  4  5  4  3  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  5  5  4  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4  6  8  6  3  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  4  3  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
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
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  7  8  6  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  3  7 14 17 15  9  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  2  0  0  0  0  0     0  0  0  0  4  8 11 11  7  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  3  4  3  0  0  0  0  0     0  0  0  0  4  7  6  4  2  0  0  0  0     0  0  0  0  3  4  3  3  2  0  0  0  0
     0  0  0  3 12 27 32 23 11  5  0  0  0     0  0  0  6 24 47 50 32 15  6  0  0  0     0  0  0  4 15 25 24 14  6  3  0  0  0
     0  0  0  6 27 55 66 53 29 12  5  0  0     0  0  2 12 49 91100 72 35 15  5  0  0     0  0  0  8 30 50 50 33 15  6  0  0  0
     0  0  0  4 15 31 41 36 20  9  4  0  0     0  0  0  7 26 51 60 47 24 10  4  0  0     0  0  0  4 15 28 30 22 10  4  0  0  0
     0  0  0  0  3  5  8  8  5  2  0  0  0     0  0  0  0  4  8 11  9  5  2  0  0  0     0  0  0  0  3  4  5  4  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  4  4  3  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  6  9  9  6  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  4  6  6  4  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  8 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  4  3  3  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  3  3  2  0  0  0  0  0     0  0  0  4  8 12 13 10  5  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  3  4  3  3  0  0  0  0     0  0  0  4 12 22 25 20 11  5  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  2  0  0  0  0     0  0  0  2  6 12 16 14  8  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  5  6  6  4  2  0  0  0  0     0  0  0  3  5  8  8  5  2  0  0  0  0     0  0  0  0  3  4  4  3  0  0  0  0  0
     0  0  2  7 20 35 40 29 14  6  0  0  0     0  0  2  8 24 46 52 37 17  7  0  0  0     0  0  0  4 13 25 28 19  9  4  0  0  0
     0  0  3 10 34 65 76 60 32 14  5  0  0     0  0  3 12 44 85100 75 38 16  6  0  0     0  0  0  7 24 47 53 39 19  8  0  0  0
     0  0  0  6 19 38 49 40 22  9  4  0  0     0  0  0  7 24 50 64 50 25 11  4  0  0     0  0  0  4 13 28 35 25 13  5  0  0  0
     0  0  0  0  3  6  9  8  5  2  0  0  0     0  0  0  0  4  8 11 10  5  2  0  0  0     0  0  0  0  3  5  6  5  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  6  6  4  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 11 12  8  4  2  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  7  8  6  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  2  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  9 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  6  5  4  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  6  4  3  3  0  0  0  0  0     0  0  0  5 10 14 13  9  4  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  4  4  3  2  0  0  0  0     0  0  0  3  9 17 18 14  7  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  2  0  0  0  0  0     0  0  0  0  3  7  9  7  4  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  6 11 13 11  6  3  0  0  0  0     0  0  0  5 11 16 15  8  4  2  0  0  0     0  0  0  3  5  8  9  5  3  0  0  0  0
     0  0  3 11 28 43 44 32 15  6  0  0  0     0  0  3 12 35 62 70 49 23  9  0  0  0     0  0  0  6 19 37 43 30 14  5  0  0  0
     0  0  0  8 28 54 64 49 25 11  5  0  0     0  0  2 10 39 82100 76 38 16  6  0  0     0  0  0  5 22 50 62 46 22  9  0  0  0
     0  0  0  3 10 23 31 24 12  5  0  0  0     0  0  0  3 15 35 49 36 18  8  3  0  0     0  0  0  3  9 21 29 21 10  4  0  0  0
     0  0  0  0  2  3  4  4  2  0  0  0  0     0  0  0  0  3  5  7  6  3  0  0  0  0     0  0  0  0  2  3  5  4  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4  8 10  7  4  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  5 11 14 10  5  3  0  0  0     0  0  0  0  0  2  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  5  7  5  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 STANDARD DEVIATION (DEGREES) OF BEAM DIVERGENCE
 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1    0.025  0.025  0.025  0.024  0.023  0.024  0.024  0.025  0.025
     2    0.024  0.025  0.025  0.024  0.023  0.024  0.024  0.024  0.025
     3    0.024  0.025  0.026  0.024  0.023  0.025  0.024  0.025  0.025
     4    0.024  0.025  0.025  0.024  0.023  0.024  0.024  0.025  0.025
     5    0.024  0.025  0.025  0.024  0.023  0.024  0.025  0.025  0.025
     6    0.024  0.025  0.024  0.024  0.023  0.025  0.025  0.025  0.025
     7    0.025  0.025  0.025  0.024  0.023  0.024  0.025  0.025  0.025
     8    0.025  0.025  0.025  0.024  0.023  0.025  0.025  0.026  0.026
     9    0.025  0.025  0.025  0.025  0.024  0.025  0.025  0.026  0.025
    10    0.025  0.025  0.025  0.025  0.024  0.025  0.025  0.026  0.026
    11    0.025  0.026  0.025  0.025  0.024  0.025  0.025  0.026  0.026
    12    0.025  0.025  0.026  0.025  0.024  0.025  0.025  0.026  0.026
    13    0.025  0.025  0.025  0.025  0.024  0.025  0.025  0.026  0.026
    14    0.025  0.025  0.026  0.025  0.024  0.025  0.026  0.026  0.026
    15    0.025  0.026  0.026  0.025  0.024  0.025  0.025  0.026  0.026
    16    0.025  0.026  0.026  0.025  0.024  0.025  0.025  0.026  0.026
    17    0.025  0.026  0.026  0.025  0.024  0.025  0.025  0.026  0.026
    18    0.024  0.025  0.026  0.025  0.024  0.025  0.025  0.026  0.026
    19    0.025  0.026  0.026  0.025  0.024  0.025  0.025  0.026  0.026
    20    0.025  0.025  0.026  0.025  0.024  0.025  0.025  0.025  0.025
    21    0.025  0.026  0.026  0.025  0.024  0.025  0.025  0.026  0.026
    22    0.025  0.026  0.026  0.025  0.024  0.025  0.025  0.025  0.026
    23    0.025  0.025  0.026  0.025  0.024  0.025  0.025  0.025  0.025
    24    0.024  0.026  0.026  0.025  0.024  0.025  0.024  0.025  0.025
    25    0.025  0.025  0.026  0.025  0.024  0.025  0.025  0.024  0.025


 STANDARD DEVIATION (DEGREES) OF REFLECTING RANGE
 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1    0.030  0.032  0.030  0.030  0.028  0.023  0.028  0.032  0.031
     2    0.028  0.029  0.028  0.030  0.028  0.024  0.027  0.029  0.030
     3    0.029  0.026  0.029  0.033  0.030  0.023  0.026  0.031  0.031
     4    0.029  0.025  0.028  0.031  0.029  0.024  0.028  0.032  0.029
     5    0.029  0.026  0.028  0.031  0.029  0.024  0.027  0.029  0.028
     6    0.030  0.032  0.029  0.030  0.029  0.025  0.027  0.029  0.030
     7    0.029  0.037  0.031  0.030  0.028  0.022  0.026  0.029  0.031
     8    0.029  0.030  0.029  0.030  0.029  0.024  0.027  0.030  0.030
     9    0.031  0.026  0.027  0.030  0.028  0.027  0.029  0.034  0.034
    10    0.029  0.027  0.028  0.030  0.029  0.024  0.027  0.030  0.029
    11    0.030  0.028  0.028  0.032  0.030  0.023  0.026  0.030  0.029
    12    0.029  0.026  0.028  0.032  0.029  0.023  0.027  0.031  0.029
    13    0.030  0.026  0.028  0.031  0.030  0.025  0.029  0.035  0.031
    14    0.030  0.026  0.026  0.029  0.028  0.025  0.030  0.035  0.031
    15    0.030  0.034  0.032  0.030  0.028  0.024  0.027  0.030  0.030
    16    0.029  0.026  0.029  0.029  0.028  0.024  0.029  0.033  0.030
    17    0.029  0.024  0.027  0.029  0.028  0.024  0.028  0.032  0.029
    18    0.028  0.026  0.028  0.030  0.028  0.022  0.025  0.030  0.028
    19    0.028  0.029  0.028  0.029  0.027  0.024  0.029  0.030  0.029
    20    0.029  0.024  0.028  0.031  0.028  0.023  0.027  0.030  0.027
    21    0.028  0.023  0.028  0.030  0.028  0.028  0.028  0.029  0.027
    22    0.029  0.023  0.027  0.030  0.028  0.024  0.027  0.031  0.028
    23    0.028  0.027  0.028  0.030  0.028  0.024  0.027  0.029  0.028
    24    0.028  0.023  0.027  0.032  0.028  0.025  0.028  0.031  0.028
    25    0.027  0.024  0.027  0.029  0.027  0.024  0.027  0.030  0.027


 ***** SUGGESTED VALUES FOR INPUT PARAMETERS *****
 BEAM_DIVERGENCE=   0.24954  BEAM_DIVERGENCE_E.S.D.=   0.02495
 REFLECTING_RANGE=  0.19821  REFLECTING_RANGE_E.S.D.=  0.02832



 NUMBER OF STRONG REFLECTIONS USED FOR LEARNING REFERENCE PROFILES

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1       46     33     56     71     56     30     57     69     55
     2       53     33     65     66     53     27     56     68     56
     3       51     33     53     66     58     32     53     68     63
     4       52     33     59     68     54     28     60     70     62
     5       47     31     57     69     58     29     48     62     52
     6       48     32     57     65     50     26     54     70     55
     7       46     36     57     66     57     32     53     68     59
     8       52     33     62     69     56     29     57     69     56
     9       51     33     55     64     56     28     50     65     58
    10       46     30     58     68     57     29     60     72     57
    11       50     30     59     67     54     28     48     65     56
    12       52     35     56     65     58     34     56     64     59
    13       50     35     66     71     57     28     57     68     55
    14       48     31     60     64     57     30     46     62     55
    15       50     29     59     76     56     34     71     70     53
    16       49     30     53     60     50     28     42     64     58
    17       50     33     61     73     56     29     60     71     56
    18       46     29     56     69     60     34     59     63     57
    19       45     34     55     55     56     23     43     66     57
    20       50     30     58     66     54     28     56     68     60
    21       52     32     60     75     56     35     61     70     58
    22       49     35     62     70     51     25     51     68     60
    23       48     25     55     64     55     27     50     60     51
    24       54     35     59     67     53     26     58     76     63
    25       46     29     54     67     59     31     57     63     56
 margin    1232    798   1452   1682   1389    730   1365   1678   1425



 NUMBER OF GRID POINTS IN REFERENCE PROFILES WITH  SIGNAL >   2.% OF MAXIMUM

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1      136    124    123    118    107    105    116    133    128
     2      121    125    124    125    108    104    112    119    121
     3      133    124    139    134    113    103    108    135    127
     4      129    120    117    130    106    108    109    135    127
     5      128    118    121    129    111    104    117    130    120
     6      129    129    117    121    107    115    118    123    125
     7      129    125    130    122    103     98    115    131    133
     8      126    124    127    122    109    103    117    136    129
     9      138    121    120    123    110    111    126    153    133
    10      129    118    124    124    109    110    121    139    131
    11      140    134    123    131    105    105    117    140    132
    12      131    120    135    135    108    103    121    140    131
    13      135    121    122    128    106    103    124    154    133
    14      134    120    125    125    106    111    126    149    137
    15      139    130    136    127    109    110    114    140    142
    16      131    120    133    128    106    109    126    145    134
    17      137    122    128    127    113    114    120    148    135
    18      118    118    130    127    107    102    111    136    128
    19      128    124    124    124    111    112    123    139    131
    20      130    120    133    133    116    108    118    136    124
    21      126    119    130    129    113    108    120    130    129
    22      130    122    122    132    114    111    119    134    133
    23      129    122    131    127    111    110    118    125    123
    24      122    122    131    135    111    109    113    135    125
    25      123    118    126    129    115    107    115    126    121
 margin     153    157    166    152    128    129    143    164    165

 MINIMUM INTEGRATION VOLUME IS   98 GRID POINTS



 CORRELATIONS BETWEEN REFERENCE PROFILES WITH CENTRAL PROFILE #1

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1     1.00   0.86   0.95   0.99   0.98   0.97   0.97   0.98   0.94
     2     1.00   0.92   0.97   0.97   0.97   0.98   0.99   0.99   0.94
     3     1.00   0.98   0.97   0.96   0.97   0.94   0.98   0.98   0.95
     4     1.00   0.98   0.98   0.98   0.97   0.97   0.98   0.98   0.95
     5     1.00   0.95   0.98   0.98   0.98   0.97   0.97   0.98   0.94
     6     1.00   0.92   0.96   0.98   0.97   0.98   0.98   0.97   0.96
     7     1.00   0.79   0.92   0.98   0.97   0.96   0.98   0.99   0.93
     8     1.00   0.89   0.93   0.98   0.97   0.97   0.96   0.96   0.97
     9     1.00   0.97   0.98   0.97   0.96   0.92   0.98   0.93   0.90
    10     1.00   0.93   0.96   0.98   0.97   0.97   0.99   0.98   0.95
    11     1.00   0.92   0.96   0.97   0.96   0.97   0.97   0.97   0.95
    12     1.00   0.94   0.97   0.96   0.95   0.97   0.97   0.97   0.95
    13     1.00   0.96   0.96   0.98   0.97   0.97   0.98   0.96   0.96
    14     1.00   0.95   0.94   0.98   0.97   0.97   0.94   0.91   0.95
    15     1.00   0.85   0.91   0.98   0.98   0.96   0.95   0.95   0.95
    16     1.00   0.96   0.90   0.96   0.94   0.95   0.94   0.93   0.95
    17     1.00   0.96   0.96   0.98   0.97   0.97   0.98   0.96   0.95
    18     1.00   0.95   0.97   0.98   0.97   0.98   0.99   0.98   0.95
    19     1.00   0.91   0.96   0.98   0.97   0.97   0.97   0.97   0.95
    20     1.00   0.97   0.96   0.98   0.98   0.97   0.97   0.97   0.96
    21     1.00   0.98   0.98   0.97   0.97   0.92   0.95   0.97   0.96
    22     1.00   0.95   0.97   0.97   0.98   0.97   0.98   0.98   0.95
    23     1.00   0.93   0.96   0.98   0.97   0.98   0.98   0.98   0.96
    24     1.00   0.96   0.96   0.98   0.97   0.98   0.95   0.98   0.97
    25     1.00   0.96   0.97   0.97   0.98   0.97   0.99   0.98   0.97
 margin    1.00   0.96   0.97   0.99   0.98   0.98   0.99   0.98   0.97



 CORRELATIONS BETWEEN REFERENCE PROFILES WITH RUN-AVERAGED PROFILES

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1     1.00   0.95   0.98   0.99   1.00   0.99   1.00   0.99   0.99
     2     1.00   0.99   1.00   0.99   0.99   0.99   1.00   0.99   0.98
     3     0.98   0.92   0.93   0.97   0.98   0.99   0.99   0.99   0.97
     4     0.99   0.97   0.98   0.99   1.00   1.00   0.99   0.99   0.99
     5     1.00   1.00   0.99   0.99   1.00   1.00   0.99   0.99   0.99
     6     0.99   0.96   0.98   0.99   0.99   0.99   0.98   0.98   0.97
     7     0.99   0.88   0.96   1.00   0.99   0.99   1.00   0.99   0.96
     8     0.99   0.98   0.99   1.00   1.00   0.99   0.99   1.00   0.99
     9     1.00   0.88   0.98   0.99   0.98   0.95   1.00   0.97   0.90
    10     1.00   0.99   0.99   0.99   1.00   0.99   0.99   0.99   0.98
    11     0.99   0.99   0.98   0.98   0.98   0.99   0.99   0.99   0.99
    12     1.00   1.00   0.99   0.99   0.99   0.99   0.99   0.99   0.99
    13     0.99   0.98   0.99   0.99   0.99   0.98   0.98   0.98   0.97
    14     1.00   1.00   0.99   0.99   0.99   1.00   0.97   0.96   0.97
    15     0.99   0.92   0.94   0.95   0.98   1.00   0.99   1.00   0.98
    16     1.00   1.00   0.95   0.95   0.95   0.98   0.99   0.99   0.99
    17     1.00   0.99   1.00   0.99   1.00   1.00   0.99   0.99   0.98
    18     1.00   1.00   1.00   1.00   0.99   0.99   1.00   1.00   1.00
    19     0.99   0.97   0.99   0.99   0.99   0.99   0.99   0.99   0.98
    20     0.99   0.97   0.99   0.99   0.99   1.00   0.99   0.99   0.99
    21     0.99   0.98   0.99   1.00   0.98   0.90   0.91   0.98   0.99
    22     1.00   0.98   0.99   0.99   0.99   1.00   1.00   1.00   0.98
    23     0.99   0.99   1.00   1.00   0.99   1.00   1.00   0.99   0.98
    24     0.99   0.97   0.99   0.99   0.98   0.98   0.99   0.99   0.98
    25     1.00   0.99   1.00   1.00   1.00   0.99   0.99   1.00   0.99
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

       -1... 0      -0.136      0.001        -2.          2
        0... 1       0.254      0.094         6.         26
        1... 2       0.465      0.097        16.         65
        2... 3       0.611      0.082        31.        105
        3... 6       0.706      0.088        80.        421
        6... 9       0.762      0.083       185.        579
        9...12       0.772      0.081       344.        707
       12...15       0.771      0.081       551.        890
       15...18       0.780      0.078       811.        924
       18...21       0.783      0.081      1152.        961
       21...24       0.785      0.078      1580.        890
       24...27       0.788      0.078      2095.        873
       27...30       0.791      0.074      2723.        784
       30...33       0.789      0.077      3377.        722
       33...36       0.792      0.079      4438.        665
       36...39       0.787      0.085      5900.        736
       39...42       0.788      0.084      8033.        684
       42...45       0.788      0.082     11651.        803
       45...48       0.792      0.080     21525.       1184
       48...51       0.800      0.097     51044.       1064

   13085 REFLECTIONS SAVED ON FILE "INTEGRATE.HKL"


 total elapsed wall-clock time for INTEGRATE      17.3 sec

 TASK     cpu time (sec)    elapsed wall-clock time (sec)
   1         273.9                    17.2
