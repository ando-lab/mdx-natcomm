
 ***** INTEGRATE ***** (VERSION Mar 15, 2019  BUILT=20190806)   6-Jan-2020

 INPUT PARAMETER VALUES
 ----------------------
 MAXIMUM_NUMBER_OF_PROCESSORS=   0  (0: automatic choice of OpenMP threads)
 number of OpenMP threads used  20
 MAXIMUM_NUMBER_OF_JOBS=   0  (0:automatic choice of forked main tasks)
 number of forked INTEGRATE tasks   1
 CLUSTER_NODES=
 SECONDS=     0
 NAME_TEMPLATE_OF_DATA_FRAMES=~/data/lysozyme/nitrate/ubatch4/lys_nitr_10_2_????.cbf CBF
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
 UNIT_CELL_CONSTANTS=    27.44    32.15    34.52  88.650 108.469 111.889
 ROTATION_AXIS= -1.00000 -0.00172 -0.00012
 OSCILLATION_RANGE=  0.100000 DEGREES
 STARTING_ANGLE= 295.000 STARTING_FRAME=       1
 X-RAY_WAVELENGTH=   0.97700 ANGSTROM
 INCIDENT_BEAM_DIRECTION=  0.000698  0.005184  1.023528
 NUMBER OF DETECTOR SEGMENTS   1
 NX=  2463 NY=  2527 QX= 0.17200 QY= 0.17200
 ORGX=   1265.49  ORGY=   1237.28
 DETECTOR_DISTANCE=   231.294
 DIRECTION_OF_DETECTOR_X-AXIS=   1.00000   0.00000   0.00000
 DIRECTION_OF_DETECTOR_Y-AXIS=   0.00000   1.00000   0.00000
 SEGMENT=    1 2463    1 2527  ! x1,x2,y1,y2 (pixel units)
    SEGMENT_ORGX=      0.00   SEGMENT_ORGY=      0.00
    SEGMENT_DISTANCE=     0.00
    DIRECTION_OF_SEGMENT_X-AXIS=  1.000000  0.000000  0.000000
    DIRECTION_OF_SEGMENT_Y-AXIS=  0.000000  1.000000  0.000000
 RESOLUTION RANGE RECORDED BY DETECTOR (ANGSTROM)    44.530     1.212
 NUMBER OF TRUSTED DETECTOR PIXELS    5219656
 MEAN CONTENTS OF TRUSTED PIXELS IN BACKGROUND TABLE    150.779
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
     1   0  1.041  5211315    0    147      90     2  0.04203  0.09854
     2   0  1.002  5211464    0    143      87     8  0.04369  0.06913
     3   0  1.002  5212002    0    143      84     5  0.04188  0.08211
     4   0  1.002  5211687    0    130      85     5  0.04151  0.07702
     5   0  1.005  5211840    0    117      84     2  0.04304  0.09857
     6   0  1.005  5212478    0    124      81     0  0.04387  0.11326
     7   0  1.005  5212330    0    116      86     3  0.04303  0.11363
     8   0  1.002  5212738    0    117      74     7  0.04789  0.09772
     9   0  1.004  5212939    0    123      74     8  0.04704  0.08427
    10   0  1.003  5212537    0    123      79     4  0.04421  0.09562
    11   0  1.002  5211530    0    124      88     3  0.04652  0.10883
    12   0  1.004  5212229    0    133      82     9  0.04304  0.09480
    13   0  1.005  5211852    0    136      85     5  0.04550  0.11153
    14   0  1.006  5211713    0    136      85     6  0.04540  0.09704
    15   0  1.003  5211417    0    134      86     2  0.04519  0.10502
    16   0  1.005  5211478    0    143      84     4  0.04675  0.09595
    17   0  1.008  5210979    0    141      88     5  0.04708  0.08061
    18   0  1.006  5210333    0    146      91     2  0.04478  0.10808
    19   0  1.005  5210906    0    142      94     5  0.04726  0.09312
    20   0  1.005  5212428    0    128      84     3  0.04489  0.12036

   480 OUT OF    481 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.57
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.433    32.125    34.509  88.665 108.470 111.895
 REC. CELL PARAMETERS   0.041651  0.033745  0.030730  83.813  70.600  67.326
 COORDINATES OF UNIT CELL A-AXIS   -14.361   -14.225    18.546
 COORDINATES OF UNIT CELL B-AXIS    14.489    25.519    13.073
 COORDINATES OF UNIT CELL C-AXIS   -21.484    21.549   -16.278
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.000     0.002     0.010
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001717 -0.000122
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000867  0.005082  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.37   1244.06
 DETECTOR ORIGIN (PIXELS) AT                     1265.24   1237.38
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.23
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     480 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   49     33     60     66     53     27     52     65     57
 SIGMAB (degree)       0.035  0.035  0.034  0.034  0.033  0.034  0.034  0.034  0.034
 SIGMAR (degree)       0.030  0.026  0.029  0.030  0.028  0.025  0.028  0.031  0.029


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    480 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  1  3  8 12  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0     0  0  0  1  6 14 19 17  9  4  2  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  3  7 10  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  3  8  8  5  2  1  0  0  0     0  0  0  2  6 13 13  7  3  1  0  0  0     0  0  0  1  4  7  7  4  1  0  0  0  0
     0  0  0  3 15 34 43 28 15  6  2  1  0     0  0  1  7 31 57 62 37 19  7  2  1  0     0  0  1  5 19 33 32 18  9  3  1  0  0
     0  0  0  5 24 54 68 53 27 11  4  1  0     0  0  1 11 47 89100 69 34 13  5  1  0     0  0  1  8 29 49 52 33 16  6  2  1  0
     0  0  0  3 13 29 36 28 15  6  2  1  0     0  0  1  6 25 49 56 39 20  8  3  1  0     0  0  0  4 15 28 31 21 10  4  1  0  0
     0  0  0  1  2  6  7  5  3  1  0  0  0     0  0  0  1  5 10 11  8  4  2  0  0  0     0  0  0  1  3  7  6  4  2  1  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  6  6  4  2  1  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 10 11  6  3  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  6  7  4  2  1  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    632 REFLECTION INTENSITIES INTEGRATED


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
    21   0  1.001  5212748    0    120      75     3  0.04521  0.12705
    22   0  1.001  5213344    0    121      67     2  0.04766  0.13800
    23   0  1.003  5212644    0    127      71     5  0.04505  0.08561
    24   0  1.002  5212011    0    137      81     0  0.04228  0.10864
    25   0  1.002  5211230    0    145      85     2  0.04161  0.09042
    26   0  1.002  5210740    0    138      94     4  0.04158  0.08700
    27   0  0.999  5211839    0    136      81     3  0.04118  0.09132
    28   0  1.004  5212269    0    129      75     4  0.04382  0.08539
    29   0  1.005  5211966    0    130      83     3  0.04594  0.10553
    30   0  1.003  5211889    0    130      82     1  0.04619  0.12165
    31   0  1.003  5211846    0    140      84     5  0.04358  0.09680
    32   0  1.004  5211097    0    135      96     7  0.04499  0.08717
    33   0  1.008  5211737    0    140      88     5  0.04372  0.08492
    34   0  1.007  5211946    0    138      89     1  0.04177  0.11284
    35   0  1.002  5211056    0    137      94     1  0.04164  0.10957
    36   0  1.003  5211979    0    144      89     1  0.04376  0.11497
    37   0  1.004  5211603    0    141      87     1  0.04191  0.08896
    38   0  1.006  5211187    0    141      92     2  0.04408  0.08883
    39   0  1.006  5211949    0    138      84     2  0.04220  0.10948
    40   0  1.007  5211242    0    143      90     3  0.04234  0.10706

   484 OUT OF    484 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.57
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.11
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.427    32.110    34.508  88.649 108.470 111.905
 REC. CELL PARAMETERS   0.041662  0.033762  0.030731  83.829  70.606  67.322
 COORDINATES OF UNIT CELL A-AXIS   -14.357   -14.226    18.539
 COORDINATES OF UNIT CELL B-AXIS    14.486    25.503    13.068
 COORDINATES OF UNIT CELL C-AXIS   -21.481    21.560   -16.264
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.017    -0.003     0.016
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001717 -0.000122
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.001094  0.004977  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.36   1243.98
 DETECTOR ORIGIN (PIXELS) AT                     1264.92   1237.45
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.15
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     484 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   51     30     54     64     55     30     55     68     58
 SIGMAB (degree)       0.034  0.035  0.035  0.034  0.033  0.034  0.034  0.034  0.033
 SIGMAR (degree)       0.031  0.027  0.029  0.030  0.029  0.025  0.028  0.031  0.031


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    484 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  1  1  0  0  0  0     0  0  0  0  0  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  1  0  0  0     0  0  0  1  4  9 13 10  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  3  2  1  0  0  0     0  0  0  1  6 16 21 18 10  4  2  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  2  1  0  0  0  0     0  0  0  1  3  8 11 10  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  3  8  8  6  2  1  0  0  0     0  0  0  1  5 13 12  7  3  1  0  0  0     0  0  0  1  3  8  7  4  2  0  0  0  0
     0  0  0  3 15 36 43 31 15  6  2  0  0     0  0  1  8 28 55 59 39 18  7  2  0  0     0  0  1  5 16 31 32 21  9  4  1  0  0
     0  0  0  5 26 61 72 55 28 11  4  1  0     0  0  1 11 48 94100 69 34 13  5  1  0     0  0  1  8 30 55 55 36 17  7  2  0  0
     0  0  0  2 13 33 41 32 15  7  2  1  0     0  0  0  6 25 52 59 41 20  8  3  1  0     0  0  1  4 17 32 34 22 11  4  2  0  0
     0  0  0  0  2  5  7  6  3  1  0  0  0     0  0  0  1  4  9 11  9  5  2  0  0  0     0  0  0  1  3  6  8  6  3  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4  7  7  4  2  1  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  2  7 13 13  8  4  2  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  8  8  6  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0
     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    634 REFLECTION INTENSITIES INTEGRATED


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
    41   0  1.007  5210174    0    146      98     9  0.04378  0.08014
    42   0  1.006  5209288    0    141     111     3  0.04226  0.11589
    43   0  1.007  5210787    0    144      92     0  0.04348  0.12592
    44   0  1.009  5211144    0    137      91     2  0.04869  0.14312
    45   0  1.006  5210068    0    127     103     4  0.04797  0.13125
    46   0  1.004  5210116    0    137      99     2  0.04535  0.13862
    47   0  1.006  5209860    0    137     107     1  0.04748  0.15475
    48   0  1.006  5210715    0    136      90     6  0.04810  0.11247
    49   0  1.006  5210819    0    144      91     1  0.04971  0.15801
    50   0  1.003  5210520    0    142      89     3  0.04610  0.09951
    51   0  1.002  5210075    0    147      95     4  0.04302  0.10155
    52   0  1.005  5210250    0    143     101     3  0.04458  0.10354
    53   0  1.006  5210547    0    147     100     3  0.04624  0.12204
    54   0  1.006  5209678    0    148     107     1  0.04496  0.13521
    55   0  1.008  5210329    0    148     103     1  0.04457  0.14614
    56   0  1.009  5211219    0    146      92     7  0.04307  0.09455
    57   0  1.009  5212100    0    134      83     8  0.04158  0.08109
    58   0  1.007  5211980    0    141      86     9  0.04299  0.07286
    59   0  1.009  5212016    0    134      86     2  0.04328  0.10714
    60   0  1.007  5211846    0    141      81     3  0.04011  0.08198

   489 OUT OF    491 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.56
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.436    32.141    34.507  88.681 108.478 111.915
 REC. CELL PARAMETERS   0.041657  0.033735  0.030736  83.784  70.582  67.298
 COORDINATES OF UNIT CELL A-AXIS   -14.362   -14.241    18.538
 COORDINATES OF UNIT CELL B-AXIS    14.490    25.531    13.085
 COORDINATES OF UNIT CELL C-AXIS   -21.485    21.540   -16.283
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.008    -0.004     0.026
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001717 -0.000122
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000409  0.005118  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.28   1244.04
 DETECTOR ORIGIN (PIXELS) AT                     1265.74   1237.32
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.25
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     489 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   52     29     57     63     56     34     54     72     57
 SIGMAB (degree)       0.034  0.035  0.035  0.034  0.033  0.034  0.034  0.033  0.034
 SIGMAR (degree)       0.030  0.027  0.029  0.031  0.029  0.027  0.029  0.031  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    489 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  1  4 11 13 11  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  1  7 16 21 17  9  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  1  4  9 12  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  8  6  3  1  0  0  0     0  0  0  1  6 12 12  7  4  1  0  0  0     0  0  0  1  4  7  6  4  2  0  0  0  0
     0  0  0  3 15 36 44 33 16  7  2  1  0     0  0  1  6 27 57 63 43 20  8  3  1  0     0  0  0  5 19 34 34 21  9  3  1  0  0
     0  0  0  5 25 54 68 54 27 12  4  1  0     0  0  1  9 45 90100 73 34 15  5  2  0     0  0  1  7 34 56 54 36 16  7  2  1  0
     0  0  0  3 13 31 39 29 16  7  3  1  0     0  0  0  5 23 52 61 42 22 10  3  1  0     0  0  0  4 18 32 33 22 11  5  2  0  0
     0  0  0  1  2  6  7  6  3  1  0  0  0     0  0  0  1  5 10 12  9  5  2  1  0  0     0  0  0  1  3  7  7  6  3  1  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  7  4  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 11 10  7  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  7  4  2  1  0  0  0     0  0  0  0  1  1  0  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    655 REFLECTION INTENSITIES INTEGRATED


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
    61   0  1.009  5210841    0    142      93     2  0.03887  0.09016
    62   0  1.008  5211435    0    146      90     4  0.04056  0.08054
    63   0  1.007  5211897    0    149      80     4  0.04110  0.07814
    64   0  1.008  5212268    0    139      78     4  0.04264  0.09640
    65   0  1.008  5212812    0    130      76     8  0.04321  0.08747
    66   0  1.009  5213489    0    136      64     3  0.04506  0.11432
    67   0  1.008  5211805    0    126      82     5  0.04191  0.08341
    68   0  1.006  5211547    0    130      85     3  0.04295  0.08776
    69   0  1.005  5211878    0    144      81     4  0.04385  0.09651
    70   0  1.005  5211818    0    135      80     5  0.04112  0.07495
    71   0  1.006  5211404    0    135      86     9  0.04228  0.06499
    72   0  1.010  5211990    0    146      82     6  0.04313  0.08665
    73   0  1.008  5211829    0    149      77     4  0.04692  0.08220
    74   0  1.009  5211361    0    143      89     6  0.03919  0.06360
    75   0  1.009  5210284    0    145      96     4  0.03948  0.06747
    76   0  1.010  5210873    0    142      91     3  0.04197  0.07846
    77   0  1.011  5211766    0    142      86     1  0.04526  0.11282
    78   0  1.009  5211746    0    132      85     5  0.04243  0.08484
    79   0  1.009  5212604    0    127      83     7  0.04001  0.08928
    80   0  1.007  5212350    0    134      79     6  0.04013  0.07192

   470 OUT OF    475 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.53
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.08
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.433    32.145    34.501  88.676 108.481 111.923
 REC. CELL PARAMETERS   0.041665  0.033732  0.030741  83.785  70.580  67.292
 COORDINATES OF UNIT CELL A-AXIS   -14.359   -14.245    18.532
 COORDINATES OF UNIT CELL B-AXIS    14.487    25.537    13.089
 COORDINATES OF UNIT CELL C-AXIS   -21.483    21.535   -16.280
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.013    -0.008     0.031
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001717 -0.000122
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000287  0.005070  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.25   1244.06
 DETECTOR ORIGIN (PIXELS) AT                     1265.87   1237.40
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.22
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     470 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   46     33     55     66     54     30     53     63     59
 SIGMAB (degree)       0.034  0.035  0.035  0.034  0.034  0.035  0.034  0.034  0.033
 SIGMAR (degree)       0.030  0.026  0.029  0.030  0.029  0.026  0.027  0.030  0.029


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    470 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  4  8 10  8  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  1  6 12 16 13  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  2  6  9  7  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  6  8  5  2  1  0  0  0     0  0  0  2  6 11 13  7  4  1  0  0  0     0  0  0  1  4  7  7  4  2  1  0  0  0
     0  0  0  3 15 31 37 27 13  6  2  1  0     0  0  1  6 28 57 61 41 20  8  3  1  0     0  0  1  5 23 38 34 22 10  4  1  0  0
     0  0  1  5 23 49 61 47 23 10  4  1  0     0  0  1 11 43 88100 72 33 15  5  1  0     0  0  1 10 34 60 60 39 18  8  3  0  0
     0  0  0  2 10 26 34 28 14  6  2  1  0     0  0  0  5 20 47 58 45 21  9  3  1  0     0  0  1  5 17 33 37 26 12  5  2  0  0
     0  0  0  1  2  4  7  6  3  1  0  0  0     0  0  0  1  4  9 11  9  5  2  1  0  0     0  0  0  1  4  7  8  6  3  1  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6  8  6  4  2  1  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  9 14 13  8  4  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  8  8  6  3  1  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    640 REFLECTION INTENSITIES INTEGRATED


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
    81   0  1.012  5212061    0    132      85     3  0.04150  0.09140
    82   0  1.009  5211358    0    141      89     3  0.04057  0.08556
    83   0  1.009  5212193    0    140      72     1  0.04386  0.08415
    84   0  1.011  5211810    0    122      80     3  0.04267  0.07752
    85   0  1.011  5212873    0    124      71     2  0.04249  0.10502
    86   0  1.009  5212955    0    113      75     3  0.04596  0.11526
    87   0  1.009  5212520    0    108      75     4  0.04274  0.09370
    88   0  1.010  5212703    0    122      74     3  0.04191  0.09546
    89   0  1.009  5212077    0    128      79     2  0.04287  0.10340
    90   0  1.008  5212099    0    123      79     1  0.04559  0.11853
    91   0  1.009  5211244    0    127      84     2  0.04326  0.10163
    92   0  1.009  5211799    0    119      81     0  0.04302  0.10903
    93   0  1.011  5212492    0    123      78     0  0.04654  0.12850
    94   0  1.009  5212493    0    124      76     0  0.04706  0.11501
    95   0  1.007  5211718    0    128      81     1  0.04321  0.09484
    96   0  1.008  5212057    0    133      79     2  0.04328  0.09119
    97   0  1.008  5212079    0    122      76     3  0.04152  0.07999
    98   0  1.011  5212683    0    124      75     2  0.04532  0.09961
    99   0  1.010  5212512    0    123      74     2  0.04668  0.09834
   100   0  1.008  5210847    0    135      91     2  0.04562  0.10244

   439 OUT OF    440 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.45
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.07
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.426    32.128    34.502  88.659 108.478 111.906
 REC. CELL PARAMETERS   0.041666  0.033744  0.030739  83.812  70.592  67.316
 COORDINATES OF UNIT CELL A-AXIS   -14.355   -14.229    18.538
 COORDINATES OF UNIT CELL B-AXIS    14.486    25.523    13.075
 COORDINATES OF UNIT CELL C-AXIS   -21.479    21.544   -16.274
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.002    -0.003     0.013
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001717 -0.000122
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000734  0.005109  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.30   1244.08
 DETECTOR ORIGIN (PIXELS) AT                     1265.34   1237.37
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.17
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     439 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     27     53     62     54     27     54     60     51
 SIGMAB (degree)       0.034  0.034  0.034  0.034  0.033  0.034  0.034  0.034  0.033
 SIGMAR (degree)       0.032  0.028  0.029  0.031  0.030  0.026  0.028  0.030  0.031


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    439 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  1  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  1  3  7  9  8  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  1  0  0  0     0  0  0  1  5 11 15 13  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  0  3  6  8  7  4  2  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  1  1  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  5  6  4  2  1  0  0  0     0  0  0  1  5  9 11  6  3  1  0  0  0     0  0  0  1  4  6  6  3  1  1  0  0  0
     0  0  0  3 12 27 33 25 13  5  2  1  0     0  0  0  7 28 53 56 38 19  7  3  1  0     0  0  1  6 20 36 33 21 10  4  1  0  0
     0  0  0  4 21 47 58 43 23  9  3  1  0     0  0  1 10 46 93100 67 33 13  5  2  0     0  0  1  9 34 61 60 38 18  7  3  1  0
     0  0  0  2 11 24 33 24 12  5  2  0  0     0  0  0  5 23 47 56 40 20  8  3  1  0     0  0  0  4 18 32 35 23 11  4  2  0  0
     0  0  0  0  2  4  6  5  3  1  0  0  0     0  0  0  1  4  7 10  8  4  2  1  0  0     0  0  0  1  3  5  6  5  3  1  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  8  8  5  2  1  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  9 14 14  8  4  2  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  7  8  5  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    614 REFLECTION INTENSITIES INTEGRATED


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
   101   0  1.008  5210974    0    142      87     2  0.03787  0.08099
   102   0  1.006  5211039    0    143      90     3  0.03980  0.08645
   103   0  1.006  5211826    0    139      84     3  0.04229  0.09207
   104   0  1.002  5210971    0    138      93     1  0.04098  0.10394
   105   0  1.005  5211529    0    147      87     1  0.04349  0.10165
   106   0  1.004  5211356    0    147      87     5  0.04519  0.08122
   107   0  1.004  5211268    0    140      84     3  0.04186  0.08901
   108   0  1.003  5212576    0    140      79     2  0.04326  0.10771
   109   0  1.005  5212127    0    145      81     3  0.04563  0.09861
   110   0  1.003  5210421    0    144      96     7  0.04124  0.07347
   111   0  1.002  5210521    0    158      97    10  0.04110  0.07731
   112   0  1.004  5211036    0    154      95     3  0.04267  0.09300
   113   0  1.001  5211067    0    137      93     2  0.04508  0.10384
   114   0  1.001  5212182    0    120      82     5  0.04337  0.09656
   115   0  1.000  5213435    0    115      67     2  0.04556  0.11786
   116   0  1.002  5213220    0    110      70     0  0.04548  0.11527
   117   0  1.002  5212912    0    126      76     0  0.04483  0.13209
   118   0  0.997  5211595    0    126      88     0  0.04269  0.12471
   119   0  0.998  5211313    0    131      87     3  0.04111  0.09310
   120   0  1.000  5211033    0    129      96     3  0.04343  0.10099

   493 OUT OF    494 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.57
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.432    32.133    34.511  88.659 108.479 111.917
 REC. CELL PARAMETERS   0.041661  0.033742  0.030731  83.809  70.590  67.305
 COORDINATES OF UNIT CELL A-AXIS   -14.358   -14.237    18.538
 COORDINATES OF UNIT CELL B-AXIS    14.491    25.525    13.078
 COORDINATES OF UNIT CELL C-AXIS   -21.485    21.553   -16.275
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.004    -0.006     0.018
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001717 -0.000122
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000644  0.004990  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.32   1244.06
 DETECTOR ORIGIN (PIXELS) AT                     1265.47   1237.51
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.23
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     493 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   49     30     57     68     54     32     58     74     62
 SIGMAB (degree)       0.034  0.035  0.035  0.034  0.033  0.035  0.034  0.034  0.033
 SIGMAR (degree)       0.030  0.027  0.029  0.031  0.030  0.026  0.028  0.031  0.031


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    493 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  1  1  1  0  1  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  1  4 10 13 10  5  2  1  0  0
     0  0  0  0  0  1  0  1  1  0  0  0  0     0  0  0  0  1  2  3  3  1  1  0  0  0     0  0  0  1  6 15 20 16  8  4  1  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  1  3  8 11  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  7  8  6  3  1  0  0  0     0  0  0  2  6 12 11  7  3  1  0  0  0     0  0  0  1  4  6  6  4  1  0  0  0  0
     0  0  0  3 15 36 44 29 15  7  2  1  0     0  0  1  7 30 59 64 39 19  8  3  1  0     0  0  1  5 20 34 35 20  9  4  1  0  0
     0  0  0  6 26 54 68 51 25 11  4  1  0     0  0  1 11 50 92100 68 33 13  5  2  1     0  0  1  9 31 53 54 35 16  6  2  1  0
     0  0  0  3 13 30 40 30 15  6  2  1  0     0  0  1  6 24 51 60 43 21  9  3  1  0     0  0  1  4 15 30 33 23 11  4  2  0  0
     0  0  0  1  2  6  8  6  4  2  0  0  0     0  0  0  1  5 11 14  9  5  2  1  0  0     0  0  0  1  3  7  8  5  3  1  0  0  0
     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  7  4  2  1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  2  6 11 11  7  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  7  5  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    640 REFLECTION INTENSITIES INTEGRATED


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
   121   0  1.000  5211975    0    128      82     3  0.04563  0.09614
   122   0  1.001  5211925    0    127      85     2  0.04551  0.12129
   123   0  1.000  5211514    0    124      88     3  0.04591  0.12048
   124   0  1.000  5211551    0    125      82     4  0.04588  0.10912
   125   0  1.000  5211976    0    124      87     4  0.04651  0.12614
   126   0  0.999  5210825    0    130      88     2  0.04684  0.11377
   127   0  0.998  5209890    0    141      98     2  0.04448  0.11500
   128   0  0.994  5210535    0    143      96     1  0.04463  0.11987
   129   0  0.997  5211069    0    131      89     1  0.04288  0.11731
   130   0  0.997  5211605    0    133      87     1  0.04431  0.13677
   131   0  0.998  5210910    0    132      89     0  0.04561  0.11737
   132   0  0.998  5210675    0    129      94     2  0.04524  0.12355
   133   0  0.997  5210806    0    136      94     1  0.04761  0.11159
   134   0  0.996  5210481    0    139      98     2  0.04625  0.13129
   135   0  0.997  5210322    0    145     103     2  0.04348  0.13459
   136   0  0.994  5210501    0    147      97     3  0.04519  0.11438
   137   0  0.994  5209925    0    158     101     4  0.04330  0.10459
   138   0  0.992  5209185    0    159     106     3  0.04479  0.09340
   139   0  0.993  5209206    0    157     105     4  0.04327  0.09791
   140   0  0.993  5209791    0    149     104     3  0.04692  0.11796

   492 OUT OF    497 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.55
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.433    32.125    34.510  88.657 108.476 111.905
 REC. CELL PARAMETERS   0.041655  0.033747  0.030730  83.816  70.596  67.318
 COORDINATES OF UNIT CELL A-AXIS   -14.359   -14.228    18.546
 COORDINATES OF UNIT CELL B-AXIS    14.490    25.520    13.068
 COORDINATES OF UNIT CELL C-AXIS   -21.484    21.551   -16.277
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.001    -0.001     0.007
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001717 -0.000122
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000801  0.005076  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.36   1244.07
 DETECTOR ORIGIN (PIXELS) AT                     1265.30   1237.40
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.23
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     492 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   52     36     57     62     49     27     57     72     64
 SIGMAB (degree)       0.034  0.035  0.035  0.034  0.033  0.034  0.034  0.034  0.034
 SIGMAR (degree)       0.030  0.028  0.029  0.031  0.031  0.029  0.028  0.030  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    492 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  1  4  9 13 10  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  3  2  1  0  0  0     0  0  0  1  6 15 21 18 10  4  2  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  2  1  0  0  0  0     0  0  0  1  4  8 11 10  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  3  6  7  5  2  1  0  0  0     0  0  0  2  6 12 10  7  3  1  0  0  0     0  0  0  1  4  7  6  4  1  0  0  0  0
     0  0  0  3 18 38 44 30 14  6  2  1  0     0  0  1  6 31 60 62 38 17  7  2  1  0     0  0  1  4 17 32 33 19  9  3  1  0  0
     0  0  0  6 27 60 71 54 28 11  4  1  0     0  0  1 11 48 95100 68 34 13  5  1  0     0  0  1  6 27 50 52 35 17  7  2  1  0
     0  0  0  3 14 32 41 31 17  7  3  1  0     0  0  0  6 26 52 60 41 22  9  3  1  0     0  0  0  4 15 29 32 21 11  5  2  0  0
     0  0  0  1  3  5  8  8  4  2  0  0  0     0  0  0  1  5  9 13 11  5  2  1  0  0     0  0  0  1  3  5  7  6  2  1  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  7  4  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0
     0  0  0  2  6 10 11  7  4  1  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  1  3  6  7  5  3  1  0  0 -1     0  0  0  1  1  1  1  1  1  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    648 REFLECTION INTENSITIES INTEGRATED


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
   141   0  0.991  5210334    0    148      97     1  0.04610  0.12599
   142   0  0.995  5209936    0    153     102     1  0.04529  0.12731
   143   0  0.993  5210190    0    147     103     3  0.04452  0.10975
   144   0  0.992  5210672    0    142      92     4  0.04288  0.09822
   145   0  0.989  5210399    0    130      94     3  0.04236  0.10257
   146   0  0.991  5210481    0    136      90     6  0.04249  0.08574
   147   0  0.991  5209707    0    138      99     5  0.04765  0.10316
   148   0  0.992  5209368    0    147     104     5  0.04378  0.10568
   149   0  0.990  5208988    0    143     112     5  0.04568  0.10349
   150   0  0.992  5210472    0    128     101     8  0.04495  0.11072
   151   0  0.990  5211977    0    127      85     3  0.04971  0.14968
   152   0  0.991  5212047    0    119      81     3  0.04633  0.14355
   153   0  0.989  5212111    0    126      79     8  0.04692  0.09335
   154   0  0.988  5211394    0    134      84     4  0.04582  0.10468
   155   0  0.989  5210349    0    139      95     2  0.04370  0.10888
   156   0  0.988  5209783    0    146     101     1  0.04266  0.11207
   157   0  0.987  5209573    0    154     107     3  0.04298  0.10213
   158   0  0.988  5210013    0    150     103     3  0.04429  0.10360
   159   0  0.986  5210036    0    143     104     1  0.04443  0.12874
   160   0  0.985  5210147    0    133     101     2  0.04406  0.11355

   494 OUT OF    496 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.60
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.15
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.420    32.123    34.498  88.639 108.467 111.870
 REC. CELL PARAMETERS   0.041659  0.033738  0.030738  83.854  70.617  67.362
 COORDINATES OF UNIT CELL A-AXIS   -14.352   -14.207    18.548
 COORDINATES OF UNIT CELL B-AXIS    14.485    25.522    13.066
 COORDINATES OF UNIT CELL C-AXIS   -21.476    21.548   -16.266
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.001     0.004    -0.013
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001717 -0.000122
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.001087  0.005248  1.023527
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.38   1244.10
 DETECTOR ORIGIN (PIXELS) AT                     1264.95   1237.21
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.12
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     494 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   51     36     60     70     62     33     52     60     53
 SIGMAB (degree)       0.034  0.034  0.035  0.034  0.033  0.034  0.034  0.034  0.033
 SIGMAR (degree)       0.033  0.031  0.030  0.030  0.029  0.027  0.029  0.031  0.031


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    494 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  0  1  1  1  0  0  0  0     0  0  0  0  2  5  7  6  4  1  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  1  0  0  0     0  0  0  1  3  8 12 11  6  3  1  0  0
     0  0  0  0  0  0  1  1  1  0  0  0  0     0  0  0  0  1  0  1  1  1  0  0  0  0     0  0  0  0  2  5  7  7  4  2  1  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  2  5  5  4  2  1  0  0  0     0  0  0  1  5 10 10  6  3  1  0  0  0     0  0  0  1  4  9  9  5  2  1  0  0  0
     0  0  0  2 10 24 31 23 12  5  2  0  0     0  0  1  6 28 54 59 39 18  7  3  1  0     0  0  1  7 29 48 47 29 13  5  2  0  0
     0  0  0  3 17 43 52 39 21  9  3  1  0     0  0  1  9 46 94100 68 34 14  5  1  0     0  0  1 10 44 79 80 52 26 10  4  1  0
     0  0  0  2  9 23 30 23 12  5  2  1  0     0  0  0  6 25 51 57 39 20  8  3  1  0     0  0  1  7 24 45 46 31 15  6  2  1  0
     0  0  0  0  1  3  6  5  2  1  0  0  0     0  0  0  1  4  8 11  9  4  2  0  0  0     0  0  0  1  4  7  9  7  3  1  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  4  3  2  1  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  9 15 15 10  4  2  1  0  0     0  0  0  0  2  2  2  1  1  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  4 14 25 26 16  8  3  1  0  0     0  0  0  1  2  3  4  3  1  1  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  2  8 15 15 11  5  2  1  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0
     0  0  0  1  2  3  3  3  1  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
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
   161   0  0.990  5210537    0    146      97     3  0.04563  0.11213
   162   0  0.989  5210380    0    142      95     5  0.04446  0.10103
   163   0  0.988  5209306    0    140     107     4  0.04280  0.10263
   164   0  0.986  5209825    0    141     101     4  0.04406  0.10343
   165   0  0.989  5210830    0    144      90     5  0.04581  0.10980
   166   0  0.989  5211758    0    135      80     2  0.04438  0.13284
   167   0  0.988  5211134    0    127      87     1  0.04561  0.13090
   168   0  0.989  5211155    0    133      85     5  0.04613  0.10003
   169   0  0.986  5211342    0    131      86     4  0.04592  0.09481
   170   0  0.986  5210800    0    142      90     4  0.04303  0.09396
   171   0  0.988  5210658    0    143      96     6  0.04426  0.07885
   172   0  0.989  5209719    0    147     103     7  0.04180  0.07709
   173   0  0.987  5210358    0    139      89     2  0.04148  0.10451
   174   0  0.985  5210616    0    134      87     1  0.04475  0.12024
   175   0  0.985  5211497    0    140      81     2  0.04500  0.10878
   176   0  0.987  5211923    0    130      74     6  0.04355  0.08150
   177   0  0.985  5211318    0    128      85     4  0.04572  0.09043
   178   0  0.984  5210962    0    129      91     2  0.04336  0.13590
   179   0  0.987  5211471    0    129      87     0  0.05271  0.16499
   180   0  0.987  5211099    0    137      93     5  0.04763  0.11072

   482 OUT OF    483 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.55
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.14
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.433    32.135    34.513  88.647 108.478 111.914
 REC. CELL PARAMETERS   0.041657  0.033738  0.030728  83.823  70.597  67.311
 COORDINATES OF UNIT CELL A-AXIS   -14.359   -14.236    18.541
 COORDINATES OF UNIT CELL B-AXIS    14.494    25.525    13.078
 COORDINATES OF UNIT CELL C-AXIS   -21.485    21.561   -16.269
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.012    -0.006     0.014
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001717 -0.000122
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000747  0.005051  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.35   1244.02
 DETECTOR ORIGIN (PIXELS) AT                     1265.37   1237.38
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.26
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     482 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   48     35     57     61     50     29     53     72     63
 SIGMAB (degree)       0.034  0.034  0.035  0.034  0.033  0.035  0.034  0.034  0.034
 SIGMAR (degree)       0.030  0.028  0.029  0.030  0.029  0.025  0.028  0.031  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    482 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  1  5 10 13  9  5  2  1  0  0
     0  0  0  0  0  0  1  0  0  0  0  0 -1     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  1  7 15 20 16  9  4  1  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  1  3  7 10 10  6  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  8  4  2  1  0  0  0     0  0  0  2  7 12 11  7  3  1  0  0  0     0  0  0  1  3  6  6  3  2  0  0  0  0
     0  0  1  4 18 37 42 29 14  6  2  1  0     0  0  1  7 31 62 62 39 18  8  3  1  0     0  0  1  5 18 34 31 19  9  3  1  0  0
     0  0  1  5 27 55 66 50 26 11  4  1  0     0  0  1 11 49 96100 70 33 14  5  2  0     0  0  1  8 30 55 52 34 16  7  3  1  0
     0  0  0  3 12 27 37 31 17  7  3  1  0     0  0  1  5 24 49 59 42 22  9  3  1  0     0  0  0  4 16 29 33 22 11  5  2  0  0
     0  0  0  1  2  6  8  7  3  1  0  0  0     0  0  0  1  4 11 13 10  4  2  1  0  0     0  0  0  1  3  6  8  5  2  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  6  6  4  2  1  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 11 11  6  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  6  7  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    634 REFLECTION INTENSITIES INTEGRATED


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
   181   0  0.988  5209893    0    130     103     8  0.04524  0.08849
   182   0  0.988  5210752    0    149      94     2  0.04895  0.14388
   183   0  0.986  5210877    0    141      91     6  0.04383  0.11159
   184   0  0.984  5210358    0    137      96     4  0.04308  0.09759
   185   0  0.983  5211219    0    138      89     5  0.04737  0.10063
   186   0  0.986  5210654    0    128      97     2  0.04415  0.12347
   187   0  0.985  5210717    0    124      92     3  0.04530  0.11743
   188   0  0.983  5211553    0    124      86     2  0.04757  0.14318
   189   0  0.983  5211526    0    123      79     2  0.04442  0.12161
   190   0  0.984  5211958    0    120      80     6  0.04415  0.09122
   191   0  0.984  5212385    0    130      76     2  0.04479  0.09254
   192   0  0.984  5211990    0    130      84     3  0.03870  0.08011
   193   0  0.982  5211901    0    125      78     2  0.03976  0.08006
   194   0  0.981  5212453    0    126      78     1  0.04556  0.11707
   195   0  0.983  5212465    0    121      78     3  0.04454  0.12799
   196   0  0.985  5211667    0    123      77     7  0.04368  0.09302
   197   0  0.983  5210740    0    132      90     7  0.04517  0.11771
   198   0  0.983  5210808    0    137      87     3  0.04255  0.11473
   199   0  0.981  5211204    0    142      82     3  0.04389  0.10849
   200   0  0.980  5211218    0    136      86     2  0.04434  0.12489

   473 OUT OF    474 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.61
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.17
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.428    32.135    34.506  88.680 108.487 111.919
 REC. CELL PARAMETERS   0.041673  0.033742  0.030738  83.780  70.572  67.294
 COORDINATES OF UNIT CELL A-AXIS   -14.355   -14.234    18.537
 COORDINATES OF UNIT CELL B-AXIS    14.492    25.528    13.075
 COORDINATES OF UNIT CELL C-AXIS   -21.483    21.538   -16.287
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.009    -0.002     0.012
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001717 -0.000122
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000422  0.005184  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.25   1244.07
 DETECTOR ORIGIN (PIXELS) AT                     1265.69   1237.26
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.21
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     473 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   48     29     56     71     55     30     52     59     56
 SIGMAB (degree)       0.034  0.035  0.035  0.034  0.033  0.034  0.034  0.034  0.034
 SIGMAR (degree)       0.031  0.028  0.030  0.031  0.030  0.026  0.028  0.030  0.029


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    473 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  1  5 11 13 10  5  2  1  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0     0  0  0  2  7 16 20 17  9  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  3  8 12  9  6  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  8 10  6  3  1  0  0  0     0  0  0  2  5 12 13  7  3  1  0  0  0     0  0  0  1  4  7  7  4  2  0  0  0  0
     0  0  0  4 17 37 41 30 15  6  2  1  0     0  0  1  8 29 60 60 38 19  8  3  1  0     0  0  1  6 20 35 32 18  8  4  1  0  0
     0  0  1  6 25 55 66 52 26 11  4  1  0     0  0  1 10 46 91100 71 34 14  5  2  0     0  0  1  8 32 54 53 35 17  7  2  1  0
     0  0  0  3 12 29 40 30 17  7  2  1  0     0  0  1  6 23 51 62 45 22  9  3  1  0     0  0  1  5 18 32 36 24 12  4  2  0  0
     0  0  0  1  2  5  8  7  3  1  1  0  0     0  0  0  1  4 10 13 10  5  2  1  0  0     0  0  0  1  4  7  9  5  2  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  2  5  7  6  3  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  1  0  0  0  0
     0  0  0  2  8 12 11  7  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  1  4  7  8  5  2  1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    624 REFLECTION INTENSITIES INTEGRATED


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
   201   0  0.983  5211271    0    129      85     4  0.04174  0.08475
   202   0  0.986  5211582    0    137      84     2  0.04565  0.11617
   203   0  0.987  5211880    0    135      82     3  0.04444  0.10057
   204   0  0.984  5210943    0    134      88     3  0.04413  0.08364
   205   0  0.985  5211507    0    128      77     2  0.04132  0.09137
   206   0  0.979  5211932    0    133      73     6  0.04364  0.07344
   207   0  0.982  5212058    0    138      73     3  0.04043  0.08285
   208   0  0.984  5210107    0    141      96     5  0.04046  0.08131
   209   0  0.983  5209852    0    152     107     1  0.04252  0.11283
   210   0  0.983  5210338    0    158      99     1  0.04714  0.11926
   211   0  0.982  5208109    0    160     120     3  0.04496  0.11021
   212   0  0.981  5207772    0    166     123     4  0.04307  0.11066
   213   0  0.986  5207603    0    171     123     8  0.04373  0.09062
   214   0  0.987  5208824    0    158     109     6  0.04534  0.09173
   215   0  0.981  5210164    0    151      97     7  0.04471  0.08707
   216   0  0.982  5211030    0    152      89     7  0.03943  0.07267
   217   0  0.983  5210766    0    140      94     0  0.04299  0.11949
   218   0  0.980  5210725    0    143      89     3  0.04189  0.09720
   219   0  0.981  5211115    0    135      84     1  0.04245  0.11677
   220   0  0.980  5210973    0    133      91     2  0.04596  0.12199

   519 OUT OF    521 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.52
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.430    32.137    34.510  88.663 108.483 111.923
 REC. CELL PARAMETERS   0.041668  0.033740  0.030733  83.799  70.583  67.296
 COORDINATES OF UNIT CELL A-AXIS   -14.357   -14.234    18.538
 COORDINATES OF UNIT CELL B-AXIS    14.495    25.530    13.074
 COORDINATES OF UNIT CELL C-AXIS   -21.484    21.549   -16.281
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.001    -0.002     0.010
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001717 -0.000122
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000495  0.004986  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.30   1244.09
 DETECTOR ORIGIN (PIXELS) AT                     1265.65   1237.54
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.24
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     519 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   55     36     61     68     64     33     59     73     60
 SIGMAB (degree)       0.034  0.035  0.035  0.034  0.033  0.034  0.034  0.034  0.034
 SIGMAR (degree)       0.030  0.027  0.029  0.031  0.030  0.027  0.028  0.030  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    519 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  1  0  0  0     0  0  0  1  5 10 13 11  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  2  8 16 21 17  9  4  2  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  4  9 12  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  1  4  8  8  5  2  1  0  0  0     0  0  0  1  6 12 10  7  3  1  0  0  0     0  0  0  1  4  7  5  4  1  0  0  0  0
     0  0  0  4 18 37 46 32 15  6  2  1  0     0  0  1  7 31 59 63 40 19  7  3  1  0     0  0  1  6 21 35 34 20 10  3  1  0  0
     0  0  1  7 29 59 72 54 27 11  4  1  0     0  0  1 12 50 92100 69 34 13  5  1  0     0  0  1 10 34 55 54 36 17  7  2  1  0
     0  0  0  3 16 32 41 30 16  7  2  1  0     0  0  1  6 27 50 58 40 20  8  3  1  0     0  0  1  5 19 30 33 21 11  4  2  0  0
     0  0  0  1  3  6  8  6  4  1  0  0  0     0  0  0  1  6 10 12  8  5  2  1  0  0     0  0  0  1  4  6  7  4  2  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  2  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  6  3  2  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 12 11  7  3  1  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  7  5  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    664 REFLECTION INTENSITIES INTEGRATED


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
   221   0  0.982  5210977    0    126      85     2  0.04687  0.11884
   222   0  0.986  5211360    0    129      77     2  0.04566  0.11781
   223   0  0.984  5210921    0    143      86     6  0.04515  0.08449
   224   0  0.984  5210639    0    145      81     6  0.04196  0.07320
   225   0  0.984  5209915    0    140      90     4  0.04271  0.08256
   226   0  0.981  5210535    0    138      94     2  0.04399  0.09831
   227   0  0.984  5211266    0    144      80     0  0.04608  0.11631
   228   0  0.981  5211016    0    137      83     0  0.04434  0.11660
   229   0  0.980  5210685    0    130      87     3  0.04071  0.09353
   230   0  0.980  5211602    0    126      80     5  0.04476  0.08979
   231   0  0.981  5212257    0    131      71     3  0.04335  0.10358
   232   0  0.983  5210799    0    122      88     6  0.04390  0.07987
   233   0  0.984  5210960    0    124      92     8  0.04523  0.09631
   234   0  0.983  5211729    0    129      81     4  0.04511  0.12938
   235   0  0.984  5212284    0    119      76     4  0.04655  0.12654
   236   0  0.981  5211943    0    130      78     5  0.04533  0.10199
   237   0  0.981  5211409    0    145      81     6  0.04317  0.07494
   238   0  0.980  5210054    0    151      92     5  0.04575  0.08559
   239   0  0.981  5209965    0    150      93     3  0.04202  0.09542
   240   0  0.981  5209420    0    128      95     6  0.04199  0.07790

   476 OUT OF    476 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.58
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.435    32.129    34.518  88.649 108.482 111.929
 REC. CELL PARAMETERS   0.041661  0.033748  0.030724  83.814  70.590  67.295
 COORDINATES OF UNIT CELL A-AXIS   -14.359   -14.236    18.542
 COORDINATES OF UNIT CELL B-AXIS    14.491    25.526    13.066
 COORDINATES OF UNIT CELL C-AXIS   -21.488    21.558   -16.280
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.000    -0.002     0.011
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001717 -0.000122
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000633  0.005043  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.34   1244.09
 DETECTOR ORIGIN (PIXELS) AT                     1265.50   1237.46
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.28
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     476 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   48     30     60     69     50     28     55     63     53
 SIGMAB (degree)       0.034  0.035  0.035  0.034  0.033  0.035  0.034  0.034  0.034
 SIGMAR (degree)       0.030  0.027  0.028  0.031  0.030  0.025  0.027  0.030  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    476 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  3  4  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  1  7 14 18 15  8  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  3  4  4  2  1  0  0  0     0  0  0  2  9 21 28 24 13  6  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  1  0  0  0     0  0  0  1  4 11 15 15  7  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  3  3  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  8 11  5  3  1  0  0  0     0  0  0  2  7 11 13  6  3  1  0  0  0     0  0  0  1  3  4  5  2  1  0  0  0  0
     0  0  0  5 24 47 53 39 19  8  2  1  0     0  0  1  9 36 61 62 42 19  7  2  1  0     0  0  0  4 16 25 24 16  7  3  1  0  0
     0  0  1  7 34 73 84 62 32 14  5  1  1     0  0  1 12 52 97100 70 34 14  5  1  1     0  0  1  6 24 42 42 29 13  5  2  0  0
     0  0  0  4 16 39 48 40 20  9  3  1  0     0  0  1  7 24 52 59 45 22 10  3  1  0     0  0  0  4 12 23 24 18  9  4  1  0  0
     0  0  0  1  3  8 11  9  5  2  1  0  0     0  0  0  1  5 11 15 11  6  2  1  0  0     0  0  0  1  2  5  6  5  2  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0 -1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0 -1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  4  4  3  1  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  7  5  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0
     0  0  0  1  2  4  4  3  2  1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     1  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    622 REFLECTION INTENSITIES INTEGRATED


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
   241   0  0.981  5210606    0    129      88     2  0.04368  0.10904
   242   0  0.981  5211426    0    132      83     5  0.04448  0.10170
   243   0  0.983  5211183    0    126      85     2  0.04782  0.12310
   244   0  0.980  5210478    0    136      98     4  0.04649  0.10647
   245   0  0.979  5211147    0    142      89     1  0.04884  0.12006
   246   0  0.979  5211026    0    137      84     7  0.04220  0.06920
   247   0  0.982  5210769    0    135      81     4  0.04217  0.07501
   248   0  0.980  5210943    0    129      81     5  0.04308  0.08177
   249   0  0.976  5211314    0    136      82     1  0.04501  0.13397
   250   0  0.977  5210813    0    137      94     2  0.04520  0.12696
   251   0  0.975  5209362    0    144     108     4  0.04497  0.10935
   252   0  0.974  5209389    0    148      97     4  0.04584  0.10024
   253   0  0.973  5210005    0    138      97     5  0.04309  0.10810
   254   0  0.977  5211669    0    142      81     5  0.04376  0.11156
   255   0  0.978  5212002    0    131      81     6  0.04508  0.09632
   256   0  0.977  5211783    0    129      87     3  0.04038  0.13294
   257   0  0.977  5211796    0    130      88     3  0.04208  0.12719
   258   0  0.977  5212670    0    129      74     1  0.04504  0.10211
   259   0  0.977  5212781    0    125      71     2  0.04629  0.12424
   260   0  0.978  5211857    0    121      82     1  0.04219  0.13237

   476 OUT OF    478 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.53
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.430    32.134    34.515  88.631 108.474 111.915
 REC. CELL PARAMETERS   0.041659  0.033739  0.030724  83.843  70.608  67.317
 COORDINATES OF UNIT CELL A-AXIS   -14.355   -14.231    18.542
 COORDINATES OF UNIT CELL B-AXIS    14.494    25.528    13.071
 COORDINATES OF UNIT CELL C-AXIS   -21.487    21.566   -16.263
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.013    -0.009     0.004
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001717 -0.000122
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000930  0.005037  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.32   1243.97
 DETECTOR ORIGIN (PIXELS) AT                     1265.09   1237.36
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.26
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     476 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   51     26     47     67     57     30     55     72     59
 SIGMAB (degree)       0.034  0.035  0.035  0.034  0.033  0.034  0.034  0.034  0.034
 SIGMAR (degree)       0.031  0.028  0.029  0.031  0.030  0.026  0.028  0.030  0.031


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    476 STRONG REFLECTIONS *****


     0  1  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     1  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  1  4 10 12 10  5  2  1  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  2  2  2  2  1  0  0  0     0  0  0  1  7 15 19 17 10  4  2  0  0
    -1  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  1  2  1  1  1  0  0  0     0  0  0  1  3  8 11 10  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  7  5  2  1  0  0  0     0  0  0  1  7 11 11  6  3  1  0  0  0     0  0  0  1  4  6  6  3  2  0  0  0  0
     0  0  0  4 16 36 40 29 15  5  2  0  0     0  0  1  8 30 59 57 37 18  6  2  0  0     0  0  1  5 17 32 30 19  9  3  1  0  0
     0  0  0  6 26 57 68 51 28 11  5  1  0     0  0  1 12 49 96100 68 35 13  5  1  0     0  0  1  9 30 54 53 35 18  7  2  1  0
     0  0  0  3 14 31 40 32 17  7  3  1  0     0  0  0  6 26 53 61 45 22  9  3  1  0     0  0  0  4 17 31 34 24 11  5  2  0  0
     0  0  0  0  2  5  8  6  4  2  0  0  0     0  0  0  1  4  9 12  8  5  2  1  0  0     0  0  0  1  3  5  7  5  2  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  6  4  2  1  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 11 11  7  4  1  0  0  0     0  0  0  1  1  1  2  1  1  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  4  7  8  5  3  1  0  0  0     0  0  0  0  1  1  2  0  1  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  1  0  0  0  0

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
   261   0  0.977  5212023    0    120      84     4  0.04536  0.10106
   262   0  0.976  5212337    0    121      79     1  0.04268  0.11502
   263   0  0.975  5211990    0    124      75     5  0.04484  0.07553
   264   0  0.977  5211682    0    123      87     3  0.04345  0.09394
   265   0  0.975  5211759    0    133      84     6  0.04079  0.09231
   266   0  0.976  5211158    0    137      87     4  0.04593  0.10035
   267   0  0.979  5210284    0    142      98     2  0.04164  0.09708
   268   0  0.980  5211028    0    146      94     2  0.04347  0.10048
   269   0  0.976  5211709    0    141      87     5  0.04229  0.08850
   270   0  0.977  5211040    0    122      85     5  0.04224  0.07248
   271   0  0.977  5211722    0    120      83     8  0.04173  0.06812
   272   0  0.978  5212425    0    121      75     4  0.04486  0.08800
   273   0  0.979  5211538    0    135      87     5  0.04535  0.09813
   274   0  0.979  5211108    0    126      92     2  0.04528  0.13099
   275   0  0.978  5210266    0    134      99     4  0.04331  0.11285
   276   0  0.975  5210358    0    137     101     5  0.04556  0.10107
   277   0  0.976  5210526    0    138      95     2  0.04787  0.11737
   278   0  0.971  5210286    0    146      97     1  0.04170  0.10783
   279   0  0.977  5209912    0    150     104     2  0.04462  0.10338
   280   0  0.979  5210895    0    158      95     2  0.04629  0.10425

   474 OUT OF    477 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.58
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.426    32.138    34.511  88.661 108.488 111.926
 REC. CELL PARAMETERS   0.041676  0.033739  0.030733  83.798  70.578  67.293
 COORDINATES OF UNIT CELL A-AXIS   -14.355   -14.235    18.533
 COORDINATES OF UNIT CELL B-AXIS    14.491    25.532    13.076
 COORDINATES OF UNIT CELL C-AXIS   -21.483    21.548   -16.285
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.006     0.001     0.014
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001717 -0.000122
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000468  0.005130  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.25   1244.03
 DETECTOR ORIGIN (PIXELS) AT                     1265.63   1237.29
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.24
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     474 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   42     35     68     67     51     30     50     64     56
 SIGMAB (degree)       0.035  0.035  0.035  0.034  0.033  0.035  0.035  0.035  0.034
 SIGMAR (degree)       0.031  0.027  0.028  0.031  0.030  0.027  0.029  0.031  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    474 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0 -1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  3  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  1  5 11 13 11  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  2  7 15 20 16  9  4  1  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  3  7 11  9  6  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  3  3  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  4  9 10  6  2  1  0  0  0     0  0  0  2  7 13 13  7  3  1  0  0  0     0  0  0  1  5  7  7  3  1  1  0  0  0
     0  0  1  4 17 37 45 33 16  7  2  1  0     0  0  1  7 30 61 67 43 20  8  3  1  0     0  0  1  6 20 33 34 21  9  4  1  0  0
     0  0  1  6 27 54 67 51 27 12  5  1  0     0  0  1 11 49 91100 69 35 15  6  2  0     0  0  1  9 33 53 51 35 17  7  3  1  0
     0  0  0  3 13 28 38 29 17  7  3  1  0     0  0  1  6 25 49 60 42 23  9  4  1  0     0  0  1  5 17 29 31 22 11  4  2  0  0
     0  0  0  1  2  6  9  8  4  2  1  0  0     0  0  0  1  5 10 13 11  6  2  1  0  0     0  0  0  1  3  6  8  6  3  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  7  7  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 11 10  7  3  1  1  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  6  5  2  1  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
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
   281   0  0.975  5209937    0    154     102     7  0.04260  0.09322
   282   0  0.982  5209785    0    153     104     5  0.04327  0.09211
   283   0  0.977  5210208    0    152      98     6  0.04386  0.08302
   284   0  0.977  5211460    0    137      87     2  0.04408  0.11632
   285   0  0.975  5211862    0    132      84     2  0.04550  0.13838
   286   0  0.974  5211074    0    127      92     3  0.04790  0.15087
   287   0  0.974  5210737    0    125      88     4  0.04556  0.12789
   288   0  0.973  5211182    0    132      86     5  0.04617  0.12249
   289   0  0.974  5211065    0    137      81     6  0.04871  0.09574
   290   0  0.973  5210191    0    137      94     6  0.04613  0.09506
   291   0  0.975  5209960    0    145      97     6  0.04391  0.09886
   292   0  0.973  5210415    0    128      99     2  0.04512  0.14386
   293   0  0.974  5211292    0    125      87     1  0.04694  0.14874
   294   0  0.974  5211075    0    120      93     1  0.04646  0.15807
   295   0  0.975  5211643    0    123      82     0  0.04913  0.15877
   296   0  0.974  5211226    0    131      89     2  0.04994  0.14539
   297   0  0.972  5209707    0    125      99     4  0.04582  0.11244
   298   0  0.970  5209251    0    135     102     5  0.04642  0.09904
   299   0  0.970  5209777    0    137      98     3  0.04822  0.12174
   300   0  0.968  5210026    0    129      94     3  0.05004  0.12225

   469 OUT OF    472 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.58
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.13
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.426    32.142    34.509  88.666 108.482 111.897
 REC. CELL PARAMETERS   0.041666  0.033728  0.030733  83.807  70.586  67.322
 COORDINATES OF UNIT CELL A-AXIS   -14.354   -14.226    18.541
 COORDINATES OF UNIT CELL B-AXIS    14.492    25.534    13.080
 COORDINATES OF UNIT CELL C-AXIS   -21.484    21.546   -16.283
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.005    -0.000     0.004
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001717 -0.000122
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000467  0.005126  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.34   1244.11
 DETECTOR ORIGIN (PIXELS) AT                     1265.72   1237.38
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.23
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     469 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   46     31     52     66     54     30     53     67     59
 SIGMAB (degree)       0.035  0.035  0.035  0.034  0.033  0.035  0.035  0.035  0.034
 SIGMAR (degree)       0.031  0.028  0.029  0.031  0.030  0.027  0.029  0.031  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    469 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  1  4  8  9  8  4  2  1  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  1  6 12 15 13  7  3  1  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  1  3  6  9  7  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  4  8  7  4  2  1  0  0  0     0  0  0  2  7 13 12  7  3  1  0  0  0     0  0  0  1  5  8  8  4  2  1  0  0  0
     0  0  0  4 17 32 37 27 14  6  2  0  0     0  0  1  7 31 58 64 43 21  8  3  1  0     0  0  1  7 24 38 40 26 12  5  2  0  0
     0  0  1  5 24 48 59 46 24 11  4  1  0     0  0  1 11 46 90100 72 37 16  6  2  0     0  0  1 11 41 67 65 44 22  9  3  1  0
     0  0  0  3 11 27 36 28 13  6  2  1  0     0  0  1  5 24 51 64 46 22 10  4  1  0     0  0  1  5 22 40 43 30 14  6  2  1  0
     0  0  0  1  3  6  6  6  4  2  0  0  0     0  0  0  1  5 10 11 10  6  3  1  0  0     0  0  0  1  5  8  8  7  3  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 10  9  6  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  4 12 18 16 11  5  2  1  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 11 12  8  4  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    623 REFLECTION INTENSITIES INTEGRATED


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
   301   0  0.971  5210867    0    128      92     0  0.04998  0.15187
   302   0  0.970  5211248    0    135      92     0  0.04413  0.16341
   303   0  0.967  5210337    0    132      96     0  0.04809  0.15904
   304   0  0.966  5210091    0    133      98     0  0.04639  0.14339
   305   0  0.967  5210673    0    132      97     0  0.05349  0.15735
   306   0  0.967  5209486    0    140     106     0  0.04573  0.15899
   307   0  0.967  5209411    0    143     102     2  0.04736  0.13903
   308   0  0.964  5209024    0    142     101     2  0.04677  0.12858
   309   0  0.965  5209622    0    152     107     2  0.04748  0.13028
   310   0  0.969  5208966    0    152     113     5  0.04582  0.11343
   311   0  0.966  5208528    0    157     119     3  0.04615  0.12226
   312   0  0.965  5208332    0    159     123     3  0.04293  0.12989
   313   0  0.963  5209595    0    160     106     3  0.04327  0.13116
   314   0  0.964  5210532    0    162      98     2  0.04761  0.12676
   315   0  0.965  5210476    0    147      95     2  0.04366  0.10979
   316   0  0.966  5210898    0    145      97     2  0.04541  0.12014
   317   0  0.961  5211077    0    129      91     4  0.04561  0.09908
   318   0  0.961  5211910    0    133      82     2  0.04394  0.12641
   319   0  0.964  5211581    0    128      86     2  0.04420  0.11837
   320   0  0.962  5211361    0    133      88     4  0.04520  0.09106

   515 OUT OF    518 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.52
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.11
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.424    32.138    34.507  88.657 108.483 111.899
 REC. CELL PARAMETERS   0.041668  0.033732  0.030735  83.816  70.589  67.324
 COORDINATES OF UNIT CELL A-AXIS   -14.353   -14.226    18.539
 COORDINATES OF UNIT CELL B-AXIS    14.491    25.530    13.079
 COORDINATES OF UNIT CELL C-AXIS   -21.481    21.549   -16.277
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.000    -0.001     0.005
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001717 -0.000122
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000660  0.005176  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.29   1244.05
 DETECTOR ORIGIN (PIXELS) AT                     1265.42   1237.25
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.21
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     515 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   58     33     68     70     60     32     55     67     58
 SIGMAB (degree)       0.034  0.035  0.035  0.034  0.033  0.034  0.034  0.034  0.034
 SIGMAR (degree)       0.032  0.027  0.028  0.030  0.029  0.028  0.028  0.030  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    515 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  1  4  8  8  6  3  1  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  1  5 11 13 11  6  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  1  3  6  8  6  3  1  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  3  5  5  4  2  1  0  0  0     0  0  0  2  6 10 10  7  3  1  0  0  0     0  0  0  2  5  8  7  4  2  1  0  0  0
     0  0  0  4 14 29 34 24 12  5  2  0  0     0  0  1  7 28 54 59 40 19  8  3  1  0     0  0  1  7 25 39 38 24 11  4  2  0  0
     0  0  0  5 22 46 54 43 21  9  3  1  0     0  0  1 11 46 89100 72 35 15  6  2  1     0  0  1 11 40 67 65 45 21  9  3  1  0
     0  0  0  3 11 25 32 25 13  6  2  1  0     0  0  1  6 23 49 61 44 23  9  4  1  0     0  0  1  6 19 38 41 27 14  6  2  1  0
     0  0  0  0  2  4  5  5  3  1  0  0  0     0  0  0  1  5 10 11  9  5  2  1  0  0     0  0  0  1  4  7  8  7  3  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  1  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 10  9  5  3  1  0  0  0     0  0  0  1  1  1  1  1  1  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  3 11 17 15 10  5  2  1  0  0     0  0  0  1  1  2  2  1  1  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  2  5 10 10  7  3  2  0  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    672 REFLECTION INTENSITIES INTEGRATED


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
   321   0  0.962  5211378    0    134      85     6  0.04357  0.07669
   322   0  0.962  5212849    0    130      72     3  0.04189  0.08939
   323   0  0.965  5212532    0    125      72     6  0.04440  0.07754
   324   0  0.966  5211869    0    126      87     8  0.04608  0.10136
   325   0  0.964  5210917    0    138      97     8  0.04376  0.09658
   326   0  0.964  5210139    0    143     106     4  0.04583  0.11657
   327   0  0.967  5210648    0    146      98     3  0.04449  0.10192
   328   0  0.966  5210687    0    138      95     2  0.04550  0.12734
   329   0  0.966  5210648    0    132      95     0  0.04591  0.12650
   330   0  0.968  5211260    0    130      90     0  0.04528  0.12342
   331   0  0.966  5212074    0    129      78     0  0.04502  0.12177
   332   0  0.964  5212436    0    113      75     0  0.04897  0.13334
   333   0  0.967  5212314    0    101      80     1  0.04723  0.14335
   334   0  0.964  5212231    0    113      79     0  0.04745  0.16239
   335   0  0.958  5211250    0    123      88     1  0.04561  0.13301
   336   0  0.959  5210676    0    133      94     3  0.04240  0.09095
   337   0  0.962  5210612    0    141      95     6  0.04327  0.08000
   338   0  0.962  5211671    0    137      86     3  0.04285  0.10576
   339   0  0.965  5211789    0    136      86     4  0.04821  0.10784
   340   0  0.965  5211804    0    132      83     4  0.04812  0.11058

   454 OUT OF    455 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.50
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.15
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.429    32.136    34.510  88.629 108.472 111.901
 REC. CELL PARAMETERS   0.041656  0.033732  0.030728  83.851  70.612  67.332
 COORDINATES OF UNIT CELL A-AXIS   -14.357   -14.228    18.541
 COORDINATES OF UNIT CELL B-AXIS    14.488    25.531    13.077
 COORDINATES OF UNIT CELL C-AXIS   -21.483    21.561   -16.265
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.006    -0.003     0.008
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001717 -0.000122
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000855  0.005005  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.33   1244.04
 DETECTOR ORIGIN (PIXELS) AT                     1265.21   1237.46
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.24
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     454 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     30     51     67     49     27     56     65     54
 SIGMAB (degree)       0.034  0.035  0.035  0.034  0.033  0.035  0.034  0.034  0.034
 SIGMAR (degree)       0.031  0.028  0.029  0.031  0.031  0.026  0.028  0.030  0.031


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    454 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  1  4  9 11  9  5  2  1  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  2  2  1  1  0  0  0     0  0  0  1  6 14 18 16  8  4  1  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  3  7 10  9  5  2  1  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  1  0  0  0  0  0  0  0  0  0  0  0     0  1  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  4  7  7  4  2  1  0  0  0     0  0  0  2  8 13 11  6  3  1  0  0  0     0  0  0  1  5  9  8  4  2  1  0  0  0
     0  0  1  5 16 34 39 28 14  6  2  0  0     0  0  1 10 34 62 61 39 19  7  3  1  0     0  0  1  7 23 40 39 24 11  4  2  0  0
     0  0  1  5 24 53 63 51 26 12  4  1  0     0  0  1 13 52 96100 70 34 15  5  1  0     0  0  1 10 37 61 62 42 21  9  3  1  0
     0  0  0  2 13 29 38 32 15  7  3  1  0     0  0  0  6 26 52 61 44 21  9  3  1  0     0  0  1  4 18 33 38 27 12  5  2  0  0
     0  0  0  0  2  5  7  6  3  1  0  0  0     0  0  0  1  5  9 12  9  5  2  1  0  0     0  0  0  1  4  6  7  5  3  1  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0 -1  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5 11 10  6  3  1  0  0  0     0  0  0  1  1  1  1  1  1  0  0  0  1     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  9 15 15 10  5  2  1  0  0     0  0  0  1  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  8  9  7  3  1  0  0  0     0  0  0  1  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    617 REFLECTION INTENSITIES INTEGRATED


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
   341   0  0.966  5210931    0    138      91     3  0.04593  0.11082
   342   0  0.966  5209606    0    139     104     6  0.04355  0.08518
   343   0  0.966  5210491    0    148      98     6  0.04423  0.09418
   344   0  0.969  5210825    0    147      91     0  0.04537  0.13913
   345   0  0.971  5211207    0    132      88     0  0.04616  0.15283
   346   0  0.968  5211914    0    136      84     1  0.04637  0.14243
   347   0  0.967  5212013    0    135      78     2  0.04445  0.10623
   348   0  0.968  5212245    0    122      77     1  0.04424  0.11236
   349   0  0.967  5210893    0    137      93     1  0.04357  0.10994
   350   0  0.966  5209988    0    144      96     5  0.04244  0.08622
   351   0  0.965  5209859    0    149      96     5  0.04388  0.09516
   352   0  0.965  5210325    0    151      94     3  0.04429  0.12069
   353   0  0.969  5210684    0    143      92     1  0.04416  0.13608
   354   0  0.967  5210897    0    132      88     8  0.04368  0.09081
   355   0  0.968  5211111    0    135      93     2  0.04708  0.13567
   356   0  0.965  5211312    0    134      88     1  0.04464  0.15637
   357   0  0.963  5210970    0    134      90     8  0.04239  0.08668
   358   0  0.963  5211863    0    137      85     7  0.04633  0.09614
   359   0  0.964  5211610    0    130      77     4  0.04464  0.08625
   360   0  0.962  5211656    0    134      80     5  0.04393  0.08354

   497 OUT OF    501 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.60
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.14
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.428    32.133    34.516  88.653 108.487 111.914
 REC. CELL PARAMETERS   0.041669  0.033741  0.030727  83.813  70.584  67.309
 COORDINATES OF UNIT CELL A-AXIS   -14.354   -14.230    18.540
 COORDINATES OF UNIT CELL B-AXIS    14.490    25.528    13.071
 COORDINATES OF UNIT CELL C-AXIS   -21.486    21.554   -16.283
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.003    -0.001     0.006
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001717 -0.000122
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000648  0.005247  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.31   1244.07
 DETECTOR ORIGIN (PIXELS) AT                     1265.46   1237.17
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.26
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     497 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   48     33     66     71     54     29     54     69     60
 SIGMAB (degree)       0.034  0.035  0.035  0.034  0.033  0.035  0.035  0.034  0.034
 SIGMAR (degree)       0.031  0.028  0.029  0.030  0.029  0.027  0.029  0.031  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    497 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  1  2  2  3  2  1  1  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  2  7 11 14 10  5  2  1  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  2  3  3  1  1  0  0  0     0  0  0  2  9 16 21 17 10  4  2  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  1  4  9 12  9  6  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
    -1  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  1  0  0  0     0  0  0  0  0  0  0  0  1  1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  4  7  8  6  2  1  0  0  0     0  0  0  1  5 11 12  7  3  1  0  0  0     0  0  0  1  3  6  6  3  2  1  0  0  0
     0  0  1  4 19 36 42 32 16  6  2  1  0     0  0  1  7 29 53 59 42 19  8  3  1  0     0  0  1  5 19 29 30 20  9  4  1  0  0
     0  0  1  6 27 56 69 54 29 13  5  2  0     0  0  1 11 46 91100 72 36 15  6  2  0     0  0  1  9 31 53 52 35 17  7  3  1  0
     0  0  0  3 12 29 42 31 17  8  3  1  0     0  0  1  6 23 48 60 42 22 10  4  1  0     0  0  1  5 17 29 32 22 11  5  2  0  0
     0  0  0  1  2  5  7  8  4  2  1  0  0     0  0  0  1  4  9 11 10  5  2  1  0  0     0  0  0  1  3  6  6  5  3  1  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  6  6  3  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 11 10  7  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  6  7  4  2  1  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    635 REFLECTION INTENSITIES INTEGRATED


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
   361   0  0.963  5211861    0    132      84     5  0.04784  0.10658
   362   0  0.962  5211179    0    148      91     4  0.04785  0.10411
   363   0  0.965  5210403    0    142     100     6  0.04470  0.08916
   364   0  0.961  5210569    0    147     100     8  0.04282  0.08667
   365   0  0.963  5210522    0    152      95     7  0.04547  0.08861
   366   0  0.965  5211369    0    131      94     4  0.04427  0.10350
   367   0  0.966  5211887    0    126      88     4  0.04224  0.11205
   368   0  0.964  5212287    0    117      82     3  0.04608  0.13340
   369   0  0.963  5212326    0    120      77     5  0.04213  0.09096
   370   0  0.964  5212236    0    127      77     6  0.04134  0.08100
   371   0  0.965  5211472    0    135      86     6  0.04326  0.10178
   372   0  0.963  5209564    0    135     106     5  0.04530  0.12086
   373   0  0.963  5209461    0    142     111     3  0.04526  0.12903
   374   0  0.963  5209563    0    152     101     1  0.04720  0.13303
   375   0  0.963  5209525    0    157     103     5  0.04520  0.11014
   376   0  0.964  5209262    0    145     108     4  0.04523  0.10595
   377   0  0.961  5210464    0    138      96     5  0.04282  0.09080
   378   0  0.963  5211372    0    134      82     2  0.05188  0.11874
   379   0  0.961  5210687    0    124      90     3  0.04695  0.12311
   380   0  0.963  5210174    0    145      95     1  0.04768  0.14604

   481 OUT OF    484 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.57
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.11
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.427    32.136    34.510  88.616 108.461 111.902
 REC. CELL PARAMETERS   0.041655  0.033731  0.030725  83.871  70.629  67.337
 COORDINATES OF UNIT CELL A-AXIS   -14.354   -14.230    18.540
 COORDINATES OF UNIT CELL B-AXIS    14.489    25.529    13.078
 COORDINATES OF UNIT CELL C-AXIS   -21.487    21.568   -16.248
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.018    -0.014     0.008
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001717 -0.000122
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.001010  0.004903  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.39   1244.03
 DETECTOR ORIGIN (PIXELS) AT                     1265.07   1237.59
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.25
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     481 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   50     33     47     67     64     33     52     63     58
 SIGMAB (degree)       0.034  0.035  0.035  0.034  0.033  0.035  0.035  0.035  0.034
 SIGMAR (degree)       0.032  0.027  0.029  0.031  0.029  0.027  0.029  0.031  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    481 STRONG REFLECTIONS *****


     0  1  1  0  0  0  0  0  0  0  0  0  0     0  1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  1  1  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  1  4  9 12  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  1  6 14 18 15  8  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  1  3  7 10  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  1  2  3  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  6  7  4  2  1  0  0  0     0  0  0  2  7 11 10  6  3  1  0  0  0     0  0  0  2  5  7  6  3  2  0  0  0  0
     0  0  0  5 16 33 41 29 15  6  2  0  0     0  0  1  9 31 58 64 40 19  8  3  1  0     0  0  1  7 22 37 37 22 11  4  1  0  0
     0  0  1  6 22 52 62 48 27 11  4  1  0     0  0  1 11 46 95100 68 35 15  5  2  0     0  0  1  9 34 60 59 38 20  8  3  1  0
     0  0  0  3 13 28 35 30 16  6  3  1  0     0  0  1  6 27 54 58 45 22  8  3  1  0     0  0  1  5 19 34 36 26 13  5  2  0  0
     0  0  0  1  2  5  9  7  4  2  0  0  0     0  0  0  1  4 11 15 11  5  2  1  0  0     0  0  0  1  3  7  9  7  3  1  0  0  0
     0  0  0  0  0  1  1  1  1  0  0 -1  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  8  8  5  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  2  9 14 13  8  4  2  1  0  0     0  0  0  0  1  2  1  1  1  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  5  8  9  6  3  1  0  0  0     0  0  0  0  1  1  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    640 REFLECTION INTENSITIES INTEGRATED


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
   381   0  0.962  5209421    0    154     101    11  0.04372  0.08309
   382   0  0.963  5209894    0    155      97     6  0.04709  0.10303
   383   0  0.965  5210254    0    156      96     8  0.04671  0.09901
   384   0  0.967  5210403    0    146      95     7  0.04677  0.10164
   385   0  0.963  5210780    0    138      95     3  0.04591  0.14959
   386   0  0.960  5210665    0    146      93     2  0.04570  0.12705
   387   0  0.961  5209218    0    145     112     4  0.04456  0.11123
   388   0  0.964  5209861    0    154      99     7  0.04323  0.07819
   389   0  0.961  5210705    0    152      90     3  0.04625  0.09513
   390   0  0.963  5210072    0    141      93     3  0.04296  0.08959
   391   0  0.961  5210208    0    145      95     3  0.04150  0.08564
   392   0  0.959  5210509    0    148      89     6  0.04183  0.08710
   393   0  0.960  5210812    0    150      88     7  0.04531  0.09360
   394   0  0.960  5210763    0    142      88     5  0.04527  0.09672
   395   0  0.960  5211112    0    144      84     9  0.04554  0.07745
   396   0  0.960  5211397    0    142      78     7  0.04535  0.08011
   397   0  0.960  5211298    0    143      84     5  0.04466  0.09787
   398   0  0.957  5210974    0    145      84     7  0.04663  0.08156
   399   0  0.958  5210716    0    137      88     9  0.04418  0.08058
   400   0  0.956  5210454    0    136      90     6  0.04420  0.08830

   510 OUT OF    512 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.52
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.13
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.420    32.131    34.503  88.648 108.483 111.881
 REC. CELL PARAMETERS   0.041668  0.033734  0.030737  83.832  70.595  67.346
 COORDINATES OF UNIT CELL A-AXIS   -14.351   -14.219    18.540
 COORDINATES OF UNIT CELL B-AXIS    14.480    25.529    13.077
 COORDINATES OF UNIT CELL C-AXIS   -21.478    21.545   -16.278
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.009     0.004     0.004
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001717 -0.000122
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000780  0.005125  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.28   1244.10
 DETECTOR ORIGIN (PIXELS) AT                     1265.26   1237.37
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.16
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     510 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   51     34     65     67     59     33     57     72     59
 SIGMAB (degree)       0.034  0.035  0.035  0.034  0.033  0.035  0.035  0.035  0.034
 SIGMAR (degree)       0.030  0.029  0.030  0.031  0.029  0.027  0.029  0.031  0.032


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    510 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  1  3  6  8  7  4  2  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  4  9 13 11  6  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  2  5  7  6  4  1  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  3  5  6  4  2  1  0  0  0     0  0  0  2  6 10 11  7  3  1  0  0  0     0  0  0  2  5  7  8  5  2  1  0  0  0
     0  0  1  3 14 27 33 24 12  5  2  0  0     0  0  1  7 31 56 60 39 19  8  3  1  0     0  0  1  8 27 46 44 28 13  5  2  0  0
     0  0  1  5 19 41 53 43 23  9  4  1  0     0  0  1 12 46 89100 71 37 15  6  2  0     0  0  2 13 45 75 76 52 25 10  4  1  0
     0  0  0  2  9 23 30 26 14  6  2  1  0     0  0  1  6 22 48 59 44 23 10  4  1  0     0  0  1  7 22 41 46 33 17  7  3  1  0
     0  0  0  0  2  4  5  5  3  1  0  0  0     0  0  0  1  4  8 11  9  5  2  1  0  0     0  0  0  1  4  8  9  7  4  2  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  2  3  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  9 14 14  9  4  2  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  1  5 14 23 23 16  8  3  1  0  0     0  0  0  0  1  3  3  2  1  0  0  0  0     0  0  0  0  1  0  1  1  0  0  0  0  0
     0  0  0  2  7 13 15 11  6  2  1  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  3  3  3  1  1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    681 REFLECTION INTENSITIES INTEGRATED


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
   401   0  0.957  5209261    0    132     107     5  0.04343  0.09683
   402   0  0.956  5209667    0    139     101     3  0.04315  0.11070
   403   0  0.957  5210088    0    141      91     3  0.04220  0.10167
   404   0  0.955  5210154    0    131      96     9  0.04509  0.08300
   405   0  0.955  5210807    0    124      92     7  0.04556  0.10129
   406   0  0.955  5211808    0    117      77     5  0.04812  0.09167
   407   0  0.957  5212147    0    109      73     2  0.04543  0.12357
   408   0  0.956  5212385    0    107      74     2  0.04364  0.12735
   409   0  0.959  5212564    0    114      73     2  0.04681  0.10941
   410   0  0.956  5211709    0    115      79     2  0.04330  0.09619
   411   0  0.954  5211433    0    126      80     5  0.04301  0.09192
   412   0  0.957  5211200    0    141      83     4  0.04641  0.09740
   413   0  0.960  5211106    0    134      87     1  0.04574  0.13005
   414   0  0.960  5210736    0    156      91     7  0.04324  0.07907
   415   0  0.957  5210328    0    167      94     4  0.04333  0.07897
   416   0  0.955  5209184    0    164     100     7  0.04290  0.07131
   417   0  0.954  5210642    0    161      89     5  0.04315  0.07485
   418   0  0.957  5210747    0    147      88     4  0.04389  0.09479
   419   0  0.956  5210526    0    138      93     3  0.04524  0.10873
   420   0  0.955  5211032    0    129      86     6  0.04790  0.10044

   484 OUT OF    486 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.55
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.424    32.132    34.511  88.629 108.476 111.921
 REC. CELL PARAMETERS   0.041670  0.033742  0.030728  83.842  70.605  67.311
 COORDINATES OF UNIT CELL A-AXIS   -14.354   -14.230    18.536
 COORDINATES OF UNIT CELL B-AXIS    14.489    25.529    13.070
 COORDINATES OF UNIT CELL C-AXIS   -21.484    21.560   -16.266
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.006    -0.005     0.006
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001717 -0.000122
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000776  0.005096  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.38   1244.02
 DETECTOR ORIGIN (PIXELS) AT                     1265.37   1237.33
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.23
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     484 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   50     29     57     71     57     31     52     65     55
 SIGMAB (degree)       0.035  0.035  0.036  0.035  0.033  0.035  0.035  0.035  0.034
 SIGMAR (degree)       0.029  0.027  0.028  0.031  0.029  0.026  0.029  0.031  0.031


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    484 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  1  2  3  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  2  7 12 14 10  6  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  3  2  1  0  0  0     0  0  0  2  9 17 21 18 10  5  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  1  4 10 13 12  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  3  3  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  5  9  8  5  3  1  0  0  0     0  0  0  2  7 14 12  6  3  1  0  0  0     0  0  0  1  4  8  6  3  1  0  0  0  0
     0  0  1  5 19 36 44 31 17  7  2  1  0     0  0  1  8 30 57 64 40 20  8  3  1  0     0  0  1  6 19 32 32 18  9  3  1  0  0
     0  0  1  6 27 56 70 58 31 13  5  2  0     0  0  1 11 49 93100 73 37 15  5  2  0     0  0  1  8 31 52 51 35 17  7  2  1  0
     0  0  0  3 14 33 43 36 19  9  3  1  0     0  0  1  7 25 54 61 47 22 10  4  1  0     0  0  1  5 16 30 30 23 10  4  1  0  0
     0  0  0  1  3  5  9  8  4  2  1  0  0     0  0  0  2  6 11 14 10  5  3  1  0  0     0  0  0  1  4  6  8  4  3  1  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4  7  6  4  2  1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 11 10  7  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  6  6  4  2  1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    633 REFLECTION INTENSITIES INTEGRATED


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
   421   0  0.954  5211982    0    132      79     5  0.04989  0.10930
   422   0  0.957  5213090    0    123      68     4  0.04765  0.12240
   423   0  0.960  5212139    0    118      79     0  0.04640  0.15512
   424   0  0.960  5212099    0    113      79     1  0.04630  0.15479
   425   0  0.960  5211449    0    123      84     9  0.04278  0.08404
   426   0  0.959  5210995    0    120      95     8  0.04432  0.09111
   427   0  0.959  5211963    0    137      85     8  0.04685  0.10467
   428   0  0.960  5211981    0    128      78     7  0.04722  0.10069
   429   0  0.959  5211785    0    140      84     8  0.04866  0.10792
   430   0  0.956  5210468    0    131      94     8  0.04840  0.09531
   431   0  0.958  5210551    0    121     101     5  0.04372  0.09846
   432   0  0.957  5212418    0    133      75     3  0.04435  0.10725
   433   0  0.958  5212705    0    131      69     3  0.04576  0.09082
   434   0  0.958  5211114    0    131      90     2  0.04142  0.09337
   435   0  0.957  5211049    0    132      89     1  0.04290  0.11125
   436   0  0.959  5211483    0    143      83     5  0.04484  0.08314
   437   0  0.956  5210995    0    147      88     4  0.04275  0.07809
   438   0  0.954  5210182    0    143      95     5  0.03889  0.07519
   439   0  0.956  5211473    0    134      81     5  0.04255  0.08434
   440   0  0.958  5212079    0    124      78     0  0.04688  0.14713

   479 OUT OF    482 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.60
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.11
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.419    32.118    34.499  88.607 108.465 111.907
 REC. CELL PARAMETERS   0.041669  0.033751  0.030735  83.877  70.628  67.334
 COORDINATES OF UNIT CELL A-AXIS   -14.349   -14.221    18.538
 COORDINATES OF UNIT CELL B-AXIS    14.480    25.521    13.061
 COORDINATES OF UNIT CELL C-AXIS   -21.479    21.559   -16.249
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.008    -0.010    -0.000
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001717 -0.000122
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.001124  0.005033  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.42   1244.04
 DETECTOR ORIGIN (PIXELS) AT                     1264.95   1237.43
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.13
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     479 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   48     33     54     64     51     27     53     70     62
 SIGMAB (degree)       0.034  0.035  0.036  0.034  0.033  0.035  0.035  0.035  0.035
 SIGMAR (degree)       0.032  0.029  0.030  0.031  0.030  0.027  0.028  0.031  0.032


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    479 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  1  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  2  5 10 12  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  1  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0     0  0  0  2  7 15 19 15  9  4  2  0  0
     0  0  0  0  1  0  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  3  7 11 10  6  2  1  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4  8  9  5  3  1  0  0  0     0  0  0  2  6 12 12  7  3  1  0  0  0     0  0  0  1  4  9  8  4  2  1  0  0  0
     0  0  1  5 17 35 39 28 15  6  2  0  0     0  0  1  9 32 58 59 38 19  8  3  1  0     0  0  1  6 22 40 40 24 11  5  2  1  0
     0  0  1  6 25 55 65 50 28 12  5  1  0     0  0  1 12 48 92100 70 37 15  6  2  0     0  0  1 10 36 61 66 44 22  9  3  1  0
     0  0  0  3 13 27 37 30 17  7  3  1  0     0  0  1  6 26 49 59 44 23 10  4  1  0     0  0  1  5 19 35 39 28 15  6  2  1  0
     0  0  0  0  2  3  6  6  3  2  0  0  0     0  0  0  1  4  8 11  9  4  2  1  0  0     0  0  0  1  4  6  8  7  3  2  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  3  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  2  6 12 11  7  3  1  0  0  0     0  0  0  1  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  3 10 16 18 12  6  2  1  0  0     0  0  0  1  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 10 12  8  4  2  1  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  3  3  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    617 REFLECTION INTENSITIES INTEGRATED


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
   441   0  0.961  5212847    0    125      70     2  0.04714  0.12961
   442   0  0.959  5211909    0    131      79     5  0.04776  0.08944
   443   0  0.954  5210890    0    139      90     9  0.04403  0.07624
   444   0  0.962  5211058    0    141      86     5  0.04323  0.08673
   445   0  0.959  5211097    0    143      85     4  0.04477  0.08463
   446   0  0.958  5211413    0    135      81     1  0.04421  0.10726
   447   0  0.960  5211722    0    129      79     3  0.04568  0.10167
   448   0  0.959  5210946    0    143      93     7  0.04284  0.11367
   449   0  0.958  5209358    0    140     108     4  0.04614  0.12579
   450   0  0.958  5209797    0    147     103     5  0.04472  0.10992
   451   0  0.958  5209974    0    149      97     2  0.04808  0.12768
   452   0  0.956  5209817    0    139     104     3  0.04365  0.11924
   453   0  0.957  5210467    0    141      95     4  0.04698  0.08873
   454   0  0.960  5211299    0    141      89     3  0.04502  0.09989
   455   0  0.959  5211555    0    143      84     2  0.04491  0.11703
   456   0  0.960  5211609    0    137      82     4  0.04269  0.07135
   457   0  0.958  5210894    0    150      83     3  0.04121  0.07152
   458   0  0.957  5210105    0    152      93     2  0.04330  0.08773
   459   0  0.960  5210066    0    153      95     3  0.04290  0.08469
   460   0  0.960  5211231    0    147      85     1  0.04389  0.12657

   508 OUT OF    509 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.55
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.427    32.141    34.511  88.636 108.475 111.886
 REC. CELL PARAMETERS   0.041656  0.033723  0.030728  83.848  70.607  67.345
 COORDINATES OF UNIT CELL A-AXIS   -14.354   -14.225    18.542
 COORDINATES OF UNIT CELL B-AXIS    14.486    25.535    13.083
 COORDINATES OF UNIT CELL C-AXIS   -21.484    21.558   -16.270
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.002    -0.002     0.006
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001717 -0.000122
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000857  0.005040  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.31   1244.02
 DETECTOR ORIGIN (PIXELS) AT                     1265.19   1237.40
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.25
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     508 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   49     35     62     68     55     30     60     74     62
 SIGMAB (degree)       0.035  0.035  0.035  0.034  0.033  0.035  0.035  0.035  0.034
 SIGMAR (degree)       0.031  0.031  0.030  0.030  0.029  0.026  0.028  0.031  0.031


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    508 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  1  1  2  2  1  1  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  1  5  8  9  7  4  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  2  6 12 14 12  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  1  3  6  8  7  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  1  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  2  4  6  6  4  2  1  0  0  0     0  0  0  2  7 10 11  7  3  1  0  0  0     0  0  0  2  5  7  7  4  2  1  0  0  0
     0  0  1  5 17 31 33 24 12  5  2  0  0     0  0  1  8 31 57 58 40 18  8  2  1  0     0  0  1  7 23 38 36 25 11  5  1  0  0
     0  0  1  7 24 48 57 44 24 10  4  1  0     0  0  1 12 47 91100 72 38 16  6  2  0     0  0  1 11 36 63 63 44 22  9  3  1  0
     0  0  0  3 11 26 33 26 15  6  3  1  0     0  0  1  6 24 51 60 43 24  9  4  1  0     0  0  1  6 20 35 38 27 14  5  2  0  0
     0  0  0  0  2  4  6  6  3  1  0  0  0     0  0  0  1  4  9 11  9  5  2  1  0  0     0  0  0  1  4  7  7  5  3  1  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 10  9  6  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3 10 16 15 10  5  2  1  0  0     0  0  0  1  2  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  9  9  7  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  2  1  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
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
   461   0  0.963  5212018    0    132      83     1  0.05034  0.14434
   462   0  0.960  5212309    0    126      75     4  0.05047  0.11866
   463   0  0.959  5211478    0    113      81     8  0.04447  0.08800
   464   0  0.960  5212657    0    122      62     7  0.04421  0.05639
   465   0  0.961  5212311    0    128      69     7  0.04694  0.07859
   466   0  0.961  5211700    0    127      76     0  0.04614  0.15010
   467   0  0.963  5211696    0    115      83     0  0.04629  0.16906
   468   0  0.962  5212609    0    108      68     0  0.04802  0.16801
   469   0  0.958  5212984    0    108      69     3  0.04598  0.11582
   470   0  0.961  5212961    0    122      67     1  0.04874  0.12107
   471   0  0.964  5212515    0    123      74     2  0.04601  0.11195
   472   0  0.961  5211147    0    128      87     1  0.04606  0.12652
   473   0  0.960  5211211    0    124      83     1  0.04791  0.13169
   474   0  0.960  5210681    0    122      89     3  0.04315  0.10979
   475   0  0.960  5210588    0    121      90     4  0.04698  0.10739
   476   0  0.960  5211773    0    132      80     1  0.04912  0.12986
   477   0  0.959  5211397    0    147      83     3  0.04894  0.10987
   478   0  0.957  5208969    0    155     113    12  0.04560  0.08369
   479   0  0.957  5208029    0    157     117     7  0.04701  0.10187
   480   0  0.959  5209611    0    153     101     7  0.04677  0.09061

   468 OUT OF    469 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.64
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.15
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.417    32.150    34.515  88.668 108.488 111.894
 REC. CELL PARAMETERS   0.041681  0.033719  0.030729  83.802  70.579  67.324
 COORDINATES OF UNIT CELL A-AXIS   -14.353   -14.226    18.528
 COORDINATES OF UNIT CELL B-AXIS    14.487    25.541    13.092
 COORDINATES OF UNIT CELL C-AXIS   -21.483    21.548   -16.293
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.012     0.012     0.013
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001717 -0.000122
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000323  0.005091  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.20   1243.98
 DETECTOR ORIGIN (PIXELS) AT                     1265.78   1237.30
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.22
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     468 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   45     29     58     67     56     27     53     64     50
 SIGMAB (degree)       0.035  0.035  0.035  0.034  0.033  0.035  0.035  0.035  0.034
 SIGMAR (degree)       0.032  0.029  0.028  0.029  0.029  0.027  0.030  0.031  0.031


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    468 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  1  1  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  1  5  8  9  7  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  2  5 11 14 12  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  3  5  8  7  4  2  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  0  0  0  0  0 -1  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  2  4  6  7  5  2  1  0  0  0     0  0  0  2  6 10 11  7  3  1  0  0  0     0  0  0  2  5  8  8  4  2  1  0  0  0
     0  0  1  4 17 32 37 27 15  6  2  1  0     0  0  1  7 30 56 62 43 22  9  3  1  0     0  0  1  6 25 40 40 27 13  5  2  0  0
     0  0  1  6 22 46 56 46 26 10  4  1  0     0  0  1 11 44 88100 76 40 16  6  2  1     0  0  1 10 38 67 69 48 24  9  4  1  0
     0  0  0  3 12 23 33 26 15  7  3  1  0     0  0  1  6 23 44 60 43 24 10  4  1  1     0  0  1  5 19 34 42 28 14  6  2  1  0
     0  0  0  0  2  5  6  6  4  1  1  0  0     0  0  0  1  4  9 12 10  6  2  1  0  0     0  0  0  1  3  7  7  6  3  1  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  1  1  2  0  0  0 -1  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  2  8 11 10  6  3  1  0  0  0     0  0  0  0  1  2  1  1  1  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0
     0  0  1  4 12 19 18 12  5  2  1  0  0     0  0  0  1  2  2  2  1  1  0  0  0  0     0  0  0  1  1  1  1  1  0  0  0  0  0
     0  0  0  2  6  9 11  7  3  1  0  0  0     0  0  0  0  1  1  2  2  0  0  0  0  0     0  0  0  1  0  0  0  1  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    624 REFLECTION INTENSITIES INTEGRATED


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
   481   0  0.957  5210436    0    143      94     3  0.04735  0.10422
   482   0  0.959  5210496    0    147      94     2  0.04915  0.13047
   483   0  0.962  5210448    0    146     103     2  0.04958  0.13887
   484   0  0.959  5209292    0    151     108     2  0.04759  0.13120
   485   0  0.959  5209678    0    151     102     2  0.04691  0.12317
   486   0  0.957  5210184    0    143      93     3  0.04778  0.10590
   487   0  0.958  5210677    0    141      82     2  0.04422  0.09637
   488   0  0.959  5210849    0    124      86     1  0.04391  0.10497
   489   0  0.960  5211868    0    131      77     5  0.04583  0.09688
   490   0  0.961  5211859    0    134      76     5  0.04514  0.08965
   491   0  0.959  5211606    0    131      81     6  0.04127  0.09321
   492   0  0.963  5211433    0    132      86     3  0.04678  0.10734
   493   0  0.966  5211912    0    128      84     3  0.04862  0.12968
   494   0  0.965  5211914    0    128      80     1  0.04392  0.13248
   495   0  0.965  5211413    0    131      91     1  0.04443  0.12224
   496   0  0.966  5212158    0    136      85     1  0.04282  0.12507
   497   0  0.966  5210851    0    139      96     6  0.04407  0.09166
   498   0  0.961  5211321    0    141      88     3  0.04272  0.11878
   499   0  0.965  5211340    0    140      91     5  0.04637  0.10891
   500   0  0.963  5210474    0    153      95     6  0.04507  0.09648

   486 OUT OF    488 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.60
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.13
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.425    32.135    34.516  88.627 108.475 111.918
 REC. CELL PARAMETERS   0.041667  0.033738  0.030723  83.846  70.608  67.315
 COORDINATES OF UNIT CELL A-AXIS   -14.355   -14.223    18.542
 COORDINATES OF UNIT CELL B-AXIS    14.486    25.538    13.062
 COORDINATES OF UNIT CELL C-AXIS   -21.487    21.556   -16.280
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.013     0.006    -0.005
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999999 -0.001717 -0.000122
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000686  0.005158  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.36   1244.09
 DETECTOR ORIGIN (PIXELS) AT                     1265.46   1237.32
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.25
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     486 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     32     57     67     56     31     56     67     54
 SIGMAB (degree)       0.035  0.035  0.035  0.034  0.034  0.035  0.035  0.035  0.035
 SIGMAR (degree)       0.030  0.026  0.028  0.031  0.030  0.026  0.028  0.030  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    486 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  2  3  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  1  1  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  2  5 12 15 12  6  3  1  0  0
     0  0  0  0  0  0  0  1  1  0  0  0  0     0  0  0  0  1  2  3  3  2  1  0  0  0     0  0  0  2  8 17 22 19 11  5  2  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  4  9 12 10  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6  8  9  6  3  1  0  0  0     0  0  0  2  8 11 11  7  3  1  0  0  0     0  0  0  1  5  7  6  3  2  1  0  0  0
     0  0  1  6 19 41 46 33 16  7  2  1  0     0  0  1 10 32 60 60 39 19  8  3  1  0     0  0  1  7 19 32 31 19  9  3  1  0  0
     0  0  1  8 29 61 72 55 30 13  5  1  0     0  0  2 16 52 94100 69 36 15  5  2  0     0  0  1 10 29 48 50 34 17  7  2  1  0
     0  0  0  3 15 33 41 34 19  8  3  1  0     0  0  1  7 27 52 60 45 24 10  4  1  0     0  0  1  4 15 27 29 22 11  5  2  0  0
     0  0  0  1  3  7  8  7  4  2  1  0  0     0  0  0  1  6 11 12 10  5  2  1  0  0     0  0  0  1  3  6  7  5  2  1  0  0  0
     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  6  5  3  1  1  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  9  9  6  3  1  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  5  6  4  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    658 REFLECTION INTENSITIES INTEGRATED


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
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  5 10 13 10  5  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  2  7 15 20 15  8  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  4  9 12  9  5  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  4  8  9  6  3  0  0  0  0     0  0  0  2  7 14 14  8  3  0  0  0  0     0  0  0  0  5  9  8  4  2  0  0  0  0
     0  0  0  5 19 38 44 30 14  5  0  0  0     0  0  0  9 33 61 64 41 19  7  0  0  0     0  0  0  7 23 38 36 22 10  3  0  0  0
     0  0  0  7 27 56 67 50 25  9  0  0  0     0  0  2 13 50 94100 70 33 12  0  0  0     0  0  0 10 35 59 57 38 17  6  0  0  0
     0  0  0  3 14 33 41 31 15  6  0  0  0     0  0  0  6 28 56 65 45 21  8  0  0  0     0  0  0  5 20 36 38 25 11  4  0  0  0
     0  0  0  0  3  7 10  7  3  0  0  0  0     0  0  0  2  6 13 16 11  5  2  0  0  0     0  0  0  0  5  9 10  6  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  6  9  8  5  3  0  0  0  0     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  9 15 13  9  4  2  0  0  0     0  0  0  0  2  3  3  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6  9  9  6  3  0  0  0  0     0  0  0  0  0  4  4  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
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
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  3  5  7  7  5  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  2  6  9 10  7  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  5  6  4  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  6  9  8  5  3  0  0  0  0     0  0  0  3  9 16 15 10  4  0  0  0  0     0  0  0  2  5 11 11  7  3  0  0  0  0
     0  0  0  5 18 30 32 20 10  4  0  0  0     0  0  0 10 35 65 69 44 21  9  0  0  0     0  0  0  7 25 45 49 31 15  6  0  0  0
     0  0  0  6 22 41 45 32 16  6  0  0  0     0  0  0 13 50 93100 69 34 14  6  0  0     0  0  0 10 38 68 72 49 24 10  5  0  0
     0  0  0  3 12 22 26 19 10  4  0  0  0     0  0  0  7 28 51 58 41 20  8  0  0  0     0  0  0  6 21 37 42 30 14  6  0  0  0
     0  0  0  0  3  5  6  4  3  0  0  0  0     0  0  0  3  7 12 14  9  4  2  0  0  0     0  0  0  3  5  9 10  7  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  4  4  3  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  5 10 11  8  5  3  0  0  0     0  0  0  0  2  3  4  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  8 15 17 12  6  4  0  0  0     0  0  0  0  3  3  4  4  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  5  8 11  8  4  3  0  0  0     0  0  0  0  0  3  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  4  4  3  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  3 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  6  8  7  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  0  0  0  0  0     0  0  0  2  5 11 15 13  7  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  2  0  0  0  0  0     0  0  0  0  4  8 11  9  5  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  2  3  3  3  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  4  6  4  2  0  0  0  0     0  0  0  0  4  8  9  6  3  0  0  0  0     0  0  0  0  3  6  6  4  3  0  0  0  0
     0  0  0  3 11 24 31 22 12  5  0  0  0     0  0  0  5 23 46 51 35 18  8  0  0  0     0  0  0  5 18 31 31 20 10  4  0  0  0
     0  0  0  5 23 48 59 47 25 11  5  0  0     0  0  2 12 48 91100 72 37 16  6  0  0     0  0  0 10 36 61 62 42 21  9  5  0  0
     0  0  0  5 18 36 44 34 19  8  5  0  0     0  0  0 10 35 66 75 55 28 12  5  0  0     0  0  0  8 26 45 47 32 16  7  0  0  0
     0  0  0  3  5 10 12  9  5  3  0  0  0     0  0  0  3  9 18 20 14  7  3  0  0  0     0  0  0  3  7 12 12  9  4  3  0  0  0
     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  0  0  0  2  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  2  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  7  7  5  3  0  0  0  0     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4 10 15 14  9  5  3  0  0  0     0  0  0  0  2  3  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4  8 12 11  7  4  3  0  0  0     0  0  0  0  0  3  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  4  4  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  4 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  4  3  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  0  0  0  0  0     0  0  0  0  4 11 17 15  9  5  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  4  4  4  0  0  0  0     0  0  0  2  6 15 24 22 13  6  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  3  0  0  0  0  0     0  0  0  0  3  8 12 10  6  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  8 10  8  4  0  0  0  0     0  0  0  0  6 12 13  8  4  2  0  0  0     0  0  0  0  4  7  7  4  2  0  0  0  0
     0  0  0  3 15 37 51 41 22  9  0  0  0     0  0  0  6 29 59 66 45 23  9  0  0  0     0  0  0  6 21 36 35 22 11  4  0  0  0
     0  0  0  5 24 57 75 61 33 14  6  0  0     0  0  0 11 46 91100 70 36 15  6  0  0     0  0  0  9 33 56 54 35 17  7  0  0  0
     0  0  0  3 13 29 38 30 16  7  0  0  0     0  0  0  6 24 46 53 36 18  7  0  0  0     0  0  0  5 17 29 29 18  9  4  0  0  0
     0  0  0  0  3  5  7  5  3  0  0  0  0     0  0  0  0  4  9 10  7  3  0  0  0  0     0  0  0  0  3  6  6  4  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  6 10  8  5  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4 10 14 12  7  4  0  0  0  0     0  0  0  0  0  2  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  5  7  7  4  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  5 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  4  0  0  0  0  0     0  0  0  0  3  7 12 12  7  5  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  5  5  0  0  0  0     0  0  0  0  4 11 20 21 13  6  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  4  0  0  0  0  0     0  0  0  0  2  5  9  9  6  5  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  4  7  5  3  0  0  0  0     0  0  0  0  4  7  9  6  3  0  0  0  0     0  0  0  0  3  5  5  3  0  0  0  0  0
     0  0  0  2 10 27 41 35 18  8  0  0  0     0  0  0  5 24 50 57 40 19  8  0  0  0     0  0  0  6 20 32 30 18  8  3  0  0  0
     0  0  0  3 18 48 70 59 32 14  6  0  0     0  0  0 10 45 91100 68 34 14  5  0  0     0  0  0 10 37 58 52 31 15  6  0  0  0
     0  0  0  2  8 23 33 26 14  6  0  0  0     0  0  0  5 22 45 49 31 15  6  0  0  0     0  0  0  5 18 29 26 15  7  3  0  0  0
     0  0  0  0  0  3  4  3  2  0  0  0  0     0  0  0  0  3  6  6  4  2  0  0  0  0     0  0  0  0  4  4  4  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  6  8  6  4  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4  9 13 10  5  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  4  6  5  3  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
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
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  4  7  7  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  6 11 11  7  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  4  7  8  5  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  4  6  4  2  0  0  0  0     0  0  0  3  9 13 12  6  3  0  0  0  0     0  0  0  4  7  9  6  4  2  0  0  0  0
     0  0  0  3  9 23 31 24 12  5  0  0  0     0  0  3  9 36 62 60 37 17  7  0  0  0     0  0  3 11 32 41 31 17  7  3  0  0  0
     0  0  0  3 15 38 52 43 22  9  0  0  0     0  0  3 13 55 99100 65 31 12  5  0  0     0  0  3 16 49 65 52 31 14  5  3  0  0
     0  0  0  3  7 21 32 26 14  6  0  0  0     0  0  0  7 28 55 59 39 19  8  0  0  0     0  0  0  8 24 36 31 18  8  3  0  0  0
     0  0  0  0  0  4  6  6  3  0  0  0  0     0  0  0  0  5 11 11  8  4  2  0  0  0     0  0  0  0  5  7  6  4  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  2  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4  6  7  5  3  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  5  9 11  8  5  3  2  0  0  0     0  0  0  0  0  3  4  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  5  6  5  4  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
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
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  4  7  8  6  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  2  3  0  0  0  0  0     0  0  0  2  6 13 15 12  7  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  2  0  0  0  0  0     0  0  0  0  3  7 10 10  6  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  4  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  4  5  3  0  0  0  0  0     0  0  0  2  6 10  9  5  2  0  0  0  0     0  0  0  3  4  6  5  3  2  0  0  0  0
     0  0  0  3 14 30 32 21 10  4  0  0  0     0  0  0  7 32 56 53 32 15  6  0  0  0     0  0  0  6 22 34 29 16  7  3  0  0  0
     0  0  0  5 25 53 61 47 25 11  6  0  0     0  0  0 11 50 95100 71 35 15  6  0  0     0  0  0  9 34 58 55 37 18  7  4  0  0
     0  0  0  3 12 28 40 36 21  9  5  0  0     0  0  0  5 23 50 65 52 28 12  5  0  0     0  0  0  4 16 32 37 27 14  6  3  0  0
     0  0  0  0  2  5  9 10  6  3  0  0  0     0  0  0  0  3  8 14 13  8  4  0  0  0     0  0  0  0  3  5  8  6  4  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  2  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  5  6  5  3  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4  8 12 11  7  4  3  0  0  0     0  0  0  0  0  3  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  4  8  8  6  3  3  0  0  0     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  2  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  8 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  2  3  2  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  2  0  0  0  0  0     0  0  0  3  7 13 14  9  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  3  4  3  2  0  0  0  0     0  0  0  3 11 21 23 17  9  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  2  0  0  0  0  0     0  0  0  2  5 11 14 11  6  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4  7  7  4  2  0  0  0  0     0  0  0  2  5  9  9  5  2  0  0  0  0     0  0  0  0  3  5  5  3  0  0  0  0  0
     0  0  0  6 21 40 42 27 13  5  0  0  0     0  0  0  7 28 54 56 35 16  6  0  0  0     0  0  0  4 16 30 31 19  9  3  0  0  0
     0  0  0  8 34 67 73 53 27 12  5  0  0     0  0  2 12 47 92100 72 36 15  6  0  0     0  0  0  7 27 53 57 40 20  8  0  0  0
     0  0  0  4 16 37 46 36 19  9  4  0  0     0  0  0  6 24 53 67 51 27 12  5  0  0     0  0  0  3 15 32 40 30 16  7  0  0  0
     0  0  0  0  2  5  9  8  5  3  0  0  0     0  0  0  0  3  8 13 12  7  3  0  0  0     0  0  0  0  3  6  8  8  4  2  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  4  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  4  7  7  5  3  0  0  0  0     0  0  0  0  0  4  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 13 14 10  5  3  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  4  8 11  8  4  3  0  0  0     0  0  0  0  0  3  4  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  2  3  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  9 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  4  4  4  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  2  0  0  0  0  0     0  0  0  3  9 15 14  9  4  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  3  9 17 18 13  6  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  3  7  8  6  3  2  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4  9 13 11  6  3  0  0  0  0     0  0  0  4 11 18 16  9  4  0  0  0  0     0  0  0  2  5 10 10  5  3  0  0  0  0
     0  0  2  9 28 49 50 31 14  6  0  0  0     0  0  2 11 39 72 76 49 22  9  0  0  0     0  0  0  6 22 43 47 31 14  6  0  0  0
     0  0  0  7 29 57 63 45 22  9  0  0  0     0  0  0 10 43 89100 71 35 15  5  0  0     0  0  0  6 26 56 65 47 23 10  0  0  0
     0  0  0  2  9 22 28 20 11  5  0  0  0     0  0  0  3 15 35 46 35 18  8  0  0  0     0  0  0  3  9 23 32 25 13  5  0  0  0
     0  0  0  0  0  3  4  4  3  0  0  0  0     0  0  0  0  2  5  8  6  4  0  0  0  0     0  0  0  0  0  4  6  5  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  4  4  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  5 11 12  8  4  4  0  0  0     0  0  0  0  3  3  4  3  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 14 18 13  7  4  0  0  0     0  0  0  0  2  4  4  3  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  6 10  8  4  3  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  2  3  3  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 STANDARD DEVIATION (DEGREES) OF BEAM DIVERGENCE
 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1    0.025  0.025  0.024  0.024  0.023  0.024  0.024  0.024  0.024
     2    0.024  0.025  0.025  0.024  0.023  0.024  0.024  0.024  0.024
     3    0.024  0.025  0.025  0.025  0.024  0.024  0.024  0.024  0.024
     4    0.024  0.025  0.025  0.024  0.024  0.025  0.024  0.024  0.024
     5    0.024  0.024  0.024  0.024  0.023  0.024  0.024  0.024  0.024
     6    0.024  0.025  0.025  0.024  0.024  0.025  0.024  0.024  0.024
     7    0.024  0.025  0.025  0.024  0.023  0.024  0.024  0.024  0.024
     8    0.024  0.025  0.025  0.024  0.023  0.024  0.024  0.024  0.024
     9    0.024  0.025  0.025  0.024  0.024  0.025  0.025  0.024  0.024
    10    0.024  0.025  0.025  0.025  0.024  0.024  0.024  0.024  0.024
    11    0.024  0.025  0.025  0.024  0.024  0.025  0.024  0.024  0.025
    12    0.025  0.025  0.026  0.025  0.024  0.025  0.025  0.024  0.024
    13    0.024  0.025  0.025  0.025  0.023  0.024  0.025  0.024  0.024
    14    0.025  0.025  0.025  0.025  0.024  0.025  0.025  0.025  0.025
    15    0.025  0.025  0.025  0.025  0.024  0.025  0.025  0.025  0.025
    16    0.024  0.025  0.025  0.025  0.024  0.025  0.025  0.024  0.024
    17    0.024  0.025  0.025  0.024  0.024  0.025  0.025  0.025  0.025
    18    0.024  0.025  0.025  0.024  0.024  0.025  0.025  0.025  0.025
    19    0.024  0.025  0.025  0.024  0.024  0.025  0.025  0.025  0.025
    20    0.024  0.025  0.026  0.024  0.024  0.025  0.025  0.025  0.024
    21    0.024  0.025  0.026  0.025  0.024  0.025  0.025  0.025  0.025
    22    0.024  0.025  0.026  0.025  0.024  0.025  0.025  0.025  0.025
    23    0.025  0.025  0.026  0.024  0.024  0.025  0.025  0.025  0.025
    24    0.024  0.025  0.026  0.025  0.024  0.025  0.025  0.025  0.025
    25    0.025  0.025  0.026  0.025  0.024  0.025  0.025  0.025  0.025


 STANDARD DEVIATION (DEGREES) OF REFLECTING RANGE
 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1    0.029  0.026  0.028  0.031  0.028  0.023  0.027  0.031  0.029
     2    0.030  0.031  0.030  0.030  0.029  0.025  0.028  0.031  0.030
     3    0.030  0.025  0.027  0.030  0.030  0.028  0.028  0.032  0.031
     4    0.029  0.024  0.028  0.030  0.030  0.029  0.027  0.030  0.030
     5    0.029  0.028  0.028  0.030  0.029  0.025  0.026  0.029  0.029
     6    0.030  0.026  0.028  0.031  0.029  0.023  0.026  0.031  0.030
     7    0.029  0.026  0.028  0.032  0.032  0.025  0.027  0.030  0.029
     8    0.031  0.042  0.034  0.030  0.028  0.026  0.029  0.033  0.038
     9    0.028  0.027  0.028  0.029  0.028  0.023  0.026  0.032  0.030
    10    0.030  0.026  0.029  0.031  0.029  0.028  0.027  0.030  0.029
    11    0.030  0.027  0.030  0.033  0.029  0.028  0.028  0.030  0.029
    12    0.029  0.025  0.029  0.035  0.032  0.024  0.027  0.031  0.030
    13    0.029  0.028  0.028  0.032  0.029  0.024  0.026  0.030  0.030
    14    0.030  0.025  0.027  0.031  0.030  0.027  0.028  0.031  0.029
    15    0.030  0.025  0.028  0.030  0.030  0.033  0.031  0.030  0.029
    16    0.029  0.027  0.028  0.030  0.030  0.030  0.028  0.030  0.029
    17    0.029  0.032  0.029  0.031  0.030  0.024  0.027  0.030  0.032
    18    0.030  0.027  0.028  0.030  0.028  0.025  0.029  0.034  0.032
    19    0.030  0.030  0.029  0.030  0.028  0.023  0.028  0.031  0.030
    20    0.030  0.040  0.034  0.031  0.029  0.027  0.028  0.030  0.035
    21    0.029  0.025  0.027  0.030  0.028  0.024  0.030  0.034  0.031
    22    0.032  0.037  0.032  0.031  0.029  0.024  0.028  0.031  0.033
    23    0.029  0.033  0.030  0.030  0.028  0.024  0.027  0.030  0.031
    24    0.031  0.026  0.027  0.029  0.030  0.035  0.033  0.031  0.029
    25    0.030  0.024  0.028  0.034  0.031  0.024  0.027  0.030  0.030


 ***** SUGGESTED VALUES FOR INPUT PARAMETERS *****
 BEAM_DIVERGENCE=   0.24530  BEAM_DIVERGENCE_E.S.D.=   0.02453
 REFLECTING_RANGE=  0.20420  REFLECTING_RANGE_E.S.D.=  0.02917



 NUMBER OF STRONG REFLECTIONS USED FOR LEARNING REFERENCE PROFILES

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1       49     33     60     66     53     27     52     65     57
     2       51     30     54     64     55     30     55     68     58
     3       52     29     57     63     56     34     54     72     57
     4       46     33     55     66     54     30     53     63     59
     5       47     27     53     62     54     27     54     60     51
     6       49     30     57     68     54     32     58     74     62
     7       52     36     57     62     49     27     57     72     64
     8       51     36     60     70     62     33     52     60     53
     9       48     35     57     61     50     29     53     72     63
    10       48     29     56     71     55     30     52     59     56
    11       55     36     61     68     64     33     59     73     60
    12       48     30     60     69     50     28     55     63     53
    13       51     26     47     67     57     30     55     72     59
    14       42     35     68     67     51     30     50     64     56
    15       46     31     52     66     54     30     53     67     59
    16       58     33     68     70     60     32     55     67     58
    17       47     30     51     67     49     27     56     65     54
    18       48     33     66     71     54     29     54     69     60
    19       50     33     47     67     64     33     52     63     58
    20       51     34     65     67     59     33     57     72     59
    21       50     29     57     71     57     31     52     65     55
    22       48     33     54     64     51     27     53     70     62
    23       49     35     62     68     55     30     60     74     62
    24       45     29     58     67     56     27     53     64     50
    25       47     32     57     67     56     31     56     67     54
 margin    1227    797   1440   1668   1376    751   1359   1676   1440



 NUMBER OF GRID POINTS IN REFERENCE PROFILES WITH  SIGNAL >   2.% OF MAXIMUM

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1      133    122    123    131    108    106    116    125    124
     2      139    126    132    126    111    104    116    126    125
     3      127    123    126    132    114    111    119    136    133
     4      123    120    129    128    111    111    115    127    125
     5      127    126    124    124    110    108    115    122    124
     6      132    123    130    130    114    108    121    132    126
     7      128    127    129    135    123    112    121    129    121
     8      142    135    141    120    110    109    121    134    139
     9      128    126    133    123    110    109    121    131    128
    10      134    131    138    135    113    108    116    130    128
    11      130    128    133    133    113    110    119    125    127
    12      133    125    137    134    118    103    118    138    126
    13      130    125    130    135    110    106    121    126    129
    14      135    124    130    137    118    109    122    138    131
    15      131    125    130    130    112    121    130    135    126
    16      135    125    126    128    110    118    119    129    126
    17      131    136    136    126    113    110    123    127    135
    18      128    128    130    126    112    111    125    147    134
    19      134    126    138    126    111    110    129    139    128
    20      132    135    142    125    111    116    125    135    139
    21      128    125    134    132    111    112    131    147    135
    22      138    130    149    130    111    113    121    133    138
    23      131    138    138    121    110    111    125    138    135
    24      139    131    129    126    110    120    137    140    134
    25      136    124    132    138    121    116    126    139    130
 margin     159    166    169    152    137    143    157    172    170

 MINIMUM INTEGRATION VOLUME IS  103 GRID POINTS



 CORRELATIONS BETWEEN REFERENCE PROFILES WITH CENTRAL PROFILE #1

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1     1.00   0.94   0.97   0.97   0.96   0.97   0.98   0.98   0.95
     2     1.00   0.84   0.93   0.98   0.95   0.93   0.96   0.99   0.93
     3     1.00   0.97   0.98   0.97   0.96   0.93   0.98   0.93   0.92
     4     1.00   0.98   0.96   0.99   0.94   0.87   0.97   0.97   0.97
     5     1.00   0.94   0.98   0.98   0.98   0.98   0.99   0.98   0.97
     6     1.00   0.96   0.96   0.97   0.97   0.95   0.98   0.98   0.95
     7     1.00   0.94   0.97   0.97   0.95   0.97   0.98   0.99   0.95
     8     1.00   0.78   0.95   0.99   0.96   0.97   0.98   0.98   0.87
     9     1.00   0.97   0.96   0.98   0.96   0.96   0.99   0.97   0.95
    10     1.00   0.97   0.96   0.98   0.97   0.89   0.96   0.98   0.94
    11     1.00   0.99   0.94   0.95   0.97   0.93   0.94   0.98   0.95
    12     1.00   0.97   0.94   0.95   0.94   0.96   0.98   0.98   0.94
    13     1.00   0.94   0.96   0.97   0.96   0.97   0.99   0.99   0.96
    14     1.00   0.98   0.98   0.98   0.97   0.91   0.96   0.98   0.96
    15     1.00   0.98   0.96   0.97   0.97   0.89   0.95   0.98   0.96
    16     1.00   0.96   0.97   0.98   0.97   0.93   0.98   0.98   0.95
    17     1.00   0.86   0.96   0.97   0.95   0.97   0.98   0.98   0.92
    18     1.00   0.97   0.95   0.98   0.97   0.96   0.97   0.95   0.94
    19     1.00   0.88   0.94   0.98   0.96   0.96   0.98   0.98   0.93
    20     1.00   0.84   0.94   0.98   0.97   0.96   0.99   0.99   0.94
    21     1.00   0.97   0.96   0.98   0.96   0.96   0.96   0.95   0.95
    22     1.00   0.80   0.90   0.99   0.96   0.95   0.96   0.97   0.94
    23     1.00   0.88   0.94   0.98   0.97   0.97   0.99   0.98   0.94
    24     1.00   0.97   0.95   0.98   0.97   0.83   0.93   0.98   0.95
    25     1.00   0.96   0.96   0.93   0.94   0.96   0.96   0.97   0.94
 margin    1.00   0.97   0.98   0.98   0.98   0.98   0.99   0.99   0.96



 CORRELATIONS BETWEEN REFERENCE PROFILES WITH RUN-AVERAGED PROFILES

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1     0.99   1.00   1.00   0.99   0.98   0.98   1.00   1.00   0.99
     2     1.00   0.95   0.98   0.99   0.97   0.90   0.96   0.99   0.99
     3     1.00   0.94   0.99   0.98   0.97   0.97   0.99   0.97   0.95
     4     1.00   0.98   0.99   0.98   0.95   0.93   0.99   0.99   0.98
     5     0.99   0.99   0.99   1.00   1.00   1.00   1.00   1.00   0.99
     6     1.00   0.99   0.99   1.00   1.00   0.99   1.00   1.00   0.99
     7     1.00   0.99   0.99   0.99   0.97   0.95   0.98   1.00   1.00
     8     0.98   0.79   0.92   0.98   0.99   0.99   0.94   0.92   0.84
     9     1.00   0.99   1.00   1.00   0.99   0.99   0.99   0.99   0.98
    10     1.00   0.95   0.98   0.99   0.99   0.95   0.99   1.00   0.97
    11     1.00   0.94   0.93   0.97   0.99   0.98   0.98   0.99   1.00
    12     0.97   0.95   0.93   0.94   0.92   0.91   0.96   0.99   0.98
    13     1.00   1.00   0.99   0.99   0.99   0.98   1.00   1.00   1.00
    14     0.99   0.97   0.98   0.99   0.99   0.98   0.99   1.00   0.99
    15     0.99   0.98   0.99   1.00   0.97   0.88   0.94   0.98   1.00
    16     0.99   1.00   0.98   0.98   0.96   0.95   0.99   1.00   1.00
    17     1.00   0.95   0.99   1.00   0.99   0.98   1.00   0.99   0.97
    18     0.99   0.99   1.00   0.99   0.99   0.99   0.98   0.97   0.96
    19     1.00   0.96   0.97   0.99   1.00   0.98   0.99   1.00   0.99
    20     0.98   0.87   0.93   0.98   0.99   0.99   0.98   0.97   0.92
    21     0.99   0.98   0.99   0.99   1.00   0.99   0.96   0.97   0.98
    22     0.99   0.89   0.96   0.99   0.99   0.96   0.98   1.00   0.97
    23     0.99   0.95   0.97   0.99   0.99   1.00   1.00   1.00   0.98
    24     0.99   0.97   0.99   0.99   0.95   0.83   0.90   0.99   1.00
    25     0.99   0.98   0.95   0.96   0.96   0.98   0.99   1.00   0.99
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

       -1... 0      -0.144      0.068        -1.          3
        0... 1       0.256      0.087         5.         32
        1... 2       0.451      0.082        15.         61
        2... 3       0.611      0.090        33.         86
        3... 6       0.708      0.091        82.        370
        6... 9       0.760      0.082       182.        518
        9...12       0.766      0.079       342.        681
       12...15       0.769      0.083       561.        824
       15...18       0.779      0.080       837.        939
       18...21       0.781      0.078      1169.        903
       21...24       0.783      0.083      1595.        922
       24...27       0.789      0.077      2150.        867
       27...30       0.784      0.082      2778.        835
       30...33       0.787      0.078      3622.        779
       33...36       0.791      0.081      4579.        739
       36...39       0.794      0.083      6097.        719
       39...42       0.792      0.080      8369.        727
       42...45       0.788      0.083     12088.        827
       45...48       0.792      0.080     23126.       1160
       48...51       0.804      0.101     54107.       1101

   13093 REFLECTIONS SAVED ON FILE "INTEGRATE.HKL"


 total elapsed wall-clock time for INTEGRATE      17.0 sec

 TASK     cpu time (sec)    elapsed wall-clock time (sec)
   1         273.0                    16.8
