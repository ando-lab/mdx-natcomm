
 ***** INTEGRATE ***** (VERSION Mar 15, 2019  BUILT=20190806)   6-Jan-2020

 INPUT PARAMETER VALUES
 ----------------------
 MAXIMUM_NUMBER_OF_PROCESSORS=   0  (0: automatic choice of OpenMP threads)
 number of OpenMP threads used  20
 MAXIMUM_NUMBER_OF_JOBS=   0  (0:automatic choice of forked main tasks)
 number of forked INTEGRATE tasks   1
 CLUSTER_NODES=
 SECONDS=     0
 NAME_TEMPLATE_OF_DATA_FRAMES=~/data/lysozyme/nitrate/ubatch4/lys_nitr_8_3_????.cbf CBF
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
 UNIT_CELL_CONSTANTS=    27.41    32.10    34.51  88.690 108.392 111.909
 ROTATION_AXIS= -1.00000 -0.00209  0.00055
 OSCILLATION_RANGE=  0.100000 DEGREES
 STARTING_ANGLE= 130.000 STARTING_FRAME=       1
 X-RAY_WAVELENGTH=   0.97700 ANGSTROM
 INCIDENT_BEAM_DIRECTION=  0.000708  0.005408  1.023527
 NUMBER OF DETECTOR SEGMENTS   1
 NX=  2463 NY=  2527 QX= 0.17200 QY= 0.17200
 ORGX=   1265.51  ORGY=   1237.09
 DETECTOR_DISTANCE=   231.319
 DIRECTION_OF_DETECTOR_X-AXIS=   1.00000   0.00000   0.00000
 DIRECTION_OF_DETECTOR_Y-AXIS=   0.00000   1.00000   0.00000
 SEGMENT=    1 2463    1 2527  ! x1,x2,y1,y2 (pixel units)
    SEGMENT_ORGX=      0.00   SEGMENT_ORGY=      0.00
    SEGMENT_DISTANCE=     0.00
    DIRECTION_OF_SEGMENT_X-AXIS=  1.000000  0.000000  0.000000
    DIRECTION_OF_SEGMENT_Y-AXIS=  0.000000  1.000000  0.000000
 RESOLUTION RANGE RECORDED BY DETECTOR (ANGSTROM)    44.586     1.212
 NUMBER OF TRUSTED DETECTOR PIXELS    5199248
 MEAN CONTENTS OF TRUSTED PIXELS IN BACKGROUND TABLE    126.371
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
     1   0  1.064  5192058    0    135      77     4  0.03847  0.07643
     2   0  1.025  5192353    0    140      74     1  0.03900  0.07799
     3   0  1.023  5191921    0    130      84     3  0.04200  0.08775
     4   0  1.025  5192284    0    122      77     1  0.04387  0.11387
     5   0  1.023  5191919    0    124      80     1  0.04507  0.12196
     6   0  1.025  5192688    0    117      67     3  0.04213  0.08066
     7   0  1.024  5192684    0    127      70     3  0.03827  0.08101
     8   0  1.022  5191693    0    127      77     2  0.04271  0.07988
     9   0  1.020  5191620    0    123      83     0  0.04174  0.09301
    10   0  1.020  5192429    0    137      74     2  0.03980  0.07299
    11   0  1.023  5191814    0    136      78     2  0.04207  0.07173
    12   0  1.024  5191612    0    146      79     1  0.03952  0.07722
    13   0  1.028  5191543    0    144      77     2  0.04093  0.07422
    14   0  1.025  5190548    0    135      94     2  0.04234  0.08025
    15   0  1.026  5191827    0    143      78     5  0.04135  0.07092
    16   0  1.024  5192105    0    134      77     4  0.04254  0.09237
    17   0  1.025  5192129    0    119      73     2  0.04006  0.09759
    18   0  1.028  5193369    0    119      62     3  0.04599  0.09318
    19   0  1.023  5192661    0    123      72     3  0.04448  0.09127
    20   0  1.026  5191305    0    127      83     0  0.04291  0.11220

   456 OUT OF    456 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.49
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.05
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.404    32.093    34.500  88.724 108.410 111.865
 REC. CELL PARAMETERS   0.041676  0.033773  0.030730  83.784  70.644  67.345
 COORDINATES OF UNIT CELL A-AXIS   -14.372   -23.312    -0.979
 COORDINATES OF UNIT CELL B-AXIS    31.380    -5.462     3.931
 COORDINATES OF UNIT CELL C-AXIS    -1.140    12.156    32.267
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.025    -0.016     0.001
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999998 -0.002087  0.000549
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000612  0.005154  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.25   1244.05
 DETECTOR ORIGIN (PIXELS) AT                     1265.44   1237.28
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.20
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     456 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   45     26     53     67     56     27     54     68     52
 SIGMAB (degree)       0.034  0.036  0.036  0.035  0.033  0.035  0.035  0.035  0.034
 SIGMAR (degree)       0.032  0.030  0.034  0.034  0.030  0.027  0.029  0.031  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    456 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  3  7  9  7  3  2  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  1  0  0  0     0  0  0  2  5 11 16 13  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  3  6 10  8  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  2  4  6  5  2  1  0  0  0     0  0  0  1  3  9 11  7  3  1  0  0  0     0  0  0  1  3  8  9  5  2  1  0  0  0
     0  0  0  2 10 25 34 26 12  5  1  0  0     0  0  0  4 20 46 58 41 18  8  2  0  0     0  0  0  4 19 36 41 27 11  4  1  0  0
     0  0  1  5 20 45 60 47 23 10  3  1  0     0  0  1  9 37 83100 73 34 14  5  1  0     0  0  1  8 32 64 70 47 21  8  3  1  0
     0  0  1  3 11 26 34 28 13  6  2  1  0     0  0  1  6 22 47 57 41 18  7  2  1  0     0  0  1  5 19 37 40 27 11  4  1  0  0
     0  0  0  0  1  3  5  4  2  1  0  0  0     0  0  0  1  3  7 10  6  3  1  0  0  0     0  0  0  1  3  6  7  4  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  3  3  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 11 12  7  3  1  0  0  0     0  0  0  1  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4 11 20 20 12  5  2  1  0  0     0  0  0  1  2  3  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 12 12  7  3  1  0  0  0     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  1  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    634 REFLECTION INTENSITIES INTEGRATED


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
    21   0  1.031  5191627    0    134      83     3  0.03927  0.09030
    22   0  1.030  5191450    0    131      85     2  0.04353  0.09735
    23   0  1.025  5191931    0    125      75     6  0.04339  0.06310
    24   0  1.027  5192542    0    119      74     4  0.04308  0.08723
    25   0  1.020  5191695    0    117      83     3  0.04772  0.12127
    26   0  1.025  5191741    0    116      80     2  0.04730  0.12924
    27   0  1.027  5191804    0    121      80     6  0.04576  0.08831
    28   0  1.025  5191761    0    133      80     7  0.04432  0.08357
    29   0  1.025  5190951    0    133      85     4  0.04038  0.08307
    30   0  1.025  5190467    0    134      92     1  0.04323  0.08200
    31   0  1.025  5191135    0    144      86     1  0.04650  0.09774
    32   0  1.024  5191393    0    135      83     2  0.04718  0.11238
    33   0  1.026  5191607    0    129      78     4  0.04229  0.08997
    34   0  1.024  5192086    0    121      75     3  0.04570  0.08768
    35   0  1.025  5193088    0    119      65     4  0.04366  0.07697
    36   0  1.022  5193168    0    112      66     2  0.04122  0.09029
    37   0  1.018  5193450    0    124      64     5  0.04469  0.06602
    38   0  1.021  5193440    0    119      67     1  0.04062  0.07976
    39   0  1.022  5193156    0    125      66     0  0.04201  0.09073
    40   0  1.020  5193316    0    132      61     1  0.03985  0.06906

   448 OUT OF    448 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.50
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.395    32.096    34.494  88.767 108.381 111.841
 REC. CELL PARAMETERS   0.041677  0.033766  0.030731  83.757  70.660  67.358
 COORDINATES OF UNIT CELL A-AXIS   -14.374   -23.302    -0.970
 COORDINATES OF UNIT CELL B-AXIS    31.382    -5.482     3.912
 COORDINATES OF UNIT CELL C-AXIS    -1.140    12.149    32.264
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.042    -0.013     0.012
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999998 -0.002087  0.000549
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000052  0.005107  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.24   1244.07
 DETECTOR ORIGIN (PIXELS) AT                     1266.17   1237.36
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.22
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     448 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     33     59     66     49     27     51     57     53
 SIGMAB (degree)       0.035  0.037  0.036  0.035  0.034  0.035  0.035  0.035  0.035
 SIGMAR (degree)       0.035  0.030  0.033  0.034  0.032  0.028  0.031  0.032  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    448 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  1     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  1  4  9 12 10  6  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  3  2  1  0  0  0     0  0  0  2  6 15 22 18 11  4  2  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  1  3  9 13 11  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0 -1  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  3  7  9  6  3  1  0  0  0     0  0  0  1  5 11 13  9  4  1  0  0  0     0  0  0  1  5  9  9  6  3  1  0  0  0
     0  0  0  3 14 29 38 29 14  6  2  0  0     0  0  0  5 23 50 59 39 17  7  2  0  0     0  0  0  5 20 38 41 25 11  4  1  0 -1
     0  0  1  6 21 53 70 51 26 11  4  1  0     0  0  1 10 40 86100 68 32 13  5  1  0     0  0  2 10 34 66 68 43 19  7  2  0  0
     0  0  1  4 14 33 44 32 17  7  2  0  0     0  0  1  7 26 54 60 41 19  8  2  0  0     0  0  1  7 21 40 38 25 11  4  1  0  0
     0  0  0  1  3  6  9  6  3  1  0  0  0     0  0  0  1  5 10 12  7  3  1  0  0  0     0  0  0  2  5  7  7  4  2  0  0  0  0
     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0
     0  0  0  1  2  3  3  2  1  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  2  7 12 12  7  3  1  0  0  0     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  1  4 11 20 20 12  5  2  0  0  0    -1  0  0  1  1  3  2  1  1  0  0  0  0     0  0  0  0  1  0  1  1  0  0  0  0  0
     0  0  1  3  7 12 11  7  3  1  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  2  2  2  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    602 REFLECTION INTENSITIES INTEGRATED


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
    41   0  1.022  5193372    0    134      65     2  0.03859  0.07777
    42   0  1.022  5193264    0    138      64     2  0.04155  0.07335
    43   0  1.021  5192495    0    130      70     4  0.04091  0.05361
    44   0  1.021  5191511    0    138      86     2  0.04216  0.07734
    45   0  1.021  5191791    0    147      77     3  0.04185  0.06944
    46   0  1.020  5191900    0    129      73     0  0.04087  0.08022
    47   0  1.019  5192216    0    131      77     0  0.04221  0.08072
    48   0  1.019  5192333    0    126      76     1  0.03795  0.08132
    49   0  1.018  5192664    0    124      76     3  0.04162  0.08451
    50   0  1.021  5192699    0    131      69     1  0.03822  0.07581
    51   0  1.020  5192384    0    137      72     2  0.03942  0.07470
    52   0  1.017  5192496    0    141      70     4  0.04084  0.05792
    53   0  1.019  5192699    0    138      73     2  0.04459  0.08501
    54   0  1.017  5193334    0    144      66     2  0.04372  0.08154
    55   0  1.016  5192485    0    146      78     3  0.03875  0.07533
    56   0  1.014  5191989    0    146      83     3  0.04252  0.06760
    57   0  1.013  5191971    0    160      79     2  0.04055  0.07241
    58   0  1.014  5191109    0    146      91     1  0.04198  0.08324
    59   0  1.014  5191700    0    144      81     2  0.04125  0.07668
    60   0  1.010  5192840    0    145      70     0  0.04358  0.09914

   498 OUT OF    498 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.54
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.05
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.417    32.094    34.497  88.677 108.432 111.894
 REC. CELL PARAMETERS   0.041666  0.033776  0.030734  83.817  70.637  67.328
 COORDINATES OF UNIT CELL A-AXIS   -14.374   -23.325    -0.993
 COORDINATES OF UNIT CELL B-AXIS    31.382    -5.440     3.950
 COORDINATES OF UNIT CELL C-AXIS    -1.140    12.149    32.266
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.012    -0.034    -0.011
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999998 -0.002087  0.000549
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.001035  0.005032  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.38   1244.04
 DETECTOR ORIGIN (PIXELS) AT                     1265.02   1237.43
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.23
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     498 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   51     29     53     67     58     29     59     76     61
 SIGMAB (degree)       0.035  0.036  0.036  0.035  0.033  0.034  0.035  0.035  0.034
 SIGMAR (degree)       0.032  0.031  0.034  0.033  0.029  0.025  0.029  0.032  0.032


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    498 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  1  3  7 10  8  4  2  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  1  6 13 19 15  8  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  1  0  0  0     0  0  0  1  3  8 11  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  2  4  6  4  2  1  0  0  0     0  0  0  1  4  8 10  7  3  1  0  0  0     0  0  0  1  3  7  8  5  2  1  0  0  0
     0  0  0  3 11 27 35 27 12  5  2  0  0     0  0  1  5 21 47 54 39 17  7  2  0  0     0  0  1  4 16 36 38 24 10  4  1  0  0
     0  0  1  5 21 49 66 49 24 10  3  1  0     0  0  1  9 37 82100 70 33 13  5  1  0     0  0  1  8 29 59 67 43 19  7  2  0  0
     0  0  1  3 13 28 37 27 14  5  2  0  0     0  0  1  6 22 47 57 39 19  7  2  0  0     0  0  1  5 19 35 37 24 11  4  1  0  0
     0  0  0  0  2  4  6  4  2  1  0  0  0     0  0  0  1  4  8  9  6  3  1  0  0  0     0  0  0  1  4  6  6  4  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  1  1  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 13 12  7  3  1  0  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4 11 20 19 12  5  2  1  0  0     0  0  0  1  2  3  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  3  7 12 11  7  3  1  0  0  0     0  0  0  1  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  2  2  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    656 REFLECTION INTENSITIES INTEGRATED


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
    61   0  1.013  5192447    0    146      72     0  0.04642  0.08409
    62   0  1.016  5191358    0    136      80     6  0.04099  0.05665
    63   0  1.013  5192152    0    136      76     2  0.03932  0.06124
    64   0  1.011  5192946    0    143      68     2  0.04263  0.07301
    65   0  1.010  5192018    0    136      76     2  0.04293  0.08470
    66   0  1.013  5191769    0    128      79     2  0.04327  0.08272
    67   0  1.011  5192245    0    140      74     5  0.04387  0.06787
    68   0  1.009  5192634    0    138      67     3  0.04093  0.06509
    69   0  1.007  5192268    0    137      70     2  0.04332  0.07418
    70   0  1.009  5192398    0    146      75     3  0.04119  0.06603
    71   0  1.007  5192236    0    140      75     2  0.03926  0.06835
    72   0  1.005  5192035    0    136      82     0  0.04025  0.07826
    73   0  1.007  5192769    0    135      73     1  0.03863  0.08165
    74   0  1.005  5193222    0    132      61     2  0.04409  0.07701
    75   0  1.010  5193173    0    132      69     3  0.04231  0.08474
    76   0  1.010  5192485    0    126      69     0  0.04044  0.09938
    77   0  1.009  5193147    0    127      63     0  0.04135  0.08674
    78   0  1.009  5193354    0    124      62     1  0.04029  0.07155
    79   0  1.011  5192802    0    125      70     3  0.04238  0.07161
    80   0  1.011  5192745    0    129      70     1  0.04409  0.09029

   463 OUT OF    463 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.51
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.07
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.415    32.102    34.489  88.670 108.395 111.906
 REC. CELL PARAMETERS   0.041661  0.033770  0.030734  83.837  70.678  67.322
 COORDINATES OF UNIT CELL A-AXIS   -14.379   -23.321    -0.982
 COORDINATES OF UNIT CELL B-AXIS    31.390    -5.441     3.955
 COORDINATES OF UNIT CELL C-AXIS    -1.142    12.140    32.261
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.007    -0.030    -0.009
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999998 -0.002087  0.000549
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.001016  0.005196  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.41   1244.08
 DETECTOR ORIGIN (PIXELS) AT                     1265.07   1237.25
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.30
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     463 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     34     64     72     51     26     50     61     54
 SIGMAB (degree)       0.035  0.036  0.037  0.035  0.033  0.035  0.035  0.035  0.034
 SIGMAR (degree)       0.031  0.028  0.033  0.033  0.030  0.024  0.028  0.031  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    463 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  1  5 11 15 13  7  3  1  0  0
     0  0  0  0  0  0  1  0  0  0  0  0 -1     0  0  0  0  1  3  5  4  3  1  1  0 -1     0  0  0  2  9 21 30 25 14  6  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  3  2  1  0  0  0     0  0  0  1  5 12 17 14  8  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  1  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  1  0  0  0     0  0  0  0  0  0  0  0  0  1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  3  6  8  6  3  1  0  0  0     0  0  0  1  4  9 11  8  3  1  0  0  0     0  0  0  1  3  6  6  4  2  0  0  0  0
     0  0  0  3 15 36 46 35 17  7  2  1  0     0  0  0  5 22 49 58 38 18  7  2  0  0     0  0  0  4 15 29 30 17  8  3  1  0  0
     0  0  1  7 29 68 84 64 33 13  5  1  0     0  0  1 11 42 88100 69 34 13  4  1  0     0  0  1  8 28 50 51 31 15  5  2  0  0
     0  0  1  4 17 38 49 37 18  7  3  1  0     0  0  1  7 26 51 59 42 18  7  3  0  0     0  0  1  6 18 29 30 19  7  3  1  0  0
     0  0  0  1  3  7  8  7  3  1  0  0  0     0  0  0  1  5 10 11  7  3  1  0  0  0     0  0  0  1  4  6  5  4  1  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0
     0  0  0  2  4  7  5  3  1  0  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  1  0  1  0  0  0  0  0  0
     0  0  0  3  7 11  9  5  2  1  0  0 -1     0  0  0  1  1  1  1  1  0  0  0  0 -1     0  0  0  0  1  0  1  0  0  0  0  0  0
     0  0  1  2  5  7  6  3  1  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    636 REFLECTION INTENSITIES INTEGRATED


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
    81   0  1.013  5192541    0    124      75     1  0.04340  0.08643
    82   0  1.009  5193530    0    117      63     0  0.04161  0.10406
    83   0  1.008  5193649    0    125      59     1  0.04198  0.09120
    84   0  1.006  5191736    0    128      77     0  0.04130  0.09294
    85   0  1.010  5191354    0    141      85     2  0.04241  0.07433
    86   0  1.006  5191666    0    146      83     2  0.04461  0.08730
    87   0  1.003  5191063    0    148      89     3  0.04501  0.08547
    88   0  1.003  5190472    0    152      95     3  0.04273  0.07575
    89   0  1.001  5190369    0    153      89     3  0.04449  0.08304
    90   0  1.005  5190811    0    150      80     1  0.04196  0.08487
    91   0  1.003  5190650    0    158      89     3  0.04180  0.07122
    92   0  1.005  5191416    0    157      89     1  0.04390  0.09548
    93   0  1.002  5191589    0    145      81     3  0.04696  0.08905
    94   0  1.005  5191382    0    142      86     3  0.04484  0.09009
    95   0  1.004  5191723    0    151      82     5  0.04350  0.08272
    96   0  1.004  5191317    0    146      84     3  0.04292  0.08603
    97   0  1.001  5191822    0    142      79     1  0.04316  0.09839
    98   0  1.002  5192428    0    146      74     2  0.04270  0.09243
    99   0  1.003  5192414    0    141      73     2  0.04078  0.09222
   100   0  1.000  5192556    0    147      77     0  0.04053  0.11644

   497 OUT OF    497 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.57
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.08
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.402    32.102    34.490  88.740 108.372 111.865
 REC. CELL PARAMETERS   0.041670  0.033764  0.030732  83.782  70.678  67.343
 COORDINATES OF UNIT CELL A-AXIS   -14.380   -23.305    -0.969
 COORDINATES OF UNIT CELL B-AXIS    31.387    -5.474     3.923
 COORDINATES OF UNIT CELL C-AXIS    -1.139    12.143    32.262
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.034    -0.015     0.012
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999998 -0.002087  0.000549
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000325  0.005116  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.26   1244.07
 DETECTOR ORIGIN (PIXELS) AT                     1265.84   1237.35
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.27
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     497 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   54     30     57     66     59     32     58     71     57
 SIGMAB (degree)       0.034  0.036  0.036  0.035  0.033  0.035  0.036  0.035  0.034
 SIGMAR (degree)       0.032  0.029  0.032  0.033  0.030  0.025  0.029  0.032  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    497 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  2  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0     0  0  0  1  4 10 14 12  6  3  1  0  0
     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  3  4  4  2  1  0  0  0     0  0  0  2  7 18 24 21 12  5  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0     0  0  0  1  4 11 15 12  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0 -1  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  2  5  7  6  3  1  0  0  0     0  0  0  1  4  7 10  8  3  1  0  0  0     0  0  0  0  3  4  5  3  1  0  0  0  0
     0  0  0  3 13 32 44 31 15  6  2  0  0     0  0  0  4 23 50 60 38 17  7  2  1  0     0  0  0  3 14 28 30 18  8  3  1  0  0
     0  0  1  5 22 58 74 57 29 12  4  1  1     0  0  1 10 39 89100 70 33 13  4  1  0     0  0  1  7 25 50 52 34 15  6  2  0  0
     0  0  1  4 15 36 44 33 17  7  2  0  0     0  0  1  8 27 54 58 39 19  7  2  1  0     0  0  1  6 17 31 30 18  9  3  1  0  0
     0  0  0  1  3  6  9  5  2  1  0  0  0     0  0  0  2  6 10 12  6  3  1  0  0  0     0  0  0  1  4  6  6  3  1  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0 -1  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  6  6  4  1  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 11 10  6  3  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4  7  6  4  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0 -1     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  1  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    668 REFLECTION INTENSITIES INTEGRATED


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
   101   0  1.003  5192415    0    141      72     0  0.04260  0.10464
   102   0  1.004  5191616    0    143      80     1  0.04217  0.08903
   103   0  1.002  5190990    0    142      86     2  0.04252  0.07319
   104   0  1.001  5191642    0    130      81     3  0.04334  0.07738
   105   0  0.999  5193292    0    130      66     3  0.04674  0.08840
   106   0  1.002  5193288    0    127      70     1  0.04833  0.09864
   107   0  0.999  5192746    0    126      70     2  0.04513  0.07719
   108   0  0.999  5191987    0    134      70     2  0.04029  0.06267
   109   0  0.998  5191389    0    147      81     5  0.03992  0.06480
   110   0  1.002  5190579    0    158      88     4  0.04443  0.07288
   111   0  0.998  5190229    0    149      95     4  0.04265  0.07462
   112   0  0.998  5190839    0    152      90     6  0.04455  0.07665
   113   0  1.001  5191231    0    147      86     5  0.04626  0.09585
   114   0  0.998  5191563    0    134      83     7  0.04297  0.07571
   115   0  1.002  5192931    0    137      70     2  0.04140  0.07944
   116   0  1.002  5192429    0    122      72     3  0.04388  0.07131
   117   0  1.000  5192824    0    124      66     1  0.04654  0.09176
   118   0  1.002  5192873    0    128      72     2  0.04682  0.09077
   119   0  1.001  5192252    0    124      71     6  0.03999  0.05747
   120   0  1.001  5191691    0    131      81     6  0.04129  0.06709

   461 OUT OF    462 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.53
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.02
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.413    32.104    34.501  88.699 108.406 111.894
 REC. CELL PARAMETERS   0.041666  0.033767  0.030727  83.803  70.655  67.323
 COORDINATES OF UNIT CELL A-AXIS   -14.383   -23.316    -0.987
 COORDINATES OF UNIT CELL B-AXIS    31.390    -5.455     3.948
 COORDINATES OF UNIT CELL C-AXIS    -1.146    12.149    32.270
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.008    -0.027     0.005
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999998 -0.002087  0.000549
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000766  0.005090  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.39   1244.11
 DETECTOR ORIGIN (PIXELS) AT                     1265.38   1237.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.30
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     461 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     33     57     66     55     28     50     61     55
 SIGMAB (degree)       0.035  0.037  0.037  0.035  0.033  0.034  0.035  0.035  0.035
 SIGMAR (degree)       0.030  0.027  0.032  0.033  0.029  0.025  0.028  0.030  0.028


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    461 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  2  1  1  0  0  0     0  0  0  1  3  8 11  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  3  2  1  0  0  0     0  0  0  2  6 14 19 16  9  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  1  0  0  0     0  0  0  1  3  8 11  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  1  2  5  7  5  3  1  0  0  0     0  0  0  1  4  9 12  8  4  1  0  0  0     0  0  0  1  3  6  7  4  2  1  0  0  0
     0  0  0  3 13 30 39 29 14  6  2  0  0     0  0  0  5 22 50 59 40 18  8  3  0  0     0  0  0  4 14 29 31 19  8  3  1  0  0
     0  0  1  6 23 53 67 51 26 10  4  1  0     0  0  1 10 41 86100 68 32 12  4  1  0     0  0  1  8 27 50 51 32 15  5  2  0  0
     0  0  1  4 14 32 40 28 14  6  2  0  0     0  0  1  7 26 53 60 37 17  7  2  0  0     0  0  1  6 18 31 31 18  8  3  1  0  0
     0  0  0  1  2  6  7  4  2  1  0  0  0     0  0  0  2  5 11 11  6  3  1  0  0  0     0  0  0  2  5  7  6  3  1  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4  6  6  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  1
     0  0  0  3  7 11 10  6  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  2  5  8  7  4  2  0  0  0  0     0  0  0  1  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    630 REFLECTION INTENSITIES INTEGRATED


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
   121   0  1.004  5192282    0    143      76     3  0.04425  0.08931
   122   0  1.001  5191830    0    143      84     3  0.04541  0.10129
   123   0  1.003  5192156    0    141      82     3  0.04123  0.08929
   124   0  1.005  5192260    0    133      78     5  0.04395  0.07879
   125   0  1.001  5192449    0    129      73     0  0.04399  0.10103
   126   0  1.001  5193125    0    135      64     1  0.04329  0.09028
   127   0  1.005  5192328    0    144      72     6  0.03935  0.06054
   128   0  1.000  5191434    0    141      82     6  0.03977  0.05565
   129   0  1.004  5193014    0    137      69     6  0.04682  0.07417
   130   0  1.001  5193403    0    131      67     3  0.04305  0.09394
   131   0  1.001  5193037    0    132      68     0  0.04361  0.12625
   132   0  1.002  5192411    0    130      78     2  0.04468  0.09681
   133   0  1.002  5192835    0    135      65     3  0.04290  0.06307
   134   0  1.003  5193086    0    140      66     4  0.04016  0.06943
   135   0  1.002  5192671    0    134      73     0  0.04010  0.08397
   136   0  1.002  5192598    0    127      68     2  0.03980  0.06846
   137   0  1.001  5193465    0    129      60     1  0.04135  0.08201
   138   0  1.000  5192811    0    131      72     7  0.04187  0.05789
   139   0  0.999  5192292    0    136      79     8  0.04459  0.05967
   140   0  0.999  5192657    0    132      70     1  0.04402  0.09535

   469 OUT OF    470 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.48
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.414    32.100    34.500  88.711 108.405 111.886
 REC. CELL PARAMETERS   0.041665  0.033771  0.030728  83.792  70.652  67.327
 COORDINATES OF UNIT CELL A-AXIS   -14.382   -23.317    -0.989
 COORDINATES OF UNIT CELL B-AXIS    31.386    -5.459     3.942
 COORDINATES OF UNIT CELL C-AXIS    -1.147    12.145    32.271
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.015    -0.035     0.008
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999998 -0.002087  0.000549
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000767  0.005132  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.36   1244.07
 DETECTOR ORIGIN (PIXELS) AT                     1265.35   1237.33
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.31
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     469 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   45     31     56     64     55     30     55     65     57
 SIGMAB (degree)       0.035  0.037  0.037  0.035  0.033  0.034  0.036  0.035  0.034
 SIGMAR (degree)       0.031  0.028  0.032  0.033  0.030  0.026  0.029  0.031  0.031


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    469 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  1  0  0  0     0  0  0  1  4  9 12 10  6  3  1  0  0
    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  3  2  1  0  0  0     0  0  0  2  7 16 22 19 10  5  2  1  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  1  4  9 13 10  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  3  6  8  6  3  2  0  0  0     0  0  0  1  4  9 12  8  4  2  0  0  0     0  0  0  1  3  7  8  4  2  0  0  0  0
     0  0  0  3 12 31 41 32 17  7  3  1  0     0  0  0  5 20 49 59 42 20  8  3  1  0     0  0  0  4 16 31 34 21  9  3  1  0  0
     0  0  1  5 23 56 70 54 28 12  4  1  0     0  0  1 10 41 90100 70 34 14  5  1  0     0  0  1  9 30 57 55 34 16  6  2  0  0
     0  0  1  4 14 32 40 29 15  6  2  0  0     0  0  1  8 27 54 58 37 18  7  2  0  0     0  0  1  8 21 35 32 19  9  3  1  0  0
     0  0  0  1  2  5  6  4  1  1  0  0  0     0  0  0  2  5 10  9  6  2  1  0  0  0     0  0  0  2  5  7  6  3  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0 -1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  8  8  4  2  1  0  0  0     0  0  0  1  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0 -1
     0  0  1  4  9 15 13  7  3  1  0  0  0     0  0  0  1  1  2  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  3  7 10  8  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  3  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  1  0  0  0  0  0

    635 REFLECTION INTENSITIES INTEGRATED


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
   141   0  1.001  5192954    0    134      68     3  0.04100  0.07141
   142   0  1.002  5191860    0    140      83     1  0.03913  0.07822
   143   0  0.999  5192239    0    137      73     2  0.04037  0.06333
   144   0  0.995  5192841    0    130      69     1  0.03984  0.07894
   145   0  0.995  5193629    0    114      64     0  0.03897  0.08584
   146   0  0.994  5192907    0    119      74     1  0.04224  0.08496
   147   0  0.993  5191999    0    123      84     0  0.04140  0.09115
   148   0  0.991  5192537    0    133      72     1  0.04196  0.07147
   149   0  0.993  5192937    0    136      71     1  0.04036  0.07375
   150   0  0.990  5192316    0    130      73     1  0.04435  0.07792
   151   0  0.991  5192221    0    133      74     2  0.04469  0.07618
   152   0  0.993  5191923    0    134      78     0  0.04026  0.07448
   153   0  0.992  5191374    0    140      83     2  0.03997  0.06821
   154   0  0.990  5190739    0    141      93     4  0.04464  0.07437
   155   0  0.989  5191019    0    157      85     3  0.04664  0.08586
   156   0  0.994  5191255    0    148      87     3  0.04336  0.08691
   157   0  0.993  5192210    0    148      75     1  0.04285  0.09948
   158   0  0.993  5191981    0    140      76     1  0.04111  0.09007
   159   0  0.993  5192318    0    122      73     0  0.04628  0.09742
   160   0  0.995  5192238    0    133      73     0  0.04246  0.08744

   477 OUT OF    478 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.52
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.03
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.412    32.097    34.484  88.686 108.411 111.901
 REC. CELL PARAMETERS   0.041671  0.033776  0.030742  83.814  70.654  67.319
 COORDINATES OF UNIT CELL A-AXIS   -14.381   -23.316    -0.989
 COORDINATES OF UNIT CELL B-AXIS    31.384    -5.449     3.951
 COORDINATES OF UNIT CELL C-AXIS    -1.143    12.142    32.256
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.009    -0.033     0.003
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999998 -0.002087  0.000549
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000840  0.005047  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.41   1244.08
 DETECTOR ORIGIN (PIXELS) AT                     1265.31   1237.45
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.25
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     477 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   50     33     57     70     50     26     52     69     61
 SIGMAB (degree)       0.035  0.037  0.038  0.035  0.033  0.035  0.036  0.035  0.035
 SIGMAR (degree)       0.032  0.028  0.032  0.033  0.030  0.026  0.029  0.031  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    477 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  2  1  1  0  0  0     0  0  0  1  3  9 12 11  7  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  1  1  2  3  3  2  1  1  0  0     0  0  0  1  6 16 23 21 12  6  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  1  0  1  2  2  1  1  0  0  0     0  0  0  1  4  9 13 11  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0 -1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  1  2  5  7  5  3  1  0  0  0     0  0  0  1  3  8 11  8  4  2  0  0  0     0  0  0  1  3  6  7  5  2  1  0  0  0
     0  0  0  3 13 32 41 32 17  8  3  0  0     0  0  0  4 22 47 56 40 20  8  3  0  0     0  0  0  3 15 28 31 20  9  4  1  0  0
     0  0  0  5 24 57 74 57 30 12  5  1  0     0  0  1 10 41 87100 69 35 13  5  1  0     0  0  1  7 27 51 52 33 16  6  2  0  0
     0  0  1  4 14 32 42 30 15  6  2  0  0     0  0  1  9 26 53 60 39 17  7  2  0  0     0  0  1  8 19 31 33 19  8  3  1  0  0
     0  0  0  1  3  6  8  5  2  1  0  0  0     0  0  0  2  6 12 12  7  3  1  0  0  0     0  0  0  3  5  8  7  3  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  1  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0 -1  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  7  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  7 12 11  7  3  1  0  0  0     0  0  0  1  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  5  7  7  4  1  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0 -1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    632 REFLECTION INTENSITIES INTEGRATED


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
   161   0  0.992  5192528    0    130      69     0  0.04219  0.09146
   162   0  0.994  5191944    0    135      79     6  0.04439  0.06935
   163   0  0.996  5192192    0    137      77     2  0.04206  0.08015
   164   0  0.993  5191622    0    137      81     3  0.04475  0.08716
   165   0  0.993  5191874    0    136      74     3  0.04450  0.07657
   166   0  0.992  5192409    0    137      67     2  0.04572  0.10354
   167   0  0.993  5192220    0    132      71     2  0.04696  0.09055
   168   0  0.990  5191929    0    118      71     1  0.04526  0.08977
   169   0  0.991  5192248    0    133      69     0  0.04539  0.10259
   170   0  0.995  5192625    0    134      66     0  0.04468  0.08518
   171   0  0.995  5191619    0    134      78     1  0.04249  0.07999
   172   0  0.995  5192008    0    139      76     2  0.04026  0.07119
   173   0  0.996  5192628    0    145      72     1  0.04360  0.08164
   174   0  0.993  5191740    0    145      83     2  0.04249  0.07597
   175   0  0.990  5191901    0    152      82     1  0.04383  0.08770
   176   0  0.989  5191868    0    141      82     1  0.04332  0.09399
   177   0  0.992  5191798    0    131      80     3  0.04165  0.08079
   178   0  0.993  5192710    0    131      69     2  0.04549  0.07580
   179   0  0.990  5193700    0    120      57     1  0.04490  0.08448
   180   0  0.991  5193220    0    112      60     1  0.04362  0.08062

   456 OUT OF    456 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.59
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.08
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.433    32.101    34.523  88.645 108.473 111.914
 REC. CELL PARAMETERS   0.041656  0.033774  0.030718  83.827  70.602  67.312
 COORDINATES OF UNIT CELL A-AXIS   -14.383   -23.338    -1.012
 COORDINATES OF UNIT CELL B-AXIS    31.388    -5.433     3.972
 COORDINATES OF UNIT CELL C-AXIS    -1.146    12.165    32.288
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.010    -0.040    -0.006
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999998 -0.002087  0.000549
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.001412  0.005091  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.50   1244.01
 DETECTOR ORIGIN (PIXELS) AT                     1264.64   1237.32
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.37
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     456 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   43     27     54     66     55     31     55     64     51
 SIGMAB (degree)       0.035  0.037  0.038  0.036  0.033  0.035  0.036  0.036  0.035
 SIGMAR (degree)       0.033  0.030  0.034  0.034  0.029  0.026  0.029  0.031  0.031


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    456 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  1  0  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  1  1  2  1  1  0  0  0     0  0  0  1  3  8 11 10  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  3  2  1  0  0  1     0  0  0  1  5 13 18 18 11  5  2  1  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  1  3  6 10  8  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  1  0  0  0  0 -1  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0
     0  0  0  0  2  5  7  6  4  1  0  0  0     0  0  0  1  4  8 10  8  4  2  0  0  0     0  0  0  1  3  7  8  6  3  1  0  0  0
     0  0  0  3 12 31 40 33 18  9  3  1  0     0  0  1  6 22 50 59 42 21 10  3  1  0     0  0  1  5 18 36 40 25 12  6  2  0 -1
     0  0  1  5 22 52 66 56 30 14  5  1  0     0  0  2 11 41 85100 74 38 16  6  1  0     0  0  2 10 34 59 64 43 21  8  3  0  0
     0  0  0  3 13 29 38 28 16  6  2  0  0     0  0  1  8 27 51 59 39 20  8  3  0  0     0  0  1  7 23 38 37 23 11  4  1  0 -1
     0  0  0  1  2  5  5  4  2  1  0  0  0     0  0  0  2  5  9  9  5  2  1  0  0  0     0  0  0  2  6  8  7  4  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  1  0  0  0  0  0 -1  0     0  0  0  0  1  1  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0 -1 -1 -1  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0 -1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0 -1  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  4  4  2  1  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  2  8 13 13  8  4  1  0 -1  0     0  0  0  1  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     1  0  1  6 14 20 19 11  5  2  1  0  0     0  0  0  1  2  3  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  2  5 12 16 12  7  3  1  0  0  0     0  0  1  1  2  2  2  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  1  2  5  5  3  1  1  0  0  0  0     0  0  1  1  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0 -1  0 -1  0  0  0  0  0  0  0     0 -1 -1 -1  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    615 REFLECTION INTENSITIES INTEGRATED


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
   181   0  0.991  5193227    0    124      64     2  0.04435  0.08576
   182   0  0.994  5192473    0    137      70     0  0.04049  0.08840
   183   0  0.995  5191560    0    133      84     0  0.03958  0.08879
   184   0  0.993  5191853    0    143      81     2  0.03727  0.07274
   185   0  0.994  5192159    0    142      75     3  0.04082  0.06795
   186   0  0.994  5192458    0    132      70     2  0.04747  0.07767
   187   0  0.995  5193595    0    136      57     0  0.04354  0.07600
   188   0  0.994  5192704    0    133      64     0  0.04193  0.07072
   189   0  0.995  5192412    0    136      70     1  0.03847  0.06208
   190   0  0.993  5192378    0    143      74     1  0.03870  0.06799
   191   0  0.991  5192093    0    135      70     1  0.04121  0.07537
   192   0  0.992  5192114    0    128      69     6  0.04292  0.06020
   193   0  0.994  5192292    0    125      72     2  0.04160  0.08376
   194   0  0.995  5192930    0    130      62     1  0.04804  0.08979
   195   0  0.995  5192634    0    134      66     1  0.04544  0.09260
   196   0  0.993  5191745    0    131      75     0  0.04067  0.06744
   197   0  0.996  5191273    0    144      82     2  0.04068  0.06851
   198   0  0.996  5191323    0    135      81     0  0.03861  0.09069
   199   0  0.998  5191378    0    131      77     3  0.04257  0.06739
   200   0  0.991  5192618    0    137      72     1  0.04359  0.09761

   469 OUT OF    470 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.52
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.05
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.417    32.100    34.499  88.699 108.413 111.896
 REC. CELL PARAMETERS   0.041664  0.033772  0.030730  83.799  70.647  67.319
 COORDINATES OF UNIT CELL A-AXIS   -14.382   -23.321    -0.985
 COORDINATES OF UNIT CELL B-AXIS    31.387    -5.450     3.945
 COORDINATES OF UNIT CELL C-AXIS    -1.144    12.153    32.267
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.010    -0.024     0.002
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999998 -0.002087  0.000549
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000815  0.005160  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.40   1244.06
 DETECTOR ORIGIN (PIXELS) AT                     1265.33   1237.28
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.33
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     469 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   48     33     60     65     54     28     53     67     56
 SIGMAB (degree)       0.035  0.038  0.038  0.036  0.033  0.035  0.036  0.035  0.034
 SIGMAR (degree)       0.030  0.027  0.032  0.032  0.029  0.027  0.029  0.031  0.029


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    469 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  1  0  0  0     0  0  0 -1  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  1  2  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  2  1  1  1  0  0     0  0  0  1  4  9 13 11  7  4  2  1  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  3  3  2  1  1  0  0     0  0  0  2  6 15 22 20 13  6  3  1  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  1  4  9 12 11  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0 -1  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  1  2  5  8  6  4  2  1  0  0     0  0  0  1  3  8 11  8  4  2  1  0  0     0  0  0  1  3  5  6  4  2  1  0  0  0
     0  0  0  3 13 31 42 34 18  8  3  1  0     0  0  0  5 22 47 56 41 20  9  3  1  0     0  0  0  4 15 27 29 18  9  3  1  0  0
     0  0  1  6 23 54 73 58 32 13  5  1  0     0  0  1 12 42 87100 70 35 14  5  1  0     0  0  2 10 29 50 51 33 15  5  2  0  0
     0  0  1  4 15 34 41 31 16  7  2  0  0     0  0  2 10 31 58 59 39 18  7  2  0 -1     0  0  2  9 23 35 30 18  8  3  1  0  0
     0  0  0  1  3  6  7  5  2  1  0  0  0     0  0  0  3  8 11 12  6  3  1  0 -1  0     0  0  0  3  7  8  6  3  1  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0 -1  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  6  6  4  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  3  7 11 10  6  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  1     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  3  6  9  7  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    620 REFLECTION INTENSITIES INTEGRATED


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
   201   0  0.996  5192328    0    132      72     4  0.04211  0.07680
   202   0  0.994  5192491    0    132      71     1  0.04191  0.07412
   203   0  0.996  5192342    0    138      70     2  0.04319  0.06481
   204   0  0.996  5192281    0    154      68     0  0.04143  0.06992
   205   0  0.998  5191377    0    143      77     1  0.04070  0.06697
   206   0  0.996  5191134    0    150      86     3  0.04062  0.06629
   207   0  0.999  5191395    0    152      91     0  0.03971  0.09404
   208   0  0.999  5191926    0    152      77     0  0.04474  0.07878
   209   0  0.999  5191850    0    149      79     1  0.04143  0.06799
   210   0  1.001  5191828    0    150      76     0  0.04282  0.06907
   211   0  1.001  5191471    0    141      81     2  0.04241  0.07139
   212   0  1.003  5191752    0    148      78     4  0.04167  0.05935
   213   0  1.003  5191012    0    145      87     1  0.04059  0.07230
   214   0  1.004  5190299    0    143      90     1  0.04213  0.06520
   215   0  1.005  5191146    0    144      86     0  0.04234  0.08063
   216   0  1.006  5191285    0    150      82     2  0.04166  0.06289
   217   0  1.005  5191215    0    153      79     1  0.04382  0.07607
   218   0  1.004  5191248    0    150      88     2  0.04487  0.08817
   219   0  1.006  5191727    0    134      80     1  0.04187  0.07716
   220   0  1.006  5192908    0    125      67     2  0.04109  0.07014

   493 OUT OF    493 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.50
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.06
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.407    32.095    34.490  88.758 108.385 111.873
 REC. CELL PARAMETERS   0.041669  0.033775  0.030736  83.753  70.656  67.327
 COORDINATES OF UNIT CELL A-AXIS   -14.382   -23.311    -0.973
 COORDINATES OF UNIT CELL B-AXIS    31.383    -5.467     3.913
 COORDINATES OF UNIT CELL C-AXIS    -1.142    12.148    32.260
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.036    -0.019     0.011
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999998 -0.002087  0.000549
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000352  0.005119  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.27   1244.04
 DETECTOR ORIGIN (PIXELS) AT                     1265.80   1237.32
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.30
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     493 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   51     35     60     73     60     31     54     65     59
 SIGMAB (degree)       0.035  0.038  0.039  0.036  0.034  0.035  0.037  0.036  0.035
 SIGMAR (degree)       0.030  0.026  0.031  0.033  0.030  0.026  0.029  0.031  0.029


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    493 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  2  1  0  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  1  2  2  2  1  0  0  0     0  0  0  1  4  9 14 13  9  5  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  3  3  2  1  0  0     0  0  0  2  6 15 22 21 14  7  3  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  1  0  0  0     0  0  0  1  4  9 13 12  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  1  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0
     0  0  0  1  2  5  9  7  4  2  1  0  0     0  0  0  1  4  9 12  9  4  2  1  0  0     0  0  0  1  3  5  6  4  2  1  0  0  0
     0  0  0  4 13 32 43 36 20  9  3  1  0     0  0  1  6 23 51 60 43 22  8  3  1  0     0  0  0  5 15 28 32 20 10  3  1  0  0
     0  0  1  7 24 57 72 58 33 14  6  2  0     0  0  2 13 44 91100 68 35 14  5  1  0     0  0  2 10 30 52 52 32 16  6  2  0  0
     0  0  1  4 17 36 43 33 17  7  3  1  0     0  0  2 10 33 59 61 40 18  7  2  0  0     0  0  2  9 24 35 33 19  8  3  1  0  0
     0  0  0  1  4  8  8  6  3  1  0  0  0     0  0  1  3  9 13 12  7  3  1  0  0  0     0  0  1  3  7  9  7  4  1  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0 -1  0     0  0  0  1  1  1  1  1  0  0  0 -1  0     0  0  0  0  1  1  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0 -1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4  6  6  3  2  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  7 10  9  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  3  6  7  6  3  1  0  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0 -1  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    636 REFLECTION INTENSITIES INTEGRATED


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
   221   0  1.006  5193259    0    121      69     2  0.04048  0.08094
   222   0  1.010  5193625    0    124      64     1  0.04509  0.08681
   223   0  1.009  5192098    0    127      77     0  0.04074  0.07313
   224   0  1.011  5192093    0    128      78     1  0.03745  0.06670
   225   0  1.012  5192760    0    131      71     1  0.04110  0.07548
   226   0  1.011  5193372    0    138      62     2  0.04205  0.08060
   227   0  1.013  5192508    0    137      70     3  0.04056  0.08009
   228   0  1.011  5191844    0    144      81     3  0.04124  0.07847
   229   0  1.009  5191566    0    153      82     5  0.04068  0.07258
   230   0  1.012  5191640    0    141      80     7  0.04206  0.05323
   231   0  1.017  5192342    0    141      76     1  0.04289  0.09113
   232   0  1.017  5193030    0    127      64     2  0.04304  0.09533
   233   0  1.020  5192624    0    128      70     1  0.04103  0.08911
   234   0  1.019  5192094    0    139      71     0  0.04301  0.08411
   235   0  1.020  5192347    0    145      72     1  0.04252  0.08791
   236   0  1.021  5192357    0    135      71     2  0.04229  0.08932
   237   0  1.021  5192527    0    138      70     3  0.04056  0.07751
   238   0  1.023  5192456    0    130      74     7  0.03991  0.05973
   239   0  1.026  5193229    0    126      65     3  0.03749  0.06258
   240   0  1.030  5192745    0    120      69     6  0.03893  0.05404

   462 OUT OF    464 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.52
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.412    32.093    34.495  88.716 108.415 111.887
 REC. CELL PARAMETERS   0.041671  0.033779  0.030735  83.782  70.639  67.323
 COORDINATES OF UNIT CELL A-AXIS   -14.378   -23.317    -0.984
 COORDINATES OF UNIT CELL B-AXIS    31.381    -5.452     3.934
 COORDINATES OF UNIT CELL C-AXIS    -1.143    12.153    32.263
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.019    -0.024     0.001
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999998 -0.002087  0.000549
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000738  0.005064  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.33   1244.07
 DETECTOR ORIGIN (PIXELS) AT                     1265.37   1237.41
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.30
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     462 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     31     53     64     51     26     52     68     60
 SIGMAB (degree)       0.035  0.038  0.039  0.037  0.034  0.035  0.037  0.036  0.034
 SIGMAR (degree)       0.030  0.027  0.031  0.032  0.029  0.027  0.030  0.031  0.029


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    462 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  2  1  1  0  0  0     0  0  0  1  3  8 12 12  8  4  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  3  2  1  1  0  0     0  0  0  1  6 13 20 18 12  6  3  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  2  1  1  0  0  0     0  0  0  1  3  7 11 10  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  1  3  6  7  6  4  2  1  0  0     0  0  0  1  5  9 11  8  4  2  0  0  0     0  0  0  1  3  6  7  5  2  1  0  0  0
     0  0  0  3 14 30 39 34 18  8  3  1  0     0  0  1  6 24 49 59 45 21  9  3  1  0     0  0  1  5 17 30 34 24 10  4  1  0  0
     0  0  1  6 23 50 66 53 30 13  5  2  0     0  0  2 14 44 89100 70 36 14  5  1  0     0  0  2 12 32 54 54 36 16  6  2  0  0
     0  0  1  4 13 28 37 27 15  7  2  1  0     0  0  2 11 30 56 59 36 18  7  2  0  0     0  0  2 10 24 35 33 17  8  3  1  0  0
     0  0  0  1  3  5  5  4  2  1  0  0  0     0  0  1  3  8 11 10  7  3  1  0  0  0     0  0  1  3  7  7  6  3  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  8  5  2  1  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  3  8 12 11  7  3  1  0  0  0     0  0  0  1  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  3  6  8  7  3  2  1  0  0  0     0  0  0  1  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    614 REFLECTION INTENSITIES INTEGRATED


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
   241   0  1.029  5192809    0    116      68     6  0.04251  0.06350
   242   0  1.032  5193843    0    140      57     6  0.04037  0.07000
   243   0  1.031  5192107    0    148      81     5  0.04200  0.07733
   244   0  1.033  5191004    0    148      90     5  0.04227  0.07159
   245   0  1.034  5191054    0    159      91     5  0.04426  0.07030
   246   0  1.035  5191485    0    153      82     6  0.04182  0.06221
   247   0  1.033  5191059    0    142      88     5  0.03999  0.06133
   248   0  1.040  5191391    0    136      83     4  0.04186  0.07259
   249   0  1.037  5192820    0    128      66     1  0.04362  0.09823
   250   0  1.042  5192821    0    119      68     4  0.04354  0.08715
   251   0  1.041  5192855    0    118      64     1  0.04392  0.08876
   252   0  1.045  5192521    0    134      69     4  0.04247  0.07555
   253   0  1.048  5191661    0    136      74     4  0.04102  0.07467
   254   0  1.045  5190925    0    136      86     2  0.03950  0.07909
   255   0  1.049  5191556    0    144      80     3  0.04359  0.08299
   256   0  1.049  5191794    0    141      75     5  0.04241  0.07230
   257   0  1.047  5191376    0    150      82     6  0.04290  0.07314
   258   0  1.049  5191946    0    146      76     1  0.04039  0.09063
   259   0  1.053  5191361    0    131      84     6  0.03955  0.05440
   260   0  1.060  5191718    0    131      79     6  0.04002  0.06070

   489 OUT OF    491 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.56
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.11
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.410    32.093    34.509  88.750 108.409 111.878
 REC. CELL PARAMETERS   0.041673  0.033778  0.030723  83.750  70.633  67.322
 COORDINATES OF UNIT CELL A-AXIS   -14.379   -23.315    -0.980
 COORDINATES OF UNIT CELL B-AXIS    31.382    -5.459     3.918
 COORDINATES OF UNIT CELL C-AXIS    -1.144    12.160    32.275
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.029    -0.018     0.005
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999998 -0.002087  0.000549
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000539  0.005070  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.25   1244.05
 DETECTOR ORIGIN (PIXELS) AT                     1265.54   1237.38
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.29
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     489 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   49     30     55     70     60     29     54     69     61
 SIGMAB (degree)       0.035  0.038  0.039  0.037  0.034  0.036  0.037  0.037  0.035
 SIGMAR (degree)       0.030  0.027  0.030  0.032  0.030  0.028  0.031  0.031  0.029


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    489 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0 -1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  1  1  0  0  0     0  0  0  1  3  6  9  9  7  4  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  2  1  1  0  0     0  0  0  1  4 10 16 15 10  6  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  1  0  0  0  0     0  0  0  1  2  6  9  8  6  3  1  0 -1
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0 -1  0  0  0  0  0  0     0  0  0  0  0  0 -1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0 -1  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  1  0  0 -1  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  2  5  8  6  4  2  1  0  0     0  0  0  1  4  7 12  9  5  2  1  0  0     0  0  0  1  3  5  8  6  3  1  0  0  0
     0  0  0  3 11 25 36 29 18  9  4  1  0     0  0  1  6 22 47 59 41 23 10  4  1  0     0  0  1  5 17 33 38 25 13  5  1  0  0
     0  0  1  5 20 45 59 48 28 14  6  2  0     0  0  2 12 41 85100 71 36 16  6  2  0     0  0  2 12 36 62 64 41 19  8  3  1  0
     0  0  1  3 12 27 34 28 16  7  3  1  0     0  0  1  8 28 54 58 40 20  9  3  1  0     0  0  2 10 28 42 40 24 11  4  1  0  0
     0  0  0  1  2  5  5  4  2  1  0 -1  0     0  0  0  2  7 11  9  6  3  1  0  0  0     0  0  1  3  7  9  7  4  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  1  1  1  0  0  0  0 -1  0     0  0  0  0  1  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0     0  0 -1  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  8  9  6  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  5 11 16 15  9  4  2  0  0  0     0  0  0  1  1  2  1  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  4  9 12 10  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  3  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0 -1  0  0  0  0  0  0  0  0  0  0

    645 REFLECTION INTENSITIES INTEGRATED


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
   261   0  1.063  5193377    0    136      62     1  0.04116  0.08719
   262   0  1.062  5193454    0    125      59     2  0.04183  0.08734
   263   0  1.068  5193179    0    118      62     5  0.04321  0.07660
   264   0  1.070  5193650    0    124      57     5  0.04314  0.06018
   265   0  1.067  5192363    0    123      74     7  0.04528  0.06040
   266   0  1.070  5192088    0    143      77     8  0.04523  0.06335
   267   0  1.072  5191790    0    144      78     4  0.04444  0.07178
   268   0  1.073  5191716    0    140      80     5  0.04194  0.06238
   269   0  1.072  5191960    0    140      75     2  0.04466  0.08336
   270   0  1.076  5191550    0    140      73     3  0.04033  0.05783
   271   0  1.077  5191182    0    143      83     4  0.04260  0.06119
   272   0  1.078  5190442    0    154      91     4  0.04130  0.07030
   273   0  1.078  5190639    0    154      92     2  0.04424  0.08618
   274   0  1.079  5191491    0    144      88     9  0.03982  0.06181
   275   0  1.082  5192445    0    129      77     5  0.03966  0.07521
   276   0  1.084  5192447    0    125      74     4  0.04168  0.07636
   277   0  1.082  5192292    0    127      74     2  0.04052  0.08268
   278   0  1.085  5192087    0    130      78     2  0.04020  0.07262
   279   0  1.084  5193418    0    128      67     1  0.04146  0.08557
   280   0  1.089  5193180    0    130      69     0  0.04586  0.10082

   472 OUT OF    472 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.56
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.411    32.094    34.495  88.728 108.413 111.886
 REC. CELL PARAMETERS   0.041672  0.033778  0.030736  83.770  70.637  67.320
 COORDINATES OF UNIT CELL A-AXIS   -14.379   -23.316    -0.986
 COORDINATES OF UNIT CELL B-AXIS    31.382    -5.455     3.930
 COORDINATES OF UNIT CELL C-AXIS    -1.145    12.150    32.264
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.022    -0.029     0.003
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999998 -0.002087  0.000549
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000720  0.005099  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.31   1244.02
 DETECTOR ORIGIN (PIXELS) AT                     1265.37   1237.32
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.29
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     472 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   50     32     59     67     52     30     52     66     55
 SIGMAB (degree)       0.035  0.038  0.039  0.037  0.034  0.035  0.036  0.036  0.035
 SIGMAR (degree)       0.033  0.027  0.030  0.032  0.029  0.028  0.031  0.031  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    472 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  1  0  0  0  0     0  0  0  0  1  2  2  3  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  2  1  1  0  0     0  0  0  1  4  8 11 11  7  4  2  1  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  3  3  2  2  1  0  0     0  0  0  1  6 13 19 19 12  6  3  1  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     1  0  0  0  1  1  1  2  1  1  0  0  0     0  0  0  1  3  8 11  9  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0
     0  0  0  1  3  7  8  7  4  2  1  0  0     0  0  0  1  5 11 12  9  5  2  1  0  0     0  0  0  1  3  7  7  5  2  1  0  0  0
     0  0  0  4 15 30 40 32 19  8  3  1  0     0  0  1  6 26 52 60 42 22  9  4  1  0     0  0  1  5 18 31 32 20 10  4  1  0  0
     0  0  1  6 26 52 68 56 31 14  6  2  0     0  0  2 12 48 89100 70 35 15  6  1  0     0  0  2 12 37 56 54 34 16  7  2  0  0
     0  0  1  4 14 31 38 29 15  8  3  0  0     0  0  2 10 32 57 59 38 18  8  2  0  0     0  0  2 11 27 38 33 19  8  3  1  0 -1
     0  0  0  1  3  5  6  4  2  1  0  0  0     0  0  1  3  8 12 10  5  3  1  0 -1  0     0  0  1  3  7  9  6  3  2  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  1  1  1  1  0  0  0 -1  0     0  0  0  1  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0 -1  0  0  0  0  0  0  0  0     0 -1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0    -1  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  7  7  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  4 10 13 11  7  3  1  0  0  0     0  0  0  1  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  4  8  9  7  4  2  1  0  0  0     0  0  0  1  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  1  3  3  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    637 REFLECTION INTENSITIES INTEGRATED


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
   281   0  1.091  5192389    0    125      75     2  0.04189  0.08030
   282   0  1.093  5192049    0    131      73     2  0.03871  0.06653
   283   0  1.100  5191344    0    130      83     4  0.03868  0.05663
   284   0  1.099  5192629    0    137      74     1  0.03838  0.08376
   285   0  1.102  5193589    0    126      70     2  0.04406  0.10108
   286   0  1.102  5192713    0    122      72     0  0.04426  0.11665
   287   0  1.104  5192358    0    118      75     0  0.04232  0.10883
   288   0  1.104  5192260    0    121      75     7  0.04556  0.07362
   289   0  1.102  5192418    0    119      74     3  0.04400  0.11281
   290   0  1.108  5193397    0    130      63     3  0.04643  0.10080
   291   0  1.107  5192066    0    134      77     3  0.04483  0.09409
   292   0  1.111  5191413    0    137      84     9  0.04203  0.06717
   293   0  1.110  5191618    0    141      80     7  0.03945  0.06207
   294   0  1.110  5192463    0    123      71     8  0.04392  0.06319
   295   0  1.114  5192299    0    128      76     4  0.04039  0.11152
   296   0  1.117  5192270    0    127      81     4  0.04496  0.10988
   297   0  1.115  5192406    0    130      76     8  0.04041  0.06437
   298   0  1.116  5193526    0    132      64     3  0.04112  0.07392
   299   0  1.116  5192673    0    137      74     2  0.04253  0.09050
   300   0  1.118  5191700    0    139      82     0  0.04264  0.08862

   451 OUT OF    451 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.54
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.13
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.408    32.091    34.488  88.746 108.393 111.883
 REC. CELL PARAMETERS   0.041672  0.033781  0.030739  83.760  70.651  67.319
 COORDINATES OF UNIT CELL A-AXIS   -14.380   -23.313    -0.977
 COORDINATES OF UNIT CELL B-AXIS    31.379    -5.458     3.922
 COORDINATES OF UNIT CELL C-AXIS    -1.146    12.149    32.257
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.024    -0.023     0.007
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999998 -0.002087  0.000549
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000515  0.005065  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.30   1244.13
 DETECTOR ORIGIN (PIXELS) AT                     1265.63   1237.48
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.28
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     451 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   48     30     52     62     55     29     53     62     53
 SIGMAB (degree)       0.036  0.039  0.040  0.037  0.034  0.036  0.037  0.036  0.035
 SIGMAR (degree)       0.031  0.027  0.030  0.032  0.030  0.029  0.032  0.031  0.029


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    451 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  1  0  0  0  0     0  0  0  0  0  0  1  1  1  0  0  0  0     0  0  0  0  1  2  3  3  2  1  1  0  0
    -1  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  1  2  2  2  1  1  0  0     0  0  0  1  4 10 14 14  9  5  2  1  0
     0  0  0  0  0  1  0  0  1  0  0  0  0     0  0  0  0  1  2  3  3  3  2  1  0  0     0  0  0  2  7 16 22 22 15  7  3  1  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  1  1  1  1  1  0  0  0     0  0  0  1  4  8 12 11  7  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  0  0  0  0  0  0  0  0  0  0     0  1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  3  8  9  8  4  2  1  0  0     0  0  0  1  5 11 11  9  5  2  1  0  0     0  0  0  1  4  7  7  5  2  1  0  0  0
     0  0  0  4 15 36 48 37 22 10  4  1  0     0  0  0  6 24 54 65 45 24 10  4  1  0     0  0  0  5 16 32 34 22 11  4  1  0  0
     0  0  1  7 28 59 74 61 35 16  6  2  0     0  0  2 14 49 91100 72 38 15  6  2  0     0  0  2 12 37 55 53 34 16  6  2  0  0
     0  0  1  6 17 33 41 32 16  8  3  1 -1     0  0  3 13 34 55 57 37 18  8  2  0  0     0  0  2 12 27 35 32 18  8  3  1  0  0
     0  0  0  1  4  7  7  4  3  1  0  0  0     0  0  0  4  9 13 11  6  3  1  0  0  0     0  0  0  4  7  9  7  3  1  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0 -1  0     0  0  0  0  1  1  1  1  0  0  0 -1  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  7  7  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0 -1     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  4  9 12 10  6  3  1  0  0  0     0  0  0  1  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  4  7  8  6  4  1  0  0  0  0     0  0  0  1  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  3  1  1  0  0  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    613 REFLECTION INTENSITIES INTEGRATED


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
   301   0  1.118  5190840    0    146      88     5  0.03811  0.05939
   302   0  1.125  5190603    0    145      93     2  0.03962  0.08145
   303   0  1.119  5191578    0    139      81     2  0.03806  0.06929
   304   0  1.127  5192208    0    140      76     2  0.03978  0.09061
   305   0  1.131  5191330    0    146      85     2  0.03990  0.07270
   306   0  1.130  5190882    0    153      88     2  0.04019  0.08438
   307   0  1.132  5191176    0    149      86     3  0.04123  0.07993
   308   0  1.129  5191206    0    156      84     3  0.04140  0.07323
   309   0  1.129  5191476    0    149      82     3  0.04084  0.07933
   310   0  1.135  5191299    0    144      86     1  0.04176  0.09508
   311   0  1.139  5191463    0    143      84     1  0.03974  0.09731
   312   0  1.141  5192305    0    127      77     3  0.04055  0.08929
   313   0  1.139  5192371    0    116      74     1  0.04313  0.10557
   314   0  1.144  5192579    0    125      70     0  0.04645  0.09140
   315   0  1.145  5191794    0    121      76     2  0.03858  0.06802
   316   0  1.148  5192106    0    119      75     3  0.03761  0.07027
   317   0  1.152  5192568    0    129      73     3  0.04060  0.08655
   318   0  1.151  5192635    0    136      76     3  0.04379  0.09641
   319   0  1.152  5191734    0    136      80     7  0.04082  0.06623
   320   0  1.151  5191756    0    137      85     5  0.04299  0.08806

   480 OUT OF    482 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.53
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.413    32.092    34.504  88.706 108.420 111.897
 REC. CELL PARAMETERS   0.041673  0.033782  0.030728  83.788  70.637  67.315
 COORDINATES OF UNIT CELL A-AXIS   -14.380   -23.317    -0.987
 COORDINATES OF UNIT CELL B-AXIS    31.380    -5.449     3.941
 COORDINATES OF UNIT CELL C-AXIS    -1.145    12.158    32.270
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.011    -0.021     0.004
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999998 -0.002087  0.000549
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000837  0.005181  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.36   1244.07
 DETECTOR ORIGIN (PIXELS) AT                     1265.25   1237.26
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.30
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     480 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   50     31     54     67     52     30     56     69     59
 SIGMAB (degree)       0.036  0.039  0.041  0.038  0.035  0.036  0.037  0.037  0.035
 SIGMAR (degree)       0.030  0.027  0.029  0.032  0.030  0.029  0.032  0.031  0.029


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    480 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  1  1  0  0  0     0  0  0  0  0  2  2  2  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0 -1     0  0  0  0  0  1  1  2  2  1  1  0  0     0  0  0  1  3  7 11 11  8  5  2  1  1
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  1  2  3  2  2  1  0  0     0  0  0  1  5 11 17 17 12  7  3  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  2  1  1  0  0  0     0  0  0  1  3  7 10  9  6  3  1  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0 -1  0 -1 -1  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0
     0  0  0  0  2  5  7  6  4  2  1  0  0     0  0  0  1  4  8 10  8  5  2  1  0  0     0  0  0  1  3  6  7  6  3  1  0  0  0
     0  0  0  3 12 29 38 32 20 10  4  1  1     0  0  1  7 24 51 58 44 23 10  4  1  0     0  0  1  6 18 33 35 25 13  6  2  0  0
     0  0  1  6 23 47 62 53 30 15  6  2  1     0  0  2 14 50 89100 73 36 16  6  2  0     0  0  2 13 38 60 59 40 19  8  3  1  0
     0  0  1  4 14 28 36 28 15  7  3  1  0     0  0  2 12 35 59 62 41 19  8  2  0  0     0  0  3 13 30 42 37 23 10  4  1  0 -1
     0  0  0  1  3  6  6  4  2  1  0  0  0     0  0  1  3  9 15 13  7  3  1  0  0  0     0  0  1  4  9 11  8  4  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  1  1  1  0  0  0  0 -1  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0 -1  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0  0 -1 -1  0  0  0


     0  0  0 -1  0  0 -1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  8  8  6  3  1  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  4 10 14 13  9  4  1  0  0  1     0  0  0  1  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  5  9 10  8  5  2  1  0  0  0     0  0  0  1  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  2  3  3  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    652 REFLECTION INTENSITIES INTEGRATED


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
   321   0  1.160  5191881    0    132      81     2  0.04293  0.08235
   322   0  1.160  5192759    0    129      65     4  0.04106  0.07102
   323   0  1.160  5192662    0    128      68     1  0.04365  0.08945
   324   0  1.164  5192605    0    131      72     1  0.04097  0.07947
   325   0  1.165  5192471    0    130      71     4  0.03881  0.06299
   326   0  1.168  5192310    0    122      71     2  0.03828  0.06859
   327   0  1.171  5192132    0    133      74     2  0.04201  0.08697
   328   0  1.171  5192340    0    140      71     0  0.04168  0.08800
   329   0  1.177  5191234    0    135      80     3  0.04046  0.06752
   330   0  1.175  5191968    0    124      74     3  0.04243  0.07382
   331   0  1.181  5192691    0    138      67     0  0.04344  0.11257
   332   0  1.180  5192040    0    133      76     4  0.04146  0.07410
   333   0  1.178  5190766    0    134      91     4  0.03843  0.06594
   334   0  1.177  5191831    0    141      80     0  0.03967  0.08832
   335   0  1.179  5192357    0    136      75     1  0.04061  0.08427
   336   0  1.182  5191922    0    139      76     3  0.03805  0.05996
   337   0  1.183  5191799    0    140      80     5  0.03883  0.06239
   338   0  1.185  5191355    0    153      86     5  0.04037  0.06812
   339   0  1.189  5191821    0    152      82     1  0.03993  0.09081
   340   0  1.190  5191744    0    133      80     1  0.04094  0.08088

   484 OUT OF    484 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.50
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.419    32.098    34.499  88.699 108.426 111.899
 REC. CELL PARAMETERS   0.041665  0.033776  0.030733  83.792  70.633  67.315
 COORDINATES OF UNIT CELL A-AXIS   -14.383   -23.323    -0.992
 COORDINATES OF UNIT CELL B-AXIS    31.385    -5.449     3.945
 COORDINATES OF UNIT CELL C-AXIS    -1.145    12.154    32.267
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.010    -0.031     0.004
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999998 -0.002087  0.000549
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000902  0.005181  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.42   1244.03
 DETECTOR ORIGIN (PIXELS) AT                     1265.23   1237.22
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.37
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     484 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   50     31     58     69     57     30     53     67     60
 SIGMAB (degree)       0.036  0.039  0.041  0.038  0.034  0.036  0.037  0.036  0.035
 SIGMAR (degree)       0.032  0.028  0.030  0.032  0.030  0.030  0.033  0.032  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    484 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  1  1  0  0  0     0  0  0  0  1  2  4  4  2  1  1  0  0
     0  0  0  0  0  0  0  0  1  1  0  0  0     0  0  0  0  0  1  3  4  3  2  1  0  0     0  0  0  1  5 11 16 16 12  6  3  1  0
     0  0  0  0  0  0  0  1  1  0  0  0  0     0  0  0  0  1  2  4  5  4  2  1  1  0     0  0  0  2  8 17 26 26 17  9  4  2  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  3  2  1  0  0  0     0  0  0  1  5 10 15 13  9  4  2  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0 -1 -1 -1  0  0  0  0     0  0  0  0  0  0 -1  0 -1  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0 -1 -1  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  1  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0
     0  0  0  1  3  6  8  6  4  2  1  0  0     0  0  0  1  5  9 11  8  5  2  1  1  0     0  0  0  1  3  6  6  5  3  1  0  0  0
     0  0  0  4 16 35 45 36 22 10  4  2  0     0  0  1  7 26 51 60 43 23 11  4  2  0     0  0  1  5 17 29 32 22 10  5  2  0  0
     0  0  1  7 29 60 78 64 35 17  7  3  0     0  0  3 16 53 93100 72 35 16  6  2  0     0  0  2 12 35 51 50 33 15  6  2  1  0
     0  0  1  5 19 38 44 32 17  8  4  1  0     0  0  3 14 38 63 60 37 18  8  3  1  0     0  0  3 14 30 39 31 17  7  3  1  0  0
     0  0  0  1  4  6  8  5  3  1  0  0  0     0  0  1  5 10 13 12  7  3  1  0  0  0     0  0  2  5  9 10  8  3  1  0  0 -1  0
     0  0  0  0  0  1  1  0  0  0 -1 -1  0     0  0  0  0  1  1  1  0  0  0  0 -1  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0 -1  0  0  0  0  0 -1 -1  0  0  0


     0  0  0  0 -1 -1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  4  6  6  4  2  1  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  4  8 10 10  6  2  1  0  0 -1     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  2  4  8  9  6  3  1  1  0  0  0     0  0  0  1  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  1  2  3  2  2  1  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0 -1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    638 REFLECTION INTENSITIES INTEGRATED


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
   341   0  1.190  5192907    0    129      67     4  0.04300  0.08687
   342   0  1.192  5193042    0    139      68     3  0.04080  0.08939
   343   0  1.195  5191559    0    140      78     0  0.04179  0.08483
   344   0  1.192  5191059    0    145      89     0  0.04091  0.08777
   345   0  1.198  5191539    0    148      85     3  0.04299  0.07528
   346   0  1.195  5191433    0    136      81     2  0.04114  0.07056
   347   0  1.201  5191444    0    143      84     3  0.04077  0.07533
   348   0  1.206  5191999    0    143      77     0  0.03930  0.08243
   349   0  1.203  5191279    0    131      82     2  0.04104  0.08330
   350   0  1.203  5191477    0    124      77     3  0.04241  0.07333
   351   0  1.203  5191986    0    118      76     2  0.04354  0.09194
   352   0  1.205  5192673    0    122      72     1  0.04426  0.11808
   353   0  1.205  5192398    0    115      75     2  0.04247  0.09601
   354   0  1.203  5191582    0    125      86     2  0.04336  0.10757
   355   0  1.204  5191652    0    139      86     1  0.04117  0.10118
   356   0  1.203  5191429    0    146      83     1  0.04218  0.08616
   357   0  1.206  5191397    0    139      82     1  0.04052  0.08555
   358   0  1.206  5192389    0    136      69     2  0.04147  0.07797
   359   0  1.206  5191363    0    128      80     1  0.04344  0.09002
   360   0  1.208  5191897    0    126      78     3  0.04122  0.08609

   473 OUT OF    475 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.54
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.411    32.089    34.495  88.736 108.415 111.895
 REC. CELL PARAMETERS   0.041676  0.033786  0.030737  83.757  70.630  67.308
 COORDINATES OF UNIT CELL A-AXIS   -14.380   -23.316    -0.982
 COORDINATES OF UNIT CELL B-AXIS    31.379    -5.450     3.927
 COORDINATES OF UNIT CELL C-AXIS    -1.147    12.160    32.260
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.016    -0.018     0.004
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999998 -0.002087  0.000549
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000645  0.005123  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.35   1244.00
 DETECTOR ORIGIN (PIXELS) AT                     1265.50   1237.27
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.30
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     473 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   50     32     55     69     58     32     50     67     52
 SIGMAB (degree)       0.036  0.040  0.041  0.039  0.035  0.036  0.037  0.037  0.035
 SIGMAR (degree)       0.032  0.027  0.029  0.032  0.031  0.033  0.034  0.032  0.029


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    473 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  1  0  0  0     0  0  0  0  1  2  4  3  3  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  3  2  2  1  0  0     0  0  0  1  4 10 14 14 11  6  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  4  3  2  1  0  0     0  0  0  2  6 14 20 21 16  8  4  2  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  2  1  1  0  0     0  0  0  1  4  8 12 12  8  4  2  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0 -1  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  1  3  7  9  7  4  2  1  0  0     0  0  0  2  5  9 11  9  5  2  1  0  0     0  0  0  1  3  5  7  6  3  1  0  0  0
     0  0  1  5 16 32 41 35 22 10  4  1  0     0  0  1  9 26 48 55 42 23 10  4  1  0     0  0  1  6 17 30 31 23 11  5  2  0  0
     0  0  1  8 27 54 70 56 34 16  7  3  0     0  0  3 16 50 89100 67 35 15  6  2  0     0  0  2 13 36 58 57 36 18  7  2  1  0
     0  0  1  7 19 35 40 32 17  8  3  1  0     0  0  3 16 38 59 59 38 18  8  3  1  0     0  0  3 14 30 40 37 20  9  4  1  0  0
     0  0  0  2  5  9 10  5  2  1  0  0  0     0  0  1  5 11 16 14  7  3  1  0 -1  0     0  0  1  4  8 11  8  4  2  0  0 -1  0
     0  0  0  1  1  1  1  1  0  0  0  0  0     0  0  0  1  2  2  1  1  0  0  0  0  0     0  0  0  1  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  1  2  1  1  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4  7  7  5  2  1  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  4 10 14 12  7  3  1  0  0  0     0  0  0  1  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  5 10 10  9  4  2  1  0  0  0     0  0  0  1  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  3  4  2  1  0  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    629 REFLECTION INTENSITIES INTEGRATED


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
   361   0  1.208  5192219    0    128      70     3  0.04480  0.06946
   362   0  1.210  5193027    0    122      63     2  0.04211  0.07158
   363   0  1.210  5192911    0    127      65     2  0.03884  0.08092
   364   0  1.213  5191916    0    128      78     4  0.03977  0.06745
   365   0  1.210  5190726    0    142      91     8  0.03958  0.06318
   366   0  1.209  5191314    0    142      89     3  0.04207  0.09907
   367   0  1.212  5191803    0    149      81     4  0.04290  0.09904
   368   0  1.211  5191020    0    152      92     5  0.04458  0.08207
   369   0  1.210  5190754    0    143      94     6  0.04084  0.07251
   370   0  1.211  5191454    0    139      85     2  0.04352  0.10173
   371   0  1.211  5192513    0    134      73     2  0.04577  0.10153
   372   0  1.216  5192057    0    122      75     1  0.04210  0.08180
   373   0  1.213  5192274    0    119      74     0  0.04511  0.10971
   374   0  1.214  5192641    0    125      70     1  0.04493  0.11602
   375   0  1.213  5190963    0    131      85     2  0.04256  0.09385
   376   0  1.213  5190354    0    132      90     6  0.04109  0.07171
   377   0  1.215  5190839    0    142      88     2  0.04543  0.09976
   378   0  1.216  5191411    0    133      85     1  0.04319  0.10830
   379   0  1.217  5191457    0    137      87     5  0.04347  0.09414
   380   0  1.217  5191429    0    134      87     4  0.04269  0.09579

   463 OUT OF    466 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.57
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.08
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.414    32.089    34.490  88.683 108.451 111.908
 REC. CELL PARAMETERS   0.041681  0.033786  0.030745  83.796  70.611  67.308
 COORDINATES OF UNIT CELL A-AXIS   -14.378   -23.319    -1.000
 COORDINATES OF UNIT CELL B-AXIS    31.376    -5.440     3.958
 COORDINATES OF UNIT CELL C-AXIS    -1.150    12.158    32.255
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.007    -0.033     0.001
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999998 -0.002087  0.000549
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.001351  0.005073  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.38   1244.06
 DETECTOR ORIGIN (PIXELS) AT                     1264.60   1237.40
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.28
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     463 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   49     27     52     68     56     31     54     64     56
 SIGMAB (degree)       0.036  0.039  0.042  0.039  0.034  0.036  0.037  0.036  0.035
 SIGMAR (degree)       0.033  0.029  0.030  0.032  0.031  0.033  0.034  0.032  0.031


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    463 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  1  1  0  0  0  0     0  0  0  0  0  1  1  1  2  1  0  0  0     0  0  0  0  1  3  4  4  3  2  0  0  0
     0  0  0  0  0  0  0  1  1  1  0  0  0     0  0  0  0  1  2  3  4  4  2  1  0  0     0  0  0  1  4 10 16 16 11  6  2  1  0
     0  0  0  0  0  0  1  1  1  0  0  0  0     0  0  0  0  1  2  4  5  4  2  1  0  0     0  0  0  1  7 15 24 25 17  9  4  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  3  2  1  0  0  0     0  0  0  1  4  8 13 13  9  4  2  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0 -1  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  3  7  8  6  4  1  0  0  0     0  0  0  1  4 10  9  7  4  2  1  0  0     0  0  0  1  3  7  7  5  3  1  1  0  0
     0  0  1  5 16 33 44 36 21  9  4  1  0     0  0  2  9 28 51 57 42 22 10  4  1  1     0  0  1  7 19 32 34 24 12  5  2  1  0
     0  0  1  7 29 57 74 60 35 16  7  3  1     0  0  3 17 53 94100 71 37 16  6  2  1     0  0  3 14 39 59 56 38 18  8  3  1  0
     0  0  1  5 17 31 41 32 19  8  3  1  0     0  0  4 15 39 58 59 37 19  8  3  1  0     0  1  4 14 32 40 35 20 10  4  1  0  0
     0  0  0  1  3  6  7  5  2  1  0  0  0     0  0  2  5 10 13 12  7  3  1  0  0  0     0  0  2  5  9 10  8  4  1  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0 -1  0     0  0  0  0  1  1  1  0  0  0  0 -1  0    -1  0  0  1  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0 -1  0  0  0  0  0  0  0  0  0  0     0  0 -1  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  3  3  3  1  1  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  8  9  7  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  4 10 14 12  8  4  2  0  0 -1     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  5  9 10  7  4  2  1  0  0 -1     0  0  0  1  1  1  1  1  0  0  0  0 -1     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  1  2  3  3  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
    -1  0  0  1  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    624 REFLECTION INTENSITIES INTEGRATED


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
   381   0  1.217  5191220    0    134      86     3  0.04412  0.09699
   382   0  1.221  5190529    0    142      94     2  0.04647  0.10144
   383   0  1.217  5190645    0    145      88     3  0.04636  0.09604
   384   0  1.216  5190257    0    150      93     5  0.04183  0.08274
   385   0  1.220  5190429    0    143      90     3  0.04327  0.08191
   386   0  1.221  5190896    0    137      92     6  0.04451  0.08966
   387   0  1.221  5191418    0    131      94     3  0.04213  0.12183
   388   0  1.221  5191735    0    135      82     4  0.03751  0.07496
   389   0  1.219  5191364    0    147      90     4  0.03718  0.06871
   390   0  1.224  5190895    0    144      92     6  0.04472  0.06826
   391   0  1.222  5191761    0    147      84     0  0.04057  0.10704
   392   0  1.224  5191256    0    145      84     2  0.04316  0.10746
   393   0  1.226  5190519    0    139     100     1  0.03909  0.09529
   394   0  1.225  5191836    0    133      82     2  0.04101  0.06892
   395   0  1.226  5192613    0    138      74     1  0.04026  0.09456
   396   0  1.226  5192198    0    131      78     0  0.03847  0.08767
   397   0  1.226  5191669    0    134      82     1  0.04158  0.08328
   398   0  1.227  5190862    0    129      92     3  0.03986  0.08713
   399   0  1.226  5191556    0    143      89     0  0.04474  0.10416
   400   0  1.230  5191909    0    143      75     0  0.04355  0.10343

   485 OUT OF    485 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.61
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.14
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.406    32.082    34.502  88.738 108.412 111.901
 REC. CELL PARAMETERS   0.041686  0.033795  0.030730  83.754  70.631  67.301
 COORDINATES OF UNIT CELL A-AXIS   -14.377   -23.311    -0.985
 COORDINATES OF UNIT CELL B-AXIS    31.372    -5.446     3.925
 COORDINATES OF UNIT CELL C-AXIS    -1.149    12.156    32.269
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.016    -0.023     0.002
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999998 -0.002087  0.000549
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000643  0.005081  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.35   1244.07
 DETECTOR ORIGIN (PIXELS) AT                     1265.50   1237.40
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.25
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     485 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   50     34     54     66     55     33     53     66     64
 SIGMAB (degree)       0.036  0.039  0.042  0.039  0.035  0.036  0.036  0.037  0.035
 SIGMAR (degree)       0.033  0.028  0.029  0.030  0.030  0.035  0.037  0.033  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    485 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  1  0  0  0  0     0  0  0  0  1  2  3  2  2  1  0  0  0
     0  0  0  0  0  0  0  1  1  0  0  0  0     0  0  0  0  1  1  2  2  2  1  0  0  0     0  0  0  1  3  7 10 10  7  4  1  0  0
     0  0  0  0  0  0  1  0  1  0  0  0  0     0  0  0  0  1  1  2  3  3  2  1  0  0     0  0  0  2  5 11 16 16 11  6  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  1  0  0  0     0  0  0  1  3  6  9  9  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0 -1  0  0  0 -1 -1  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  1  1  0  0  0 -1  0     0  0  0  0  0  0  0  1  1  0  0  0  0
     0  0  0  1  3  5  7  6  3  1  0  0  0     0  0  0  1  4  8 11  9  6  2  1  0  0     0  0  0  1  3  6  7  6  4  2  1  0  0
     0  0  0  4 14 28 34 30 18  9  3  1  0     0  0  1  6 23 48 57 47 26 12  5  2  0     0  0  1  6 19 37 40 30 15  7  2  1  0
     0  0  1  7 24 48 60 51 29 14  6  2  0     0  0  3 15 47 85100 76 40 18  8  2  1     0  0  3 15 41 67 71 48 23 10  4  1  0
     0  0  1  5 14 27 33 26 15  7  3  1  0     0  1  3 12 30 50 56 37 21  9  3  1 -1     0  1  4 14 29 43 43 24 12  5  2  0  0
     0  0  0  2  4  6  6  4  2  1  0  0  0     0  0  1  5  9 12 10  6  3  1  0  0  0     0  0  2  4  8 11  8  5  2  1  0  0  0
     0  0  0  0  1  1  1  0  0  0 -1 -1  0     0  0  0  1  2  2  1  0  0  0  0 -1  0     0  0  0  1  1  1  1  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0 -1 -1  0  0  0


     0 -1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  2  1  1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 11 10  7  4  2  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  7 16 22 20 12  5  2  1  0  0     0  0  0  1  2  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  2  8 15 17 14  7  3  1  0  0  0     0  0  0  1  2  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  2  4  5  3  1  1  0  0  0  0     0  0  0  1  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
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
   401   0  1.228  5190922    0    144      86     1  0.03980  0.07056
   402   0  1.227  5191031    0    139      89     2  0.04054  0.07910
   403   0  1.230  5191758    0    136      86     0  0.03965  0.09663
   404   0  1.231  5192266    0    122      83     2  0.04301  0.10392
   405   0  1.231  5192326    0    119      75     8  0.04437  0.08250
   406   0  1.231  5192061    0    123      79     3  0.04653  0.11043
   407   0  1.230  5191895    0    118      79     3  0.04258  0.09103
   408   0  1.232  5191772    0    138      82     3  0.03994  0.07856
   409   0  1.233  5190140    0    145     100     0  0.03917  0.09255
   410   0  1.232  5190672    0    148      96     3  0.03932  0.08154
   411   0  1.234  5191318    0    154      90     2  0.04079  0.08372
   412   0  1.234  5190542    0    151      96     3  0.04045  0.07840
   413   0  1.235  5191679    0    142      80     3  0.04203  0.06942
   414   0  1.237  5191957    0    137      79     4  0.04091  0.07597
   415   0  1.240  5193178    0    120      66     2  0.04074  0.09903
   416   0  1.242  5192887    0    125      66     1  0.04271  0.10497
   417   0  1.238  5192117    0    128      73     5  0.04161  0.07433
   418   0  1.237  5191826    0    129      77     2  0.04080  0.07670
   419   0  1.239  5191546    0    130      80     5  0.04167  0.07162
   420   0  1.238  5192059    0    130      75     4  0.04183  0.07702

   467 OUT OF    471 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.63
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.407    32.085    34.502  88.748 108.402 111.896
 REC. CELL PARAMETERS   0.041681  0.033792  0.030728  83.749  70.639  67.304
 COORDINATES OF UNIT CELL A-AXIS   -14.376   -23.313    -0.979
 COORDINATES OF UNIT CELL B-AXIS    31.375    -5.446     3.916
 COORDINATES OF UNIT CELL C-AXIS    -1.146    12.156    32.269
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.024    -0.019    -0.002
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999998 -0.002087  0.000549
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000504  0.005134  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.31   1244.06
 DETECTOR ORIGIN (PIXELS) AT                     1265.65   1237.32
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.29
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     467 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   46     29     56     64     54     32     60     66     45
 SIGMAB (degree)       0.036  0.040  0.042  0.039  0.035  0.036  0.036  0.036  0.035
 SIGMAR (degree)       0.033  0.027  0.028  0.031  0.030  0.036  0.037  0.034  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    467 STRONG REFLECTIONS *****


     0 -1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  1  1  0  0  0     0  0  0  0  1  2  3  3  2  1  0  0  0
     0  0  0  0  1  0  0  0  1  1  0  0  0     0  0  0  0  1  1  2  3  3  2  1  0  0     0  0  0  1  4  8 12 12  9  5  2  0  0
     0  0  0  0  1  0  0  1  1  0  0  0  0     0  0  0  0  1  1  3  4  4  2  1  0  0     0  0  0  2  6 12 19 20 14  7  3  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  2  2  1  0  0  0     0  0  0  1  3  7 10 11  7  4  2  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0 -1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0
     0  0  0  1  3  7  7  6  3  2  0  0  0     0  0  0  2  6 11 12 10  5  3  1  0  0     0  0  0  2  4  8  8  7  4  2  1  0  0
     0  0  1  4 14 30 41 33 20 10  4  1  0     0  0  1  8 27 54 64 47 26 12  5  2  0     0  0  1  7 22 38 41 28 15  7  3  1  0
     0  0  1  8 27 49 65 57 31 15  7  2  1     0  0  4 18 53 89100 76 38 18  8  3  2     0  0  4 16 46 69 67 46 22 10  4  1  0
     0  0  1  6 18 30 35 29 16  8  3  1  0     0  0  4 15 38 55 56 38 19 10  4  1  0     0  0  4 14 33 47 40 24 11  5  2  0  0
     0  0  1  2  5  8  7  4  3  1  0  0  0     0  0  2  5 10 13 10  6  3  1  0  0  0     0  0  1  4  8 11  8  4  2  0  0  0  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  1  1  1  1  0  0  0  0 -1  0     0  0  0  0  1  1  1  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0 -1  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  2  1  1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  3  7 11 10  6  3  2  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  2  7 18 22 18 11  5  2  1  0  0     0  0  0  1  2  3  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  1
     0  0  2  7 16 19 13  7  3  1  0  0  0     0  0  0  1  3  3  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  2  4  4  3  1  1  0  0  0  0     0  0  0  1  1  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    618 REFLECTION INTENSITIES INTEGRATED


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
   421   0  1.240  5191639    0    135      80     4  0.04078  0.07682
   422   0  1.242  5191938    0    141      78     6  0.04039  0.08116
   423   0  1.244  5192046    0    144      75     2  0.04269  0.08973
   424   0  1.247  5191077    0    140      88     3  0.03971  0.09189
   425   0  1.244  5191258    0    140      83     3  0.04236  0.07198
   426   0  1.245  5190713    0    135      90     2  0.03864  0.07646
   427   0  1.244  5191657    0    139      79     0  0.04176  0.09030
   428   0  1.243  5191535    0    139      80     4  0.04015  0.08739
   429   0  1.243  5191462    0    143      76     2  0.04143  0.09032
   430   0  1.248  5191210    0    147      86     6  0.04024  0.06886
   431   0  1.248  5191207    0    137      83     5  0.04091  0.06987
   432   0  1.245  5191444    0    145      83     4  0.04202  0.08142
   433   0  1.249  5190702    0    150      91     1  0.04139  0.08629
   434   0  1.246  5190643    0    154      96     1  0.04144  0.07859
   435   0  1.246  5191183    0    144      85     1  0.03950  0.07645
   436   0  1.247  5191657    0    139      82     2  0.04215  0.08510
   437   0  1.248  5191499    0    129      82     1  0.04484  0.09099
   438   0  1.249  5192034    0    123      75     5  0.04206  0.07793
   439   0  1.252  5192285    0    125      72     5  0.04268  0.08552
   440   0  1.250  5191985    0    129      77     2  0.04408  0.08997

   489 OUT OF    493 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.51
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.13
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.414    32.089    34.490  88.720 108.430 111.899
 REC. CELL PARAMETERS   0.041676  0.033787  0.030743  83.767  70.620  67.308
 COORDINATES OF UNIT CELL A-AXIS   -14.379   -23.320    -0.987
 COORDINATES OF UNIT CELL B-AXIS    31.378    -5.445     3.932
 COORDINATES OF UNIT CELL C-AXIS    -1.145    12.158    32.255
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.015    -0.024    -0.001
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999998 -0.002087  0.000549
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000792  0.005059  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.37   1244.08
 DETECTOR ORIGIN (PIXELS) AT                     1265.33   1237.43
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.32
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     489 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     32     56     66     55     29     58     72     66
 SIGMAB (degree)       0.036  0.040  0.042  0.039  0.035  0.035  0.036  0.036  0.035
 SIGMAR (degree)       0.032  0.028  0.029  0.031  0.031  0.035  0.037  0.034  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    489 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  1  1  0  0  0  0     0  0  0  0  0  1  1  2  1  1  0  0  0     0  0  0  0  1  3  4  4  2  1  0  0  0
     0  0  0  0  0  0  1  1  2  1  0  0  0     0  0  0  0  1  2  4  4  4  2  1  0  0     0  0  0  2  5 12 18 17 11  5  2  1  0
     0  0  0  0  0  0  1  1  2  1  0  0  0     0  0  0  0  1  2  5  7  6  3  1  0  0     0  0  0  2  8 18 28 28 19  9  4  1  1
     0  0  0  0  0  0  0  1  1  1  0  0  0     0  0  0  0  1  1  3  3  3  2  1  0  0     0  0  0  2  6 10 15 15  9  5  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  1  0  0  0  0     0  0  0  0  1  2  3  2  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0 -1 -1  0     0  0  0  0  0  1  1  1  0  0  0 -1  0     0  0  0  0  0  0  1  1  0  0  0  0  0
     0  0  0  1  4  8  9  6  3  1  0  0  0     0  0  0  2  6 11 11  8  4  2  1  0  0     0  0  0  1  3  6  6  5  3  2  1  0  0
     0  0  1  6 19 40 51 39 22 10  4  1  0     0  0  1 10 29 56 64 45 22 10  4  1  1     0  0  1  6 18 32 34 24 12  5  2  1  1
     0  0  2 10 31 63 80 65 38 18  7  2  1     0  0  3 17 50 92100 71 37 18  7  2  1     0  0  3 13 34 55 52 35 17  8  3  1  1
     0  0  2  8 23 36 43 34 20 10  4  1  0     0  1  4 16 39 55 55 37 19  9  4  1  0     0  0  3 11 27 34 30 18  9  4  1  0  0
     0  0  1  3  6  8  9  6  4  2  0  0  0     0  0  2  6 11 13 11  6  4  2  0  0  0     0  0  1  4  7  8  6  3  2  1  0  0  0
     0  0  0  1  1  1  1  0  0  0  0 -1  0     0  0  0  1  1  1  1  1  0  0  0 -1  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0 -1 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0 -1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4  7  7  5  3  1  0  0  0     0  0  0  0  0  0  1  0  0  0 -1  0  1     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  4  9 12 10  7  4  2  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  1     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  4  9  9  6  4  2  1  0  0  0     0  0  0  1  1  1  1  1  1  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0
     0  0  0  1  2  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    656 REFLECTION INTENSITIES INTEGRATED


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
   441   0  1.252  5191628    0    129      81     1  0.04075  0.09151
   442   0  1.254  5191541    0    138      80     2  0.04260  0.08101
   443   0  1.255  5192189    0    146      77     3  0.04152  0.08955
   444   0  1.253  5191497    0    144      83     6  0.04092  0.06884
   445   0  1.252  5190484    0    139      93     8  0.04202  0.07550
   446   0  1.251  5191277    0    138      87     4  0.04311  0.09575
   447   0  1.252  5192102    0    123      76     3  0.04073  0.08613
   448   0  1.255  5192301    0    120      75     3  0.03866  0.08934
   449   0  1.256  5192971    0    127      64     2  0.03938  0.08128
   450   0  1.254  5192036    0    136      74     1  0.04388  0.09599
   451   0  1.257  5191409    0    134      82     1  0.04385  0.09199
   452   0  1.257  5190778    0    137      86     4  0.04079  0.07835
   453   0  1.259  5191600    0    133      75     3  0.04080  0.08373
   454   0  1.259  5191830    0    126      80     4  0.03928  0.08434
   455   0  1.257  5191549    0    125      82     2  0.04265  0.11958
   456   0  1.262  5191887    0    136      80     6  0.03791  0.08620
   457   0  1.260  5191199    0    137      83     8  0.03835  0.06453
   458   0  1.259  5190579    0    148      90     3  0.04364  0.08940
   459   0  1.260  5189801    0    156     102     8  0.04270  0.08555
   460   0  1.259  5190154    0    154     100     6  0.04099  0.09506

   476 OUT OF    480 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.53
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.410    32.085    34.511  88.728 108.424 111.903
 REC. CELL PARAMETERS   0.041683  0.033792  0.030724  83.758  70.622  67.301
 COORDINATES OF UNIT CELL A-AXIS   -14.377   -23.316    -0.985
 COORDINATES OF UNIT CELL B-AXIS    31.376    -5.442     3.925
 COORDINATES OF UNIT CELL C-AXIS    -1.144    12.164    32.276
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.020    -0.017    -0.003
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999998 -0.002087  0.000549
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000689  0.005086  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.31   1244.03
 DETECTOR ORIGIN (PIXELS) AT                     1265.40   1237.35
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.30
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     476 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   51     33     61     70     60     30     48     63     52
 SIGMAB (degree)       0.036  0.039  0.042  0.039  0.035  0.035  0.036  0.036  0.035
 SIGMAR (degree)       0.032  0.028  0.029  0.030  0.030  0.035  0.038  0.034  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    476 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  1  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0
     0  0  0  0  0  0  0  1  1  0  0  0  0     0  0  0  0  0  1  2  2  2  1  0  0  0     0  0  0  1  3  7 11 11  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  3  3  2  1  0  0     0  0  0  2  5 11 17 17 11  6  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  2  1  1  0  0  0     0  0  0  1  3  6  9  9  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0 -1  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  1  3  6  7  5  2  1  0  0  0     0  0  0  2  6 10 11  8  4  2  1  0  0     0  0  0  1  4  7  9  6  3  2  1  0  0
     0  0  1  5 16 30 37 30 16  7  3  1  0     0  0  1  9 29 54 60 43 23 11  4  1  0     0  0  1  7 22 40 43 30 17  8  3  1  0
     0  0  1  8 23 48 60 50 29 14  5  2  1     0  0  4 17 48 89100 73 39 19  7  3  1     0  0  3 15 42 71 72 48 24 12  4  1  0
     0  0  1  6 15 28 33 27 15  8  3  1  0     0  0  3 14 32 54 56 39 20 10  3  1  1     0  0  3 13 31 45 41 25 12  6  2  1  0
     0  0  0  1  4  7  6  5  3  1  0  0  0     0  0  1  4 10 13 10  7  3  1  0  0  0     0  0  1  4  8 11  8  5  2  1  0  0  0
     0  0  0  0  0  1  1  0  0  0  0 -1  0     0  0  0  1  1  2  1  1  0  0  0 -1  0     0  0  0  1  1  1  1  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0 -1 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  3  2  2  1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  2  7 11 12  9  6  2  1  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  2  6 15 22 20 13  7  3  1  0  0     0  0  0  1  2  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  5 13 16 12  7  3  1  0  0  0     0  0  0  1  2  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  4  2  1  1  0  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    639 REFLECTION INTENSITIES INTEGRATED


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
   461   0  1.264  5190758    0    148      94     8  0.04673  0.09863
   462   0  1.262  5190702    0    141      97     9  0.04236  0.10025
   463   0  1.260  5190627    0    138      93     7  0.04194  0.08690
   464   0  1.260  5189951    0    134     100     7  0.04121  0.08426
   465   0  1.262  5190450    0    150     100     7  0.04261  0.09607
   466   0  1.264  5190057    0    150     104     9  0.04284  0.09722
   467   0  1.262  5189699    0    142     100    10  0.04116  0.06796
   468   0  1.263  5189880    0    144     101     7  0.04206  0.08603
   469   0  1.263  5190919    0    137      86     5  0.04356  0.09585
   470   0  1.266  5190567    0    134      93     7  0.04268  0.11754
   471   0  1.262  5190184    0    138     100     4  0.04240  0.12468
   472   0  1.260  5189953    0    148      99     8  0.04330  0.10644
   473   0  1.262  5190005    0    139      95    10  0.04271  0.08515
   474   0  1.264  5190355    0    136      99     7  0.04351  0.11174
   475   0  1.261  5191053    0    140      91     7  0.04192  0.09633
   476   0  1.261  5190515    0    128      94     8  0.04339  0.09221
   477   0  1.263  5191205    0    130      86     6  0.04233  0.10920
   478   0  1.262  5190821    0    132      93     5  0.04289  0.12579
   479   0  1.263  5190030    0    130      95     4  0.04295  0.10860
   480   0  1.265  5190512    0    143      89     3  0.04347  0.10878

   500 OUT OF    517 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.63
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.17
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.409    32.084    34.488  88.720 108.429 111.902
 REC. CELL PARAMETERS   0.041685  0.033792  0.030745  83.766  70.621  67.304
 COORDINATES OF UNIT CELL A-AXIS   -14.376   -23.315    -0.982
 COORDINATES OF UNIT CELL B-AXIS    31.374    -5.441     3.932
 COORDINATES OF UNIT CELL C-AXIS    -1.145    12.164    32.251
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.011    -0.013    -0.003
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999998 -0.002087  0.000549
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000869  0.005101  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.32   1244.07
 DETECTOR ORIGIN (PIXELS) AT                     1265.18   1237.37
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.26
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     500 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   50     38     58     65     54     32     56     68     63
 SIGMAB (degree)       0.036  0.040  0.042  0.039  0.036  0.035  0.036  0.036  0.036
 SIGMAR (degree)       0.034  0.029  0.030  0.032  0.031  0.036  0.039  0.034  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    500 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  0  1  1  1  0  0  0  0     0  0  0  0  1  3  4  3  2  1  0  0  0
     0  0  0  0  0  0  1  1  1  1  0  0  0     0  0  0  0  1  2  4  4  3  2  1  0  0     0  0  0  1  6 11 17 16  9  5  2  0  0
     0  0  0  0  0  1  1  1  1  1  0  0  0     0  0  0  0  1  3  5  6  5  3  1  1  0     0  0  0  2  8 17 26 26 17  8  3  1  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  1  1  1  3  3  2  1  0  0  0     0  0  0  2  5 10 13 14  9  5  2  1  0
     0  1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0 -1  0     0  0  0  1  1  2  3  3  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  4  8  9  7  3  1  0  0  0     0  0  0  2  6 11 11  9  5  2  1  0  0     0  0  0  1  4  7  8  7  4  2  1  0  0
     0  0  1  5 19 36 46 37 21  9  4  1  1     0  0  1  9 32 56 62 44 23 11  4  1  1     0  0  1  7 23 37 39 26 15  7  3  1  1
     0  0  1  9 30 57 73 62 36 18  7  2  1     0  0  3 17 53 92100 72 38 18  7  3  1     0  0  2 13 39 63 61 41 21 10  4  1  0
     0  0  2  8 20 36 40 32 19 10  4  1  0     0  0  4 15 37 58 54 36 20  9  4  1  0     0  0  3 11 28 38 33 20 10  4  2  1  0
     0  0  1  3  7  9  9  6  3  2  1  0  0     0  0  2  6 12 14 12  6  3  1  0  0  0     0  0  1  4  8  8  7  4  2  0  0  0  0
     0  0  0  0  1  1  1  0  0  0  0 -1  0     0  0  0  1  2  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1 -1  0     0  0  0  0  0  0  0  0  0 -1 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  3  3  2  1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 10 11  8  5  2  1  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  5 13 18 16 10  5  2  1  0  0     0  0  0  1  2  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  5 11 12  9  5  3  1  0  0 -1     0  0  0  1  2  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  3  2  1  0  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    659 REFLECTION INTENSITIES INTEGRATED


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
   481   0  1.264  5190881    0    138      87     3  0.04640  0.11450
   482   0  1.264  5190917    0    131      88     4  0.04206  0.13591
   483   0  1.264  5191634    0    120      80     6  0.03984  0.09410
   484   0  1.264  5192156    0    113      74     4  0.04193  0.09356
   485   0  1.263  5191856    0    117      80     4  0.04021  0.11122
   486   0  1.265  5191826    0    122      75     5  0.04286  0.09730
   487   0  1.263  5191407    0    115      81     9  0.03995  0.06517
   488   0  1.263  5191224    0    125      84     5  0.04226  0.08823
   489   0  1.268  5191199    0    135      83     3  0.04007  0.08459
   490   0  1.266  5190235    0    144      93     4  0.04161  0.07573
   491   0  1.270  5190056    0    146      96     6  0.03980  0.07321
   492   0  1.269  5190200    0    150      97     7  0.03997  0.07565
   493   0  1.268  5190356    0    151      94     6  0.03889  0.08274
   494   0  1.270  5189803    0    141     106     6  0.04114  0.09208
   495   0  1.270  5190701    0    137      94     3  0.04395  0.10177
   496   0  1.267  5191200    0    132      84     4  0.04035  0.09587
   497   0  1.265  5190668    0    135      87     4  0.04202  0.09654
   498   0  1.265  5190049    0    151      97     4  0.03922  0.10356
   499   0  1.269  5189224    0    149     109    11  0.04427  0.08411
   500   0  1.270  5189670    0    155      97     6  0.04362  0.09151

   470 OUT OF    480 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.59
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.15
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.417    32.092    34.495  88.718 108.436 111.898
 REC. CELL PARAMETERS   0.041673  0.033783  0.030740  83.767  70.614  67.309
 COORDINATES OF UNIT CELL A-AXIS   -14.379   -23.322    -0.992
 COORDINATES OF UNIT CELL B-AXIS    31.381    -5.444     3.937
 COORDINATES OF UNIT CELL C-AXIS    -1.148    12.160    32.260
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.009    -0.028    -0.001
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.039
 LAB COORDINATES OF ROTATION AXIS -0.999998 -0.002087  0.000549
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000825  0.005065  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.37   1244.07
 DETECTOR ORIGIN (PIXELS) AT                     1265.29   1237.41
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.36
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     470 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   52     29     53     67     57     30     57     63     49
 SIGMAB (degree)       0.036  0.039  0.041  0.039  0.035  0.035  0.035  0.035  0.035
 SIGMAR (degree)       0.034  0.030  0.030  0.030  0.030  0.036  0.039  0.035  0.031


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    470 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0
     0  0  0  0  0  0  1  1  1  0  0  0  0     0  0  0  0  0  1  2  2  1  1  0  0  0     0  0  0  0  1  3  5  4  2  1  0  0  0
     0  0  0  0  0  0  1  2  2  1  0  0  0     0  0  0  0  1  2  5  7  4  2  1  0  0     0  0  0  2  8 15 21 19 11  5  2  0  0
     0  0  0  0  0  0  1  2  2  1  0  0  0     0  0  0  0  2  4  7  8  5  3  1  0  0     0  0  0  3 10 22 31 29 18  8  3  1  0
     0  0  0  0  0  0  1  1  1  0  0  0  0     0  0  0  0  1  2  4  4  3  1  0  0  0     0  0  0  2  6 12 16 14  9  4  2  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  1  2  3  3  3  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0 -1  0  0 -1  0  0  0  0     0  0  0  0  0 -1  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0 -1  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0 -1 -1  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  2  5  9  9  7  3  1  0  0  0     0  0  0  2  6 11 11  8  4  2  0  0  0     0  0  0  1  4  6  6  5  3  1  0  0  0
     0  0  1  7 24 45 52 38 20  9  3  1  0     0  0  1 10 32 60 66 44 22 10  4  1  1     0  0  1  6 20 34 37 25 12  6  2  1  0
     0  0  2 10 36 71 81 61 34 16  6  2  1     0  0  3 16 53 96100 68 36 17  7  2  1     0  0  2 11 35 58 55 36 19  8  3  1  0
     0  0  1  8 22 40 39 30 18  8  3  1  0     0  0  3 13 33 55 49 32 19  8  3  2  1     0  0  2  9 23 36 28 17 10  5  1  1  0
     0  0  1  3  7  9  8  6  3  2  0  0  0     0  0  1  4  9 11 10  7  4  2  0  0  0     0  0  1  2  5  7  6  4  2  1  0  0  0
     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  1  1  1  1  0  0  0  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0 -1 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0 -1 -1  0 -1 -1 -1  0  0  0  0  0     0  0 -1 -1  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  7  7  6  3  1  1  0  0     0  0  0  0  0  1  0  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  4 10 13 11  7  4  2  0  0  0     0  0  0  1  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  4  8  9  6  4  2  1  0  0  0     0  0  0  1  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  2  1  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    638 REFLECTION INTENSITIES INTEGRATED


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
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  4  0  0  0  0  0     0  0  0  0  5 11 14 12  6  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  2  3  4  4  4  0  0  0  0     0  0  0  2  8 17 22 19 10  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  0  4 10 13 10  5  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  7  8  6  3  0  0  0  0     0  0  0  2  6 11 12  8  4  0  0  0  0     0  0  0  0  4  7  8  5  3  0  0  0  0
     0  0  0  4 17 37 44 32 15  6  0  0  0     0  0  0  7 28 57 62 42 18  7  0  0  0     0  0  0  5 20 37 37 24 10  4  0  0  0
     0  0  0  7 28 59 71 53 25 10  0  0  0     0  0  2 12 47 91100 69 31 11  4  0  0     0  0  2  9 34 59 60 38 16  6  0  0  0
     0  0  0  4 16 36 43 29 13  5  0  0  0     0  0  0  7 29 58 63 39 17  6  0  0  0     0  0  0  6 22 38 38 22  9  3  0  0  0
     0  0  0  0  4  8  9  5  3  0  0  0  0     0  0  0  2  7 14 14  7  3  0  0  0  0     0  0  0  2  5 10  9  4  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 10  9  6  3  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3 10 16 15  9  4  2  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  7 11 10  5  3  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  2 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  4  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  6  7  6  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  2  0  0  0  0  0  0     0  0  0  3  5  9 11 10  6  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  3  4  6  7  5  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  3  3  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  3  3  0  0  0  0     0  0  0  0  0  0  4  3  4  0  0  0  0
     0  0  0  0  3  5  7  7  5  3  0  0  0     0  0  0  0  5 12 17 15  9  5  0  0  0     0  0  0  0  5 11 16 13  8  6  0  0  0
     0  0  0  4 13 27 33 26 15  7  0  0  0     0  0  0  7 26 56 65 47 25 11  6  0  0     0  0  0  5 20 41 46 32 17  8  6  0  0
     0  0  2  6 24 47 54 40 21  9  0  0  0     0  0  3 14 52 97100 66 33 14  6  0  0     0  0  3 12 38 63 60 38 18  7  0  0  0
     0  0  2  7 17 32 34 23 11  5  0  0  0     0  0  5 16 43 68 63 38 17  7  0  0  0     0  0  3 13 32 42 35 19  9  3  0  0  0
     0  0  3  6  8 10  9  5  3  0  0  0  0     0  0  3  9 16 21 16  9  4  0  0  0  0     0  0  3  6 11 12  9  4  3  0  0  0  0
     0  0  0  0  4  3  0  0  0  0  0  0  0     0  0  0  4  4  3  3  0  0  0  0  0  0     0  0  0  0  2  2  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  4  5  6  6  6  0  0  0  0     0  0  0  0  0  0  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4  5 10 11  9  7  5  0  0  0     0  0  0  0  3  3  3  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  5 10 14 12  8  5  4  0  0  0     0  0  0  0  4  5  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  6  8  9  7  4  3  0  0  0  0     0  0  0  0  4  4  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  6  6  5  4  3  0  0  0  0  0     0  0  0  0  0  3  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  3 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  3  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  4  0  0  0  0  0     0  0  0  0  3  7 11 12  9  7  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  4  6  0  0  0  0     0  0  0  4  6 13 20 21 16 10  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  4  5 10 14 13  9  6  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  3  3  3  4  3  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  4  6  6  5  3  0  0  0     0  0  0  0  3  8 10  9  6  5  4  0  0     0  0  0  0  3  6  6  6  5  4  0  0  0
     0  0  0  4 11 25 34 32 22 13  8  0  0     0  0  0  6 22 43 50 41 25 14  8  0  0     0  0  0  6 19 30 30 22 13  7  5  0  0
     0  0  3  9 26 51 68 60 40 22 12  0  0     0  0  6 20 54 91100 76 45 23 10  0  0     0  0  5 19 45 64 59 40 21 10  6  0  0
     0  0  4 11 25 43 50 41 25 13  8  0  0     0  3  8 25 55 79 78 53 29 14  6  0  0     0  0  8 24 46 57 48 28 14  7  4  0  0
     0  0  3  7 10 13 13 10  6  3  0  0  0     0  0  5 11 21 26 22 14  7  3  0  0  0     0  0  6 11 18 19 14  8  4  0  0  0  0
     0  0  0  3  3  3  2  0  0  0  0  0  0     0  0  0  4  4  4  3  0  0  0  0  0  0     0  0  0  0  3  3  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  4  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4  6  8  7  5  4  2  0  0  0     0  0  0  0  0  3  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  8 15 18 14  9  5  3  0  0  0     0  0  0  0  6  4  4  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0 10 16 18 12  7  3  2  0  0  0     0  0  0  0  7  5  4  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  8  7  6  4  3  0  0  0  0  0     0  0  0  0  0  4  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  4 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  2  3  4  3  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  4  2  0  0  0  0     0  0  0  2  5 12 20 18 11  6  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  5  6  5  0  0  0  0     0  0  0  3  6 17 29 29 19  9  5  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  3  2  0  0  0  0     0  0  0  3  4  9 15 15 10  6  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4  9 11  8  4  2  0  0  0     0  0  0  3  7 13 13  9  4  2  0  0  0     0  0  0  3  6  8  7  4  2  0  0  0  0
     0  0  0  5 19 42 56 44 24 11  5  0  0     0  0  3 10 33 62 68 46 23 10  5  0  0     0  0  3  9 26 40 37 23 11  4  3  0  0
     0  0  3  8 26 59 82 69 40 19  9  0  0     0  0  4 15 47 88100 72 38 17  7  0  0     0  0  4 14 36 56 54 35 17  7  4  0  0
     0  0  2  6 14 31 42 36 21 11  6  0  0     0  0  4 10 27 47 53 37 20  9  5  0  0     0  0  5 10 21 31 29 18  9  4  0  0  0
     0  0  0  2  4  6  7  6  4  3  0  0  0     0  0  0  4  6  9  9  7  4  3  0  0  0     0  0  0  0  5  7  6  3  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  5  9 11  9  5  3  0  0  0  0     0  0  0  0  0  3  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  6 12 16 13  7  3  2  0  0  0     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  5  7  9  7  4  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  5  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  5 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  4  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  0  3  9 15 12  6  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  4  4  3  0  0  0  0     0  0  0  0  4 13 22 20 10  5  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  3  0  0  0  0  0     0  0  0  0  3  6 11 10  5  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  8 10  7  3  0  0  0  0     0  0  0  0  5 10 12  8  3  0  0  0  0     0  0  0  0  3  6  5  3  2  0  0  0  0
     0  0  0  3 12 34 47 34 16  6  0  0  0     0  0  0  4 22 53 61 39 17  7  0  0  0     0  0  0  4 16 32 31 18  8  3  0  0  0
     0  0  0  3 17 52 73 54 27 11  0  0  0     0  0  0  7 36 86100 64 29 12  0  0  0     0  0  0  6 26 53 53 31 14  5  0  0  0
     0  0  0  2  8 25 36 26 13  5  0  0  0     0  0  0  3 18 45 51 31 14  6  0  0  0     0  0  0  3 14 29 29 16  7  3  0  0  0
     0  0  0  0  0  4  5  4  3  0  0  0  0     0  0  0  0  3  7  8  5  3  0  0  0  0     0  0  0  0  4  5  5  3  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  5  7  6  4  3  0  0  0  0     0  0  0  0  0  2  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4  7 12 10  6  3  2  0  0  0     0  0  0  0  0  3  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  4  6  5  3  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  4  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  6 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  2  6  8  8  5  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  5  7  9  0  0  0  0  0     0  0  0  0  4  8 12 11  8  6  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  4  5  5  0  0  0  0     0  0  0  0  4  9 14 12  7  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  3  0  0  0  0  0     0  0  0  0  3  5  6  5  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  8 11  9  5  0  0  0  0     0  0  0  0  4 11 13  9  4  2  0  0  0     0  0  0  0  3  5  6  4  3  0  0  0  0
     0  0  0  2 10 29 39 30 15  6  0  0  0     0  0  0  4 23 56 63 41 18  7  0  0  0     0  0  0  4 17 33 33 20  8  3  0  0  0
     0  0  0  4 18 44 54 39 19  7  0  0  0     0  0  3 11 47 97100 62 28 10  0  0  0     0  0  0 11 39 67 61 35 15  5  0  0  0
     0  0  0  4 10 23 27 18  8  4  0  0  0     0  0  2  8 32 58 56 31 14  5  0  0  0     0  0  2  9 31 47 39 20  8  4  0  0  0
     0  0  0  0  4  4  4  3  0  0  0  0  0     0  0  0  3  7 12 11  6  2  0  0  0  0     0  0  0  3  8 11  9  5  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  0  0  0  0  0  0  0     0  0  0  0  0  2  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  2  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  5  7  6  4  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  5 11 15 12  7  4  3  0  0  0     0  0  0  0  3  3  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  6 12 13  9  5  3  0  0  0  0     0  0  0  0  4  4  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  6  6  4  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  7 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  5  8  6  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  5 10 10  0  0  0  0     0  0  0  0  3  5 11 14 12  8  0  0  0
     0  0  0  0  0  0  3  0  0  0  0  0  0     0  0  0  0  0  3  5  7 10  0  0  0  0     0  0  0  3  6 12 19 21 15  8  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  4  0  0  0  0  0     0  0  0  3  5  8 11 11  7  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  4  6  4  3  0  0  0     0  0  0  0  2  4  6  5  3  2  0  0  0     0  0  0  0  2  3  4  3  2  0  0  0  0
     0  0  0  2  7 19 32 32 20  9  5  0  0     0  0  0  3 13 34 49 39 21  9  3  0  0     0  0  0  3 11 24 30 20 10  4  0  0  0
     0  0  0  6 22 49 65 57 32 14  6  0  0     0  0  3 13 43 85100 75 38 16  6  0  0     0  0  2 12 38 64 65 43 20  8  3  0  0
     0  0  0  7 17 32 37 30 15  7  0  0  0     0  0  3 13 35 57 58 40 19  8  0  0  0     0  0  3 14 34 46 41 25 11  4  0  0  0
     0  0  0  2  4  5  5  3  2  0  0  0  0     0  0  0  3  6  9  9  5  2  0  0  0  0     0  0  0  4  7  8  6  3  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  5  6  7  5  3  0  0  0  0     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  7 13 19 17 10  5  3  0  0  0     0  0  0  0  4  4  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  9 15 16 13  7  3  2  0  0  0     0  0  0  0  7  4  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  8  5  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  8 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  4  6  0  0  0  0     0  0  0  2  5  9 11 11  7  5  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  3  4  5  6  0  0  0  0     0  0  0  3  9 18 24 23 15  7  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  4  0  0  0  0  0     0  0  0  3  6 12 16 14  9  5  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  4  4  4  2  0  0  0  0     0  0  0  0  3  5  7  5  3  0  0  0  0     0  0  0  0  3  4  4  3  2  0  0  0  0
     0  0  0  4 14 28 36 30 17  7  0  0  0     0  0  0  6 20 39 50 39 20  8  3  0  0     0  0  0  4 12 25 31 22 11  4  0  0  0
     0  0  2  9 30 59 74 62 35 15  6  0  0     0  0  3 14 44 84100 76 39 16  6  0  0     0  0  3 10 29 52 58 42 20  8  3  0  0
     0  0  0  6 19 38 47 37 20  8  0  0  0     0  0  3 10 30 54 61 43 21  9  0  0  0     0  0  3  8 22 34 35 23 11  4  0  0  0
     0  0  0  2  3  6  7  5  3  0  0  0  0     0  0  0  3  5  8  9  5  3  0  0  0  0     0  0  0  0  4  6  5  3  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  4  7  8  6  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  5  9 14 15 10  5  3  0  0  0     0  0  0  0  3  3  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  6  8 10  9  6  3  0  0  0  0     0  0  0  0  0  3  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  9 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  2  2  0  0  0  0  0     0  0  0  3  6 11 11  8  4  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  3  8 16 19 14  8  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  4  7  9  7  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4  7  8  6  3  0  0  0  0     0  0  0  3  6 12 14 10  5  2  0  0  0     0  0  0  2  4  9 12  9  4  0  0  0  0
     0  0  0  7 22 40 44 30 15  6  0  0  0     0  0  0  9 30 60 70 47 23  9  0  0  0     0  0  0  5 17 37 45 31 15  6  0  0  0
     0  0  0  7 30 60 68 48 24  9  0  0  0     0  0  0 10 43 89100 69 33 13  4  0  0     0  0  0  6 25 51 57 38 18  7  0  0  0
     0  0  0  3 12 26 31 21 10  4  0  0  0     0  0  0  4 18 38 43 27 12  5  0  0  0     0  0  0  3 10 20 23 14  6  3  0  0  0
     0  0  0  0  3  3  4  3  0  0  0  0  0     0  0  0  0  3  5  6  3  2  0  0  0  0     0  0  0  0  2  3  3  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  5  5  5  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  4 10 12  9  4  3  0  0  0     0  0  0  0  0  3  3  4  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 11 13  9  4  3  0  0  0     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  5  5  4  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  2  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 STANDARD DEVIATION (DEGREES) OF BEAM DIVERGENCE
 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1    0.024  0.025  0.026  0.024  0.023  0.025  0.025  0.025  0.024
     2    0.025  0.026  0.026  0.025  0.024  0.025  0.026  0.025  0.025
     3    0.024  0.026  0.026  0.025  0.023  0.024  0.025  0.025  0.025
     4    0.025  0.026  0.027  0.025  0.024  0.025  0.025  0.025  0.025
     5    0.024  0.026  0.026  0.025  0.023  0.025  0.026  0.025  0.024
     6    0.025  0.026  0.027  0.025  0.023  0.024  0.026  0.025  0.025
     7    0.025  0.026  0.027  0.025  0.023  0.024  0.026  0.026  0.025
     8    0.025  0.026  0.028  0.025  0.023  0.025  0.026  0.025  0.025
     9    0.025  0.027  0.029  0.026  0.023  0.025  0.027  0.026  0.025
    10    0.025  0.027  0.029  0.026  0.024  0.025  0.026  0.026  0.025
    11    0.025  0.027  0.029  0.026  0.024  0.025  0.027  0.027  0.025
    12    0.025  0.027  0.029  0.027  0.024  0.025  0.027  0.026  0.025
    13    0.025  0.027  0.029  0.027  0.024  0.026  0.027  0.027  0.025
    14    0.025  0.027  0.030  0.027  0.024  0.025  0.027  0.026  0.025
    15    0.026  0.028  0.031  0.027  0.024  0.026  0.027  0.026  0.025
    16    0.025  0.029  0.031  0.028  0.025  0.026  0.027  0.027  0.025
    17    0.026  0.029  0.032  0.028  0.024  0.025  0.027  0.027  0.025
    18    0.025  0.029  0.032  0.029  0.025  0.026  0.027  0.027  0.025
    19    0.026  0.030  0.033  0.029  0.024  0.025  0.027  0.026  0.025
    20    0.025  0.029  0.033  0.029  0.025  0.026  0.027  0.027  0.025
    21    0.026  0.030  0.033  0.030  0.025  0.025  0.026  0.027  0.026
    22    0.026  0.030  0.033  0.030  0.025  0.025  0.026  0.026  0.025
    23    0.026  0.029  0.033  0.030  0.025  0.025  0.026  0.026  0.025
    24    0.026  0.030  0.033  0.030  0.025  0.025  0.026  0.026  0.026
    25    0.026  0.029  0.032  0.029  0.025  0.025  0.025  0.026  0.025


 STANDARD DEVIATION (DEGREES) OF REFLECTING RANGE
 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1    0.031  0.032  0.035  0.033  0.030  0.026  0.029  0.031  0.030
     2    0.035  0.027  0.034  0.035  0.033  0.038  0.033  0.032  0.030
     3    0.032  0.038  0.036  0.033  0.029  0.023  0.029  0.032  0.033
     4    0.032  0.028  0.033  0.036  0.033  0.027  0.032  0.034  0.030
     5    0.031  0.028  0.033  0.034  0.029  0.023  0.028  0.032  0.030
     6    0.028  0.026  0.031  0.033  0.028  0.022  0.028  0.030  0.028
     7    0.030  0.028  0.033  0.033  0.029  0.024  0.029  0.033  0.031
     8    0.031  0.027  0.032  0.035  0.031  0.024  0.028  0.031  0.029
     9    0.033  0.042  0.037  0.034  0.030  0.025  0.029  0.030  0.033
    10    0.029  0.026  0.031  0.034  0.030  0.024  0.029  0.031  0.029
    11    0.030  0.025  0.031  0.034  0.029  0.025  0.029  0.031  0.030
    12    0.029  0.026  0.031  0.033  0.028  0.025  0.030  0.030  0.029
    13    0.030  0.025  0.030  0.032  0.030  0.031  0.033  0.031  0.029
    14    0.031  0.026  0.029  0.032  0.028  0.025  0.030  0.033  0.030
    15    0.031  0.026  0.030  0.032  0.030  0.030  0.032  0.031  0.030
    16    0.029  0.027  0.029  0.031  0.029  0.029  0.032  0.031  0.029
    17    0.031  0.026  0.030  0.032  0.030  0.028  0.033  0.034  0.031
    18    0.031  0.025  0.029  0.034  0.031  0.033  0.037  0.032  0.029
    19    0.032  0.030  0.030  0.034  0.033  0.030  0.034  0.032  0.030
    20    0.032  0.026  0.028  0.030  0.031  0.042  0.042  0.033  0.029
    21    0.033  0.026  0.027  0.030  0.031  0.042  0.041  0.034  0.031
    22    0.032  0.027  0.028  0.032  0.031  0.034  0.038  0.035  0.031
    23    0.032  0.028  0.029  0.030  0.031  0.039  0.042  0.035  0.030
    24    0.034  0.027  0.029  0.035  0.033  0.036  0.041  0.034  0.030
    25    0.034  0.030  0.030  0.032  0.030  0.036  0.041  0.037  0.032


 ***** SUGGESTED VALUES FOR INPUT PARAMETERS *****
 BEAM_DIVERGENCE=   0.26208  BEAM_DIVERGENCE_E.S.D.=   0.02621
 REFLECTING_RANGE=  0.21695  REFLECTING_RANGE_E.S.D.=  0.03099



 NUMBER OF STRONG REFLECTIONS USED FOR LEARNING REFERENCE PROFILES

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1       45     26     53     67     56     27     54     68     52
     2       47     33     59     66     49     27     51     57     53
     3       51     29     53     67     58     29     59     76     61
     4       47     34     64     72     51     26     50     61     54
     5       54     30     57     66     59     32     58     71     57
     6       47     33     57     66     55     28     50     61     55
     7       45     31     56     64     55     30     55     65     57
     8       50     33     57     70     50     26     52     69     61
     9       43     27     54     66     55     31     55     64     51
    10       48     33     60     65     54     28     53     67     56
    11       51     35     60     73     60     31     54     65     59
    12       47     31     53     64     51     26     52     68     60
    13       49     30     55     70     60     29     54     69     61
    14       50     32     59     67     52     30     52     66     55
    15       48     30     52     62     55     29     53     62     53
    16       50     31     54     67     52     30     56     69     59
    17       50     31     58     69     57     30     53     67     60
    18       50     32     55     69     58     32     50     67     52
    19       49     27     52     68     56     31     54     64     56
    20       50     34     54     66     55     33     53     66     64
    21       46     29     56     64     54     32     60     66     45
    22       47     32     56     66     55     29     58     72     66
    23       51     33     61     70     60     30     48     63     52
    24       50     38     58     65     54     32     56     68     63
    25       52     29     53     67     57     30     57     63     49
 margin    1215    780   1406   1676   1379    739   1350   1657   1413



 NUMBER OF GRID POINTS IN REFERENCE PROFILES WITH  SIGNAL >   2.% OF MAXIMUM

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1      133    131    140    130    100    108    114    121    122
     2      151    130    140    139    118    114    134    133    125
     3      134    135    142    130    104    103    114    124    126
     4      130    118    144    136    120    102    117    134    122
     5      129    124    138    134    110    104    124    130    121
     6      121    121    140    134    106    101    108    125    122
     7      131    124    134    128    106    103    122    139    124
     8      132    130    148    136    112    102    109    121    122
     9      137    139    161    138    110    109    118    127    133
    10      125    129    152    138    113    106    114    127    115
    11      124    129    148    140    116    105    126    138    122
    12      124    132    155    140    113    107    125    129    115
    13      125    130    142    136    113    115    131    137    114
    14      136    131    146    143    114    105    122    137    120
    15      138    134    154    144    118    121    136    132    116
    16      130    144    154    146    120    121    132    132    118
    17      141    148    159    146    124    118    130    140    119
    18      135    142    172    159    125    127    134    135    116
    19      142    155    170    160    128    127    134    127    119
    20      135    143    163    146    123    137    138    137    117
    21      144    148    161    154    124    133    139    137    122
    22      144    145    175    162    132    133    140    135    120
    23      141    137    166    154    119    136    140    135    119
    24      151    144    174    174    135    140    140    141    121
    25      151    148    173    164    122    135    141    137    124
 margin     167    193    209    186    148    164    175    169    156

 MINIMUM INTEGRATION VOLUME IS  100 GRID POINTS



 CORRELATIONS BETWEEN REFERENCE PROFILES WITH CENTRAL PROFILE #1

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1     1.00   0.95   0.96   0.97   0.98   0.96   0.98   0.98   0.98
     2     1.00   0.95   0.96   0.95   0.97   0.84   0.92   0.96   0.95
     3     1.00   0.85   0.92   0.97   0.96   0.95   0.96   0.98   0.94
     4     1.00   0.92   0.94   0.98   0.95   0.93   0.93   0.95   0.96
     5     1.00   0.98   0.95   0.95   0.98   0.95   0.96   0.96   0.97
     6     1.00   0.97   0.92   0.97   0.97   0.97   0.95   0.97   0.93
     7     1.00   0.96   0.90   0.98   0.98   0.96   0.95   0.94   0.97
     8     1.00   0.94   0.94   0.95   0.93   0.95   0.96   0.98   0.93
     9     1.00   0.68   0.82   0.98   0.96   0.95   0.93   0.96   0.91
    10     1.00   0.95   0.94   0.97   0.97   0.97   0.95   0.96   0.95
    11     1.00   0.96   0.94   0.96   0.99   0.93   0.95   0.94   0.94
    12     1.00   0.96   0.92   0.96   0.97   0.96   0.95   0.97   0.97
    13     1.00   0.97   0.94   0.96   0.98   0.91   0.90   0.97   0.97
    14     1.00   0.96   0.93   0.98   0.98   0.95   0.95   0.96   0.96
    15     1.00   0.95   0.91   0.97   0.98   0.91   0.94   0.97   0.94
    16     1.00   0.93   0.92   0.97   0.99   0.97   0.96   0.98   0.95
    17     1.00   0.95   0.88   0.96   0.98   0.96   0.94   0.94   0.94
    18     1.00   0.94   0.91   0.95   0.98   0.90   0.88   0.97   0.94
    19     1.00   0.86   0.86   0.96   0.94   0.96   0.94   0.97   0.93
    20     1.00   0.93   0.88   0.96   0.98   0.85   0.88   0.96   0.96
    21     1.00   0.91   0.88   0.96   0.96   0.85   0.91   0.96   0.94
    22     1.00   0.91   0.85   0.97   0.99   0.94   0.92   0.97   0.97
    23     1.00   0.91   0.89   0.96   0.98   0.90   0.85   0.96   0.97
    24     1.00   0.88   0.85   0.91   0.93   0.96   0.90   0.97   0.92
    25     1.00   0.93   0.87   0.95   0.98   0.94   0.92   0.96   0.95
 margin    1.00   0.97   0.95   0.98   0.99   0.98   0.97   0.98   0.97



 CORRELATIONS BETWEEN REFERENCE PROFILES WITH RUN-AVERAGED PROFILES

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1     0.99   0.93   0.93   0.98   0.99   0.98   0.96   0.97   0.96
     2     0.99   0.95   0.96   0.98   0.95   0.84   0.93   0.97   0.99
     3     0.99   0.86   0.93   0.98   0.99   0.96   0.97   0.98   0.94
     4     0.99   0.98   0.98   0.98   0.96   0.83   0.89   0.96   0.99
     5     1.00   0.96   0.96   0.99   1.00   0.98   0.97   0.98   0.99
     6     0.99   0.99   0.98   0.99   0.99   0.97   0.97   0.99   0.99
     7     0.99   0.99   0.97   0.98   0.99   0.98   0.96   0.98   0.99
     8     0.99   0.99   0.98   0.99   0.98   0.96   0.97   0.99   1.00
     9     0.99   0.76   0.90   0.99   0.98   0.93   0.97   0.99   0.93
    10     1.00   0.99   0.99   0.99   0.99   0.97   0.98   1.00   1.00
    11     1.00   0.95   0.99   0.99   1.00   0.99   0.98   0.98   0.97
    12     1.00   0.99   0.99   0.99   1.00   0.99   0.99   0.99   1.00
    13     0.99   0.99   0.99   1.00   0.98   0.94   0.95   0.98   0.99
    14     0.99   0.99   0.99   0.99   0.99   0.99   0.99   0.99   0.99
    15     0.99   0.94   0.98   0.99   1.00   0.97   0.99   0.99   0.97
    16     1.00   0.98   0.98   0.99   0.99   0.99   0.99   1.00   0.99
    17     0.99   0.99   0.98   1.00   1.00   0.97   0.95   0.97   0.97
    18     0.99   0.96   0.96   0.99   0.99   0.96   0.96   0.98   1.00
    19     0.99   0.95   0.97   0.98   0.96   0.96   0.98   0.99   0.98
    20     0.99   0.98   0.97   0.97   0.94   0.85   0.90   0.98   1.00
    21     0.99   0.96   0.96   0.97   0.95   0.87   0.94   0.98   0.99
    22     0.99   0.96   0.95   0.98   0.99   0.98   0.97   0.98   0.98
    23     0.99   0.97   0.96   0.98   0.99   0.92   0.90   0.96   0.96
    24     1.00   0.96   0.94   0.96   0.96   0.98   0.95   0.97   0.97
    25     0.99   0.91   0.93   0.98   0.99   0.97   0.93   0.96   0.96
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

       -1... 0      -0.135      0.062        -2.         12
        0... 1       0.220      0.102         5.         75
        1... 2       0.420      0.100        17.        129
        2... 3       0.573      0.079        34.        161
        3... 6       0.691      0.089        88.        770
        6... 9       0.745      0.078       209.       1031
        9...12       0.756      0.074       378.       1161
       12...15       0.761      0.077       604.       1186
       15...18       0.773      0.073       880.       1043
       18...21       0.770      0.074      1221.        930
       21...24       0.772      0.076      1616.        788
       24...27       0.777      0.077      2107.        688
       27...30       0.779      0.077      2613.        619
       30...33       0.779      0.079      3459.        560
       33...36       0.773      0.080      4464.        550
       36...39       0.772      0.081      5895.        484
       39...42       0.772      0.080      8049.        555
       42...45       0.769      0.081     11906.        652
       45...48       0.764      0.088     21831.        887
       48...51       0.758      0.102     49653.        715

   12996 REFLECTIONS SAVED ON FILE "INTEGRATE.HKL"


 total elapsed wall-clock time for INTEGRATE      15.0 sec

 TASK     cpu time (sec)    elapsed wall-clock time (sec)
   1         246.9                    14.9
