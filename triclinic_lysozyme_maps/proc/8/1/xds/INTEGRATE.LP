
 ***** INTEGRATE ***** (VERSION Mar 15, 2019  BUILT=20190806)   6-Jan-2020

 INPUT PARAMETER VALUES
 ----------------------
 MAXIMUM_NUMBER_OF_PROCESSORS=   0  (0: automatic choice of OpenMP threads)
 number of OpenMP threads used  20
 MAXIMUM_NUMBER_OF_JOBS=   0  (0:automatic choice of forked main tasks)
 number of forked INTEGRATE tasks   1
 CLUSTER_NODES=
 SECONDS=     0
 NAME_TEMPLATE_OF_DATA_FRAMES=~/data/lysozyme/nitrate/ubatch4/lys_nitr_8_1_????.cbf CBF
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
 UNIT_CELL_CONSTANTS=    27.43    32.13    34.52  88.633 108.435 111.905
 ROTATION_AXIS= -0.99999 -0.00285 -0.00132
 OSCILLATION_RANGE=  0.100000 DEGREES
 STARTING_ANGLE= 300.000 STARTING_FRAME=       1
 X-RAY_WAVELENGTH=   0.97700 ANGSTROM
 INCIDENT_BEAM_DIRECTION=  0.000718  0.005368  1.023527
 NUMBER OF DETECTOR SEGMENTS   1
 NX=  2463 NY=  2527 QX= 0.17200 QY= 0.17200
 ORGX=   1265.40  ORGY=   1237.08
 DETECTOR_DISTANCE=   231.326
 DIRECTION_OF_DETECTOR_X-AXIS=   1.00000   0.00000   0.00000
 DIRECTION_OF_DETECTOR_Y-AXIS=   0.00000   1.00000   0.00000
 SEGMENT=    1 2463    1 2527  ! x1,x2,y1,y2 (pixel units)
    SEGMENT_ORGX=      0.00   SEGMENT_ORGY=      0.00
    SEGMENT_DISTANCE=     0.00
    DIRECTION_OF_SEGMENT_X-AXIS=  1.000000  0.000000  0.000000
    DIRECTION_OF_SEGMENT_Y-AXIS=  0.000000  1.000000  0.000000
 RESOLUTION RANGE RECORDED BY DETECTOR (ANGSTROM)    44.459     1.212
 NUMBER OF TRUSTED DETECTOR PIXELS    5175365
 MEAN CONTENTS OF TRUSTED PIXELS IN BACKGROUND TABLE    117.543
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
     1   0  1.089  5167021    0    154      90     4  0.04647  0.10853
     2   0  1.059  5166072    0    161     106     4  0.04467  0.10813
     3   0  1.059  5165810    0    157     105     2  0.04297  0.11017
     4   0  1.053  5166651    0    153      92     1  0.04314  0.13047
     5   0  1.055  5167525    0    144      85     2  0.04700  0.14455
     6   0  1.055  5166802    0    130      85     5  0.04785  0.10945
     7   0  1.049  5166983    0    130      88     5  0.04737  0.12819
     8   0  1.052  5168200    0    128      77     6  0.04690  0.12705
     9   0  1.051  5167661    0    136      81     7  0.04672  0.12243
    10   0  1.053  5167524    0    140      83     4  0.04792  0.13264
    11   0  1.051  5166996    0    138      88     9  0.04276  0.07740
    12   0  1.054  5165840    0    142     102    11  0.04111  0.08248
    13   0  1.048  5166422    0    153      97     7  0.04334  0.10812
    14   0  1.051  5167418    0    156      88     5  0.04064  0.12053
    15   0  1.051  5167631    0    154      83     9  0.04254  0.09389
    16   0  1.051  5167608    0    152      87     5  0.04393  0.13175
    17   0  1.049  5167681    0    149      86     2  0.04478  0.16921
    18   0  1.049  5167383    0    153      86     5  0.04657  0.15196
    19   0  1.048  5168099    0    135      76     8  0.04219  0.12280
    20   0  1.050  5167731    0    147      85     9  0.04039  0.12826

   478 OUT OF    481 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.44
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.425    32.126    34.501  88.639 108.407 111.901
 REC. CELL PARAMETERS   0.041646  0.033741  0.030723  83.868  70.678  67.336
 COORDINATES OF UNIT CELL A-AXIS   -14.337   -23.355    -1.051
 COORDINATES OF UNIT CELL B-AXIS    31.417    -5.395     3.992
 COORDINATES OF UNIT CELL C-AXIS    -1.194    12.072    32.298
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.010    -0.038     0.005
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.040
 LAB COORDINATES OF ROTATION AXIS -0.999995 -0.002850 -0.001317
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000835  0.005172  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.39   1244.04
 DETECTOR ORIGIN (PIXELS) AT                     1265.29   1237.24
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.32
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     478 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   55     33     61     66     52     29     56     64     56
 SIGMAB (degree)       0.033  0.033  0.034  0.033  0.032  0.032  0.031  0.032  0.032
 SIGMAR (degree)       0.031  0.028  0.030  0.032  0.030  0.026  0.029  0.032  0.029


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    478 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  2  4  4  2  1  0  0  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  1  5 12 15 11  5  2  1  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  2  3  3  1  1  0  0  0     0  0  0  1  7 19 25 19  9  4  2  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  1  5 11 14 11  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  3  4  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  6 10 13  7  3  1  0  0  0     0  0  0  2  8 13 16  8  3  1  0  0  0     0  0  0  1  5  7  8  4  2  1  0  0  0
     0  0  0  4 19 43 50 34 16  6  2  1  0     0  0  1  7 30 63 66 40 19  7  3  1  0     0  0  0  4 17 32 32 19  8  3  1  0  0
     0  0  0  5 26 66 78 54 26 12  4  1  0     0  0  1  9 44 95100 65 31 13  5  1  0     0  0  1  6 25 48 47 30 14  6  2  1  0
     0  0  0  3 15 37 45 31 15  7  2  1  0     0  0  0  6 26 55 60 41 18  8  3  1  0     0  0  0  4 16 29 31 21  9  4  1  0  0
     0  0  0  1  4  9 11  7  4  1  0  0  0     0  0  0  1  7 14 16 10  5  2  1  0  0     0  0  0  1  4  8  9  5  3  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  6  5  3  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  9  8  5  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  6  7  4  2  1  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  1  0  0  0  0  0  0  0  0  0  0  0     0  1  0  0  0  0  0  0  0  0  0  0  0

    667 REFLECTION INTENSITIES INTEGRATED


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
    21   0  1.050  5167747    0    133      85     9  0.04264  0.09570
    22   0  1.053  5167632    0    137      87     9  0.04597  0.10900
    23   0  1.049  5167760    0    137      82     8  0.04329  0.09626
    24   0  1.047  5167340    0    142      90     5  0.04522  0.11881
    25   0  1.046  5167060    0    147      90     7  0.04348  0.09718
    26   0  1.047  5167698    0    149      80     5  0.04511  0.11212
    27   0  1.045  5166739    0    147      90     4  0.04538  0.13311
    28   0  1.046  5166123    0    140      99     4  0.04420  0.13521
    29   0  1.050  5167348    0    141      89     7  0.04636  0.10349
    30   0  1.046  5167488    0    142      83     9  0.04133  0.09246
    31   0  1.046  5166840    0    140      89     9  0.04330  0.07166
    32   0  1.049  5167054    0    133      87     9  0.04673  0.08093
    33   0  1.050  5167776    0    124      79     8  0.04411  0.07996
    34   0  1.048  5168077    0    125      76     8  0.04336  0.08901
    35   0  1.050  5168241    0    134      75     5  0.04806  0.10424
    36   0  1.048  5167517    0    130      89     4  0.04574  0.11640
    37   0  1.052  5167955    0    141      77     5  0.04346  0.09297
    38   0  1.051  5167228    0    146      84     9  0.04351  0.08893
    39   0  1.052  5167658    0    150      84     7  0.04379  0.10649
    40   0  1.052  5168165    0    140      76     5  0.04538  0.13593

   497 OUT OF    498 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.57
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.17
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.426    32.127    34.507  88.639 108.424 111.901
 REC. CELL PARAMETERS   0.041648  0.033741  0.030722  83.861  70.660  67.334
 COORDINATES OF UNIT CELL A-AXIS   -14.335   -23.358    -1.058
 COORDINATES OF UNIT CELL B-AXIS    31.418    -5.392     3.992
 COORDINATES OF UNIT CELL C-AXIS    -1.194    12.075    32.303
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.010    -0.042     0.001
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.040
 LAB COORDINATES OF ROTATION AXIS -0.999995 -0.002850 -0.001317
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000748  0.005121  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.40   1244.10
 DETECTOR ORIGIN (PIXELS) AT                     1265.41   1237.37
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.30
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     497 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     32     58     71     61     29     53     68     58
 SIGMAB (degree)       0.033  0.033  0.034  0.033  0.032  0.032  0.031  0.032  0.032
 SIGMAR (degree)       0.034  0.029  0.031  0.033  0.031  0.027  0.029  0.032  0.031


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    497 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  1  0  0  0  0     0  0  0  0  1  3  4  3  1  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  1  1  0     0  0  0  1  4 11 15 11  5  2  1  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  4  3  2  1  0  0  0     0  0  0  1  6 17 22 17  9  4  1  0  0
     0  0  0  0  0  0  1  0  1  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0     0  0  0  0  4 10 13 10  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0 -1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0 -1  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  5 11 13  7  4  1  0  0  0     0  0  0  1  7 16 17 10  4  1  0  0  0     0  0  0  1  5 10  9  6  2  1  0  0  0
     0  0  0  3 15 40 50 33 16  7  2  1  0     0  0  1  6 27 62 69 43 21  9  3  1  0     0  0  0  4 18 36 39 23 11  5  2  0  0
     0  0  0  4 24 59 69 49 24 11  4  1  0     0  0  1  8 43 94100 66 31 14  5  1  0     0  0  1  7 29 57 58 36 17  7  3  1  0
     0  0  0  2 14 36 40 28 14  6  2  0  0     0  0  0  5 26 59 61 39 19  8  3  1  0     0  0  0  4 17 38 37 22 11  4  2  0  0
     0  0  0  1  3  7  9  7  3  1  0  0  0     0  0  0  1  6 12 15 10  5  2  1  0  0     0  0  0  1  4  8 10  7  3  2  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  1  1  2  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  3  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  9  9  5  2  1  0  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0    -1  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  2  7 14 13  8  4  2  1  0  0     0  0  0  1  1  2  1  1  0  0  0  0  0     0  0  0  1  1  1  1  1  0  0  0  0  0
     0  0  0  1  4 10  9  6  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  1  1  1  0  0  0  0  0  0  0
     0  0  0  0  1  2  3  2  1  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    635 REFLECTION INTENSITIES INTEGRATED


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
    41   0  1.053  5167462    0    144      82     6  0.04596  0.12423
    42   0  1.052  5166830    0    139      90     9  0.04364  0.10568
    43   0  1.052  5166886    0    143      90     7  0.04401  0.09665
    44   0  1.052  5167879    0    140      81     9  0.04360  0.10949
    45   0  1.051  5168519    0    135      70     3  0.04585  0.10810
    46   0  1.051  5168045    0    131      75     1  0.04574  0.11195
    47   0  1.050  5166370    0    132      95     9  0.04446  0.09366
    48   0  1.055  5166719    0    148      90     9  0.04593  0.08869
    49   0  1.049  5166634    0    157      91    10  0.04474  0.07880
    50   0  1.046  5165377    0    156     107     8  0.04472  0.10897
    51   0  1.048  5166148    0    154     100     7  0.04564  0.09818
    52   0  1.051  5166324    0    137      95     6  0.04508  0.10203
    53   0  1.051  5166479    0    129      95     7  0.04514  0.10577
    54   0  1.048  5167769    0    134      81     6  0.04408  0.09890
    55   0  1.046  5166795    0    134      90     5  0.04336  0.09411
    56   0  1.050  5166731    0    135      95     5  0.04651  0.13145
    57   0  1.049  5166438    0    147      99     8  0.04465  0.14721
    58   0  1.049  5166289    0    150     101    11  0.04504  0.13631
    59   0  1.047  5166280    0    158     100    10  0.04748  0.13797
    60   0  1.052  5165479    0    159     110    11  0.04620  0.14810

   483 OUT OF    485 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.50
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.16
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.429    32.124    34.510  88.630 108.431 111.911
 REC. CELL PARAMETERS   0.041647  0.033747  0.030720  83.865  70.656  67.327
 COORDINATES OF UNIT CELL A-AXIS   -14.338   -23.359    -1.060
 COORDINATES OF UNIT CELL B-AXIS    31.416    -5.388     3.992
 COORDINATES OF UNIT CELL C-AXIS    -1.190    12.076    32.306
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.012    -0.043     0.002
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.040
 LAB COORDINATES OF ROTATION AXIS -0.999995 -0.002850 -0.001317
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000716  0.005128  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.33   1244.10
 DETECTOR ORIGIN (PIXELS) AT                     1265.39   1237.37
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.31
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     483 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   52     27     48     62     51     30     57     77     58
 SIGMAB (degree)       0.033  0.034  0.034  0.033  0.033  0.033  0.031  0.032  0.032
 SIGMAR (degree)       0.033  0.029  0.031  0.033  0.031  0.028  0.030  0.032  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    483 STRONG REFLECTIONS *****


     0 -1  0  0  0  0  0  0  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     1  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  3  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  1  4 11 15 12  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0 -1     0  0  0  0  1  3  4  3  2  1  0  0  0     0  0  0  1  7 19 25 20 10  4  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  0  5 13 17 12  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  3  3  3  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  2  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  5 10 12  7  3  1  0  0  0     0  0  0  2  8 15 16 10  4  2  0  0  0     0  0  0  1  5  9  9  6  2  1  0  0  0
     0  0  0  3 17 40 48 34 15  7  3  1  0     0  0  1  6 29 61 67 44 19  9  3  1  0     0  0  0  5 19 37 38 24 10  4  1  0  0
     0  0  0  5 25 63 73 52 26 11  4  1  0     0  0  1  9 41 93100 65 32 13  5  2  0     0  0  1  7 28 57 57 35 17  7  3  1  0
     0  0  0  2 15 38 44 29 14  6  2  1  0     0  0  1  5 26 55 63 38 18  8  3  1  0     0  0  1  4 19 36 38 21 10  5  2  0  0
     0  0  0  0  3  7  9  6  3  1  0  0  0     0  0  0  1  7 14 15 10  5  2  0  0  0     0  0  0  1  6 11 11  7  4  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  2  2  2  0  0  0  0  0     0  0  0  0  1  1  2  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  3  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0
     0  0  0  1  5  9 10  6  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0 -1     0  0  0  0  1  0  0  0  0  0  0  0  0
     0  0  0  2  8 15 15  9  4  2  1  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  1  5  9 10  6  3  1  0  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  2  3  3  2  1  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0 -1  0  0  0  0  0  0  0  0  0  0  0     0 -1  0  0  0  0  0  0  0  0  0  0  0

    658 REFLECTION INTENSITIES INTEGRATED


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
    61   0  1.050  5165662    0    156     107    11  0.04671  0.13110
    62   0  1.052  5166708    0    147      95     6  0.04928  0.13028
    63   0  1.050  5167006    0    143      91     6  0.04500  0.10936
    64   0  1.050  5166508    0    142      93     4  0.04599  0.13208
    65   0  1.047  5166291    0    136      97     1  0.04565  0.13300
    66   0  1.053  5166599    0    136      94     2  0.04452  0.14046
    67   0  1.050  5167953    0    131      82     5  0.04537  0.14241
    68   0  1.048  5168225    0    132      77     4  0.04696  0.15784
    69   0  1.054  5167803    0    124      83     5  0.04598  0.18723
    70   0  1.052  5167168    0    139      89     9  0.04717  0.11782
    71   0  1.052  5167477    0    136      89     7  0.04580  0.12786
    72   0  1.055  5167097    0    145      85     9  0.04563  0.12764
    73   0  1.053  5166461    0    152      96     9  0.04588  0.11111
    74   0  1.054  5166346    0    153      94     9  0.04347  0.08033
    75   0  1.055  5166353    0    163     101     6  0.04371  0.09876
    76   0  1.056  5166100    0    150     102     7  0.04430  0.09412
    77   0  1.057  5166645    0    151      92     4  0.04436  0.10166
    78   0  1.052  5167597    0    151      83     5  0.04659  0.09301
    79   0  1.054  5166557    0    145      91     4  0.04304  0.08426
    80   0  1.060  5166273    0    136      99     5  0.04490  0.09040

   492 OUT OF    493 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.55
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.16
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.426    32.121    34.498  88.606 108.418 111.926
 REC. CELL PARAMETERS   0.041651  0.033751  0.030727  83.892  70.677  67.320
 COORDINATES OF UNIT CELL A-AXIS   -14.337   -23.357    -1.052
 COORDINATES OF UNIT CELL B-AXIS    31.414    -5.378     4.000
 COORDINATES OF UNIT CELL C-AXIS    -1.187    12.073    32.295
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.006    -0.034    -0.004
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.040
 LAB COORDINATES OF ROTATION AXIS -0.999995 -0.002850 -0.001317
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000589  0.005263  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.32   1244.09
 DETECTOR ORIGIN (PIXELS) AT                     1265.54   1237.17
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.29
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     492 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   51     33     66     73     57     30     50     60     54
 SIGMAB (degree)       0.033  0.034  0.034  0.033  0.032  0.032  0.032  0.032  0.032
 SIGMAR (degree)       0.034  0.032  0.031  0.033  0.031  0.028  0.031  0.032  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    492 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  2  4  4  3  1  1  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  1  5 13 17 11  6  3  1  0  0
     0  0  0  0  0  0  1  1  1  0  0  0  0     0  0  0  0  1  3  4  3  2  1  0  0  0     0  0  0  1  7 18 26 21 11  4  2  1  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  1  4 12 17 13  7  3  1  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  0  1  3  3  3  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  6 13 13  8  4  2  0  0  0     0  0  0  2  8 18 16 10  5  2  0  0  0     0  0  0  1  4  9  8  5  2  1  0  0  0
     0  0  0  3 17 43 52 35 18  7  3  1  0     0  0  1  6 28 62 66 42 20  9  3  1  0     0  0  0  4 17 35 34 22 10  4  2  0  0
     0  0  0  4 25 61 77 52 27 11  5  1  0     0  0  0  8 42 91100 62 32 13  5  1  0     0  0  1  6 26 52 53 31 16  7  2  0  0
     0  0  0  3 14 37 48 32 15  7  2  1  0     0  0  0  6 24 55 65 40 19  8  3  1  0     0  0  0  5 16 33 37 23 10  4  2  0  0
     0  0  0  1  4  9  9  7  4  1  0  0  0     0  0  0  2  7 14 15 11  5  2  1  0  0     0  0  0  2  5 10 10  7  3  1  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  1  1  2  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  1  0  1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  1  4  8  7  5  2  1  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  1  1  1  0  1  0  0  0  0
     0  0  0  2  6 11 12  7  4  1  0  0  0     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  0  1  0  1  0  0  0  0  0
     0  0  0  1  4  8  9  5  2  1  0  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  1  2  3  3  2  1  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    660 REFLECTION INTENSITIES INTEGRATED


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
    81   0  1.063  5166934    0    145      91     5  0.04467  0.08817
    82   0  1.061  5166584    0    145      93     3  0.04505  0.11119
    83   0  1.062  5166184    0    149     100     6  0.04190  0.09410
    84   0  1.063  5166475    0    151      94     8  0.04321  0.08940
    85   0  1.059  5167466    0    127      88     9  0.04452  0.10137
    86   0  1.058  5168421    0    137      79     8  0.04441  0.12847
    87   0  1.056  5167977    0    127      79     7  0.04599  0.10459
    88   0  1.058  5167347    0    131      91     7  0.04215  0.12538
    89   0  1.058  5167420    0    136      86     2  0.04755  0.11493
    90   0  1.059  5167509    0    131      83     2  0.04338  0.11182
    91   0  1.056  5167789    0    129      79     7  0.04132  0.07538
    92   0  1.057  5167909    0    133      81     7  0.04276  0.08492
    93   0  1.061  5168158    0    135      80     8  0.04606  0.11234
    94   0  1.058  5168300    0    141      73     3  0.04666  0.14915
    95   0  1.061  5166574    0    137      96     2  0.04698  0.14740
    96   0  1.063  5166304    0    146      97     2  0.04468  0.17018
    97   0  1.066  5166024    0    142     100    10  0.04647  0.12399
    98   0  1.062  5167187    0    145      92    10  0.04679  0.11337
    99   0  1.065  5167962    0    146      83     9  0.04440  0.13281
   100   0  1.056  5167927    0    135      76     8  0.04178  0.11572

   472 OUT OF    473 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.48
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.15
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.434    32.124    34.503  88.612 108.441 111.914
 REC. CELL PARAMETERS   0.041641  0.033746  0.030727  83.880  70.651  67.328
 COORDINATES OF UNIT CELL A-AXIS   -14.336   -23.367    -1.068
 COORDINATES OF UNIT CELL B-AXIS    31.417    -5.381     3.997
 COORDINATES OF UNIT CELL C-AXIS    -1.188    12.067    32.302
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.014    -0.058    -0.006
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.040
 LAB COORDINATES OF ROTATION AXIS -0.999995 -0.002850 -0.001317
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000472  0.005007  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.32   1244.10
 DETECTOR ORIGIN (PIXELS) AT                     1265.70   1237.53
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.30
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     472 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   49     32     49     65     55     28     53     67     57
 SIGMAB (degree)       0.032  0.033  0.033  0.033  0.032  0.033  0.032  0.032  0.032
 SIGMAR (degree)       0.036  0.030  0.031  0.033  0.032  0.031  0.032  0.032  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    472 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
    -1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  3  4  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  1  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  1  4 11 15 12  6  2  1  0  0
     0  0  0  0  1  1  1  0  0  0  0  0 -1     0  0  0  0  1  2  4  4  2  1  0  0  0     0  0  0  1  7 18 25 21 12  5  2  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  3  1  1  0  0  0     0  0  0  0  4 12 17 13  7  3  1  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  4  5  3  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
    -1  0  0  0 -1  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0 -1 -1  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0 -1  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  3 10 11  8  3  1  0  0  0     0  0  0  2  7 16 16 10  4  2  0  0  0     0  0  0  1  4 10 10  6  2  1  0  0  0
     0  0  0  3 14 38 49 35 17  7  2  1  0     0  0  0  5 27 61 68 44 20  9  3  1  0     0  0  0  4 19 39 41 26 11  5  1  0  0
     0  0  0  4 23 60 72 54 29 12  4  1  0     0  0  1  8 43 94100 68 32 14  5  1  0     0  0  1  7 32 60 60 40 18  7  3  0  0
     0  0  0  2 14 37 42 30 14  6  2  0  0     0  0  0  4 25 59 59 38 18  8  3  1  0     0  0  0  4 19 40 39 24 11  5  2  0  0
     0  0  0  0  3 10 11  7  3  1  0  0  0     0  0  0  1  6 16 17 11  4  2  0  0  0     0  0  0  1  5 12 13  8  3  1  0  0  0
     0  0  0  0  0  1  2  1  0  0  0  0  0     0  0  0  0  1  2  3  1  1  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  4  3  2  1  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  6 12 11  7  3  1  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  2 10 17 17 11  5  2  1  0  0     0  0  0  1  2  2  3  2  1  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0
     0  0  0  2  6 13 12  7  3  1  1  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  2  4  5  3  1  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    635 REFLECTION INTENSITIES INTEGRATED


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
   101   0  1.058  5167617    0    145      88     7  0.04280  0.08887
   102   0  1.058  5167185    0    148      91     5  0.04180  0.09368
   103   0  1.056  5166668    0    144      94     8  0.03877  0.07103
   104   0  1.055  5167500    0    147      90     4  0.04030  0.09347
   105   0  1.056  5167836    0    141      88     2  0.03920  0.14495
   106   0  1.058  5168035    0    130      82     8  0.04227  0.10125
   107   0  1.058  5167778    0    136      85     9  0.04300  0.12027
   108   0  1.062  5167238    0    134      93    10  0.03884  0.10456
   109   0  1.058  5166637    0    136      98    10  0.04340  0.09488
   110   0  1.060  5167460    0    142      88     9  0.04349  0.09507
   111   0  1.060  5167956    0    139      82     9  0.04391  0.10363
   112   0  1.062  5167566    0    134      83     9  0.04111  0.06670
   113   0  1.060  5167228    0    140      84     5  0.04323  0.07486
   114   0  1.063  5167426    0    143      89     4  0.04299  0.10443
   115   0  1.062  5167181    0    149      97     8  0.04345  0.10089
   116   0  1.064  5167429    0    141      84     4  0.04215  0.09725
   117   0  1.068  5167637    0    142      84     8  0.04358  0.07895
   118   0  1.068  5167500    0    139      91     3  0.04348  0.11201
   119   0  1.066  5167587    0    134      87     8  0.04624  0.09436
   120   0  1.066  5167691    0    139      87     8  0.04458  0.09930

   506 OUT OF    509 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.52
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.16
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.428    32.123    34.505  88.621 108.434 111.908
 REC. CELL PARAMETERS   0.041648  0.033746  0.030724  83.874  70.656  67.332
 COORDINATES OF UNIT CELL A-AXIS   -14.336   -23.360    -1.054
 COORDINATES OF UNIT CELL B-AXIS    31.415    -5.387     3.995
 COORDINATES OF UNIT CELL C-AXIS    -1.187    12.082    32.299
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.011    -0.031    -0.001
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.040
 LAB COORDINATES OF ROTATION AXIS -0.999995 -0.002850 -0.001317
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000720  0.005301  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.34   1244.07
 DETECTOR ORIGIN (PIXELS) AT                     1265.39   1237.10
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.29
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     506 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   51     33     65     69     54     31     58     70     59
 SIGMAB (degree)       0.032  0.033  0.033  0.033  0.032  0.032  0.031  0.032  0.032
 SIGMAR (degree)       0.033  0.028  0.030  0.033  0.031  0.027  0.031  0.032  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    506 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  0  0  0  0     0  0  0  0  3 10 14 11  6  2  1  0  0
     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  1  2  4  3  2  1  0  0  0     0  0  0  1  5 14 21 17  9  4  1  0  0
     0  0  0  0  0  1  1  0  1  0  0  0  0     0  0  0  0  0  2  3  2  1  1  0  0  0     0  0  0  0  3 10 13 11  5  2  1  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  2  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  1  4  9 11  7  3  1  0  0  0     0  0  0  2  7 14 15 10  4  2  0  0  0     0  0  0  1  5  9  9  5  2  1  0  0  0
     0  0  0  2 14 38 47 32 16  6  2  1  0     0  0  0  5 26 62 67 40 19  8  3  1  0     0  0  1  4 18 39 40 23 10  4  2  0  0
     0  0  0  4 22 54 69 48 23 10  4  1  0     0  0  0  8 41 92100 63 30 13  5  1  0     0  0  1  7 30 59 58 35 17  8  3  1  0
     0  0  0  2 12 33 39 29 13  6  2  0  0     0  0  0  4 23 55 58 40 18  8  2  1  0     0  0  0  4 17 36 37 24 11  5  1  0  0
     0  0  0  0  3  7  9  5  3  1  0  0  0     0  0  0  1  6 14 15  8  5  2  1  0  0     0  0  0  1  5 11 11  6  3  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  3  3  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5 10 10  6  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  8 15 14  8  4  2  1  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5 10 10  6  3  1  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  3  3  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  1  0  0  0  0  0  0  0  0  0  0

    647 REFLECTION INTENSITIES INTEGRATED


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
   121   0  1.066  5167682    0    132      80     7  0.04178  0.08677
   122   0  1.073  5167925    0    129      80     3  0.04169  0.11869
   123   0  1.072  5167990    0    130      81     1  0.04443  0.14097
   124   0  1.072  5167733    0    126      83     1  0.04364  0.15863
   125   0  1.077  5167742    0    122      89     4  0.04510  0.14891
   126   0  1.076  5167926    0    129      82     0  0.04817  0.21669
   127   0  1.078  5167867    0    124      87     7  0.04857  0.15545
   128   0  1.077  5166836    0    119      93     8  0.04535  0.14551
   129   0  1.075  5166655    0    120      91     5  0.04591  0.16166
   130   0  1.076  5166392    0    132      96     9  0.04499  0.12858
   131   0  1.075  5166418    0    139      98    10  0.04509  0.11875
   132   0  1.074  5166414    0    130      93     8  0.04317  0.11316
   133   0  1.078  5166431    0    132      95    10  0.04348  0.09951
   134   0  1.082  5166982    0    122      89     9  0.04375  0.08478
   135   0  1.081  5167906    0    109      80     4  0.04393  0.16040
   136   0  1.078  5169351    0    110      65     7  0.04992  0.12035
   137   0  1.082  5169037    0    115      70     7  0.04417  0.13236
   138   0  1.082  5168424    0    112      78     8  0.04164  0.10561
   139   0  1.082  5168676    0    120      73     8  0.04460  0.09862
   140   0  1.085  5168835    0    124      73     8  0.04371  0.11282

   438 OUT OF    438 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.48
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.426    32.121    34.497  88.637 108.431 111.893
 REC. CELL PARAMETERS   0.041647  0.033746  0.030731  83.863  70.654  67.343
 COORDINATES OF UNIT CELL A-AXIS   -14.333   -23.359    -1.058
 COORDINATES OF UNIT CELL B-AXIS    31.412    -5.393     3.989
 COORDINATES OF UNIT CELL C-AXIS    -1.189    12.072    32.294
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.017    -0.045     0.001
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.040
 LAB COORDINATES OF ROTATION AXIS -0.999995 -0.002850 -0.001317
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000688  0.005057  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.38   1244.06
 DETECTOR ORIGIN (PIXELS) AT                     1265.48   1237.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.26
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     438 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   43     29     51     62     51     22     49     58     56
 SIGMAB (degree)       0.033  0.032  0.033  0.033  0.032  0.033  0.031  0.032  0.032
 SIGMAR (degree)       0.032  0.028  0.031  0.032  0.031  0.028  0.030  0.032  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    438 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  3  8 11  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  2  2  1  1  0  0  0     0  0  0  1  5 15 21 16  9  4  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  1  0  0  0     0  0  0  0  3 10 15 11  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  3  9 11  6  3  1  0  0  0     0  0  0  1  6 14 16  9  4  1  1  0  0     0  0  0  1  5 10 10  5  2  1  0  0  0
     0  0  0  2 13 35 42 29 15  6  2  0  0     0  0  0  5 25 60 64 40 20  8  3  1  0     0  0  0  4 20 40 39 23 11  4  2  0  0
     1  0  0  3 22 56 71 47 24 10  4  1  0     0  0  0  7 41 91100 62 30 13  5  1  0     0  0  0  7 30 59 60 36 17  7  3  1  0
     1  0  0  2 13 34 43 29 13  6  2  1  0     0  0  0  4 24 54 61 38 16  7  3  0  0     0  0  0  4 18 37 38 23 10  5  2  0  0
     0  0  0  0  3  7  9  6  3  1  0  0  0     0  0  0  1  6 13 14  9  5  2  0  0  0     0  0  0  1  5 10 10  6  3  1  0  0  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5 10  9  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0 -1
     0  0  0  2  8 15 14  8  4  2  1  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  6 10 10  6  3  1  0  0  0     0  0  0  0  1  1  1  1  1  0  0  0  1     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  3  3  2  1  1  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    597 REFLECTION INTENSITIES INTEGRATED


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
   141   0  1.085  5168311    0    130      80     8  0.04153  0.11579
   142   0  1.082  5167801    0    133      83     7  0.04176  0.08109
   143   0  1.083  5168284    0    135      80     8  0.04265  0.07926
   144   0  1.085  5167755    0    140      83     3  0.04432  0.09407
   145   0  1.086  5166542    0    144      94     4  0.04156  0.08684
   146   0  1.083  5166943    0    153      91     3  0.04022  0.07709
   147   0  1.090  5167328    0    150      89     3  0.04072  0.08396
   148   0  1.088  5167912    0    147      83     4  0.03754  0.08212
   149   0  1.088  5168296    0    143      88     2  0.04150  0.12998
   150   0  1.087  5168696    0    132      76     7  0.04235  0.09785
   151   0  1.093  5169271    0    113      68     5  0.04197  0.10664
   152   0  1.095  5169624    0    123      64     3  0.04096  0.11137
   153   0  1.087  5168128    0    121      87     9  0.04005  0.09208
   154   0  1.088  5168104    0    128      87     8  0.04135  0.08429
   155   0  1.090  5168168    0    143      82     9  0.04008  0.08929
   156   0  1.087  5167823    0    138      81     7  0.03855  0.06538
   157   0  1.086  5168201    0    145      83     6  0.04221  0.09071
   158   0  1.087  5168288    0    145      79     7  0.04287  0.07820
   159   0  1.089  5167160    0    139      89     4  0.03825  0.06896
   160   0  1.089  5167935    0    150      81     6  0.04176  0.07332

   512 OUT OF    512 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.59
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.427    32.124    34.483  88.608 108.416 111.912
 REC. CELL PARAMETERS   0.041644  0.033744  0.030739  83.896  70.680  67.334
 COORDINATES OF UNIT CELL A-AXIS   -14.337   -23.358    -1.050
 COORDINATES OF UNIT CELL B-AXIS    31.416    -5.386     3.998
 COORDINATES OF UNIT CELL C-AXIS    -1.183    12.066    32.282
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.014    -0.038    -0.002
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.040
 LAB COORDINATES OF ROTATION AXIS -0.999995 -0.002850 -0.001317
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000588  0.005132  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.34   1243.99
 DETECTOR ORIGIN (PIXELS) AT                     1265.56   1237.24
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.29
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     512 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   48     30     63     69     53     29     57     76     63
 SIGMAB (degree)       0.033  0.034  0.033  0.033  0.032  0.032  0.031  0.032  0.032
 SIGMAR (degree)       0.033  0.029  0.030  0.033  0.032  0.027  0.031  0.032  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    512 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  2  4  5  3  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  3  1  1  0  0  0     0  0  0  1  5 15 22 16  8  3  1  0  0
     0  0  0  0  0  0  0  1  1  0  0  0  0     0  0  0  0  1  3  5  5  3  1  0  0  0     0  0  0  1  8 24 33 26 14  6  2  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  1  5 15 19 14  8  3  1  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  1  0     0  0  0  0  1  3  5  3  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0 -1  0     0  0  0  0  1  1  1  0  0  0  0 -1  0
     0  0  0  1  6 14 14  9  4  1  1  0  0     0  0  0  2 10 17 16  9  4  1  0  0  0     0  0  0  2  6  9  8  4  1  1  0  0  0
     0  0  0  3 18 49 61 40 19  7  3  1  0     0  0  0  6 29 63 69 41 20  7  3  1  0     0  0  1  4 16 31 32 18  9  3  1  0  0
     0  0  0  5 30 72 86 61 31 12  5  1  0     0  0  1 10 46 95100 64 32 12  5  1  0     0  0  1  6 27 49 48 29 14  6  2  0  0
     0  0  0  3 18 44 50 34 17  7  3  0  0     0  0  0  7 30 60 61 39 18  7  3  0  0     0  0  0  5 19 34 31 19  8  3  1  0  0
     0  0  0  1  4 10 11  8  3  2  0  0  0     0  0  0  2  7 16 15 10  4  2  0  0  0     0  0  0  1  5 10  9  6  2  1  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0 -1  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  6  6  3  1  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0
     0  0  0  2  6 10  9  5  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  4  7  6  4  2  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  1  0  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    652 REFLECTION INTENSITIES INTEGRATED


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
   161   0  1.087  5168382    0    139      75     3  0.04088  0.09409
   162   0  1.094  5167658    0    137      85     2  0.04111  0.09884
   163   0  1.096  5167749    0    144      85     4  0.04338  0.07862
   164   0  1.095  5168092    0    141      81     5  0.04318  0.07221
   165   0  1.092  5168092    0    140      79     5  0.04131  0.07869
   166   0  1.095  5167871    0    139      79     5  0.04192  0.06615
   167   0  1.096  5168066    0    141      83     5  0.03909  0.08226
   168   0  1.095  5167837    0    142      88     6  0.04120  0.08691
   169   0  1.094  5168378    0    138      74     5  0.04207  0.07084
   170   0  1.093  5167948    0    137      84     4  0.04224  0.08144
   171   0  1.095  5168355    0    136      81     4  0.04309  0.10297
   172   0  1.092  5169097    0    134      72     3  0.04392  0.08578
   173   0  1.097  5168490    0    137      77     2  0.04072  0.08065
   174   0  1.095  5167702    0    133      82     1  0.04170  0.09855
   175   0  1.095  5168348    0    142      77     4  0.03711  0.07666
   176   0  1.095  5168201    0    138      77     5  0.04066  0.06628
   177   0  1.098  5168398    0    132      80     8  0.04012  0.06447
   178   0  1.102  5169129    0    121      73     5  0.04536  0.09579
   179   0  1.098  5169408    0    118      70     4  0.04037  0.08944
   180   0  1.101  5168936    0    126      78     4  0.04316  0.11029

   477 OUT OF    477 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.53
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.428    32.121    34.502  88.609 108.440 111.913
 REC. CELL PARAMETERS   0.041650  0.033749  0.030728  83.884  70.653  67.331
 COORDINATES OF UNIT CELL A-AXIS   -14.336   -23.360    -1.052
 COORDINATES OF UNIT CELL B-AXIS    31.413    -5.383     3.999
 COORDINATES OF UNIT CELL C-AXIS    -1.184    12.086    32.294
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.007    -0.025    -0.002
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.040
 LAB COORDINATES OF ROTATION AXIS -0.999995 -0.002850 -0.001317
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000547  0.005207  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.32   1244.12
 DETECTOR ORIGIN (PIXELS) AT                     1265.60   1237.28
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.27
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     477 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     33     58     71     55     31     53     61     54
 SIGMAB (degree)       0.033  0.034  0.034  0.033  0.032  0.033  0.032  0.032  0.032
 SIGMAR (degree)       0.033  0.028  0.031  0.034  0.032  0.028  0.030  0.032  0.031


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    477 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  1  3  4  3  1  1  0  0  0
     0  0  0  0  0  1  0  0  1  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0     0  0  0  1  4 11 14 10  6  2  1  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0     0  0  0  1  5 16 21 16  9  4  2  1  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  2  2  2  1  0  0  0  0     0  0  0  1  3 10 13 10  5  2  1  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  1  3  4  3  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  2  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  4 10 12  7  4  2  0  0  0     0  0  0  1  8 17 16  9  5  2  0  0  0     0  0  0  1  6 11  9  5  2  1  0  0  0
     0  0  0  2 13 38 46 31 16  7  3  1  0     0  0  0  5 27 65 69 41 20  9  3  1  0     0  0  0  5 21 43 40 22 10  4  1  0  0
     0  0  0  4 21 57 66 48 25 11  4  1  0     0  0  0  8 42 96100 66 31 13  5  2  0     0  0  1  8 33 64 59 36 17  7  2  1  0
     0  0  0  2 12 33 41 30 15  6  2  1  0     0  0  0  5 25 58 63 41 20  8  3  1  0     0  0  0  5 21 39 38 23 11  5  1  0  0
     0  0  0  1  3  8 11  7  3  1  0  0  0     0  0  0  1  6 15 17 10  5  2  1  0  0     0  0  0  1  5 11 10  6  3  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  3  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 11 10  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  9 16 15  8  4  2  1  0  0     0  0  0  1  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 10  9  6  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  3  3  2  1  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
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
   181   0  1.098  5168897    0    130      70     1  0.04354  0.09984
   182   0  1.100  5168792    0    130      75     3  0.04275  0.09143
   183   0  1.105  5168812    0    128      71     4  0.04742  0.09906
   184   0  1.103  5168986    0    138      68     5  0.04415  0.08475
   185   0  1.101  5167990    0    142      78     7  0.04567  0.08663
   186   0  1.099  5167523    0    134      84     6  0.04046  0.07872
   187   0  1.093  5167997    0    142      83     9  0.04229  0.08030
   188   0  1.096  5167756    0    139      93     5  0.04322  0.09874
   189   0  1.096  5166883    0    152      94     4  0.04535  0.08653
   190   0  1.098  5166130    0    145      95     1  0.04119  0.09316
   191   0  1.098  5166971    0    161      89     5  0.04223  0.07874
   192   0  1.101  5166937    0    157      90     2  0.04149  0.07786
   193   0  1.103  5166558    0    150     100     2  0.04151  0.09976
   194   0  1.102  5167829    0    148      84     0  0.04529  0.12745
   195   0  1.104  5167215    0    146      95     6  0.04492  0.12466
   196   0  1.107  5166651    0    153      98     7  0.04436  0.08670
   197   0  1.106  5166039    0    146     102     3  0.04466  0.11325
   198   0  1.108  5166853    0    150      91     5  0.04305  0.09524
   199   0  1.110  5166740    0    147      95     6  0.04343  0.09645
   200   0  1.108  5166133    0    145     102     6  0.04254  0.08304

   493 OUT OF    493 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.51
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.425    32.120    34.492  88.645 108.430 111.891
 REC. CELL PARAMETERS   0.041648  0.033747  0.030737  83.856  70.652  67.342
 COORDINATES OF UNIT CELL A-AXIS   -14.335   -23.357    -1.049
 COORDINATES OF UNIT CELL B-AXIS    31.411    -5.396     3.985
 COORDINATES OF UNIT CELL C-AXIS    -1.187    12.081    32.285
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.017    -0.029     0.005
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.040
 LAB COORDINATES OF ROTATION AXIS -0.999995 -0.002850 -0.001317
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000843  0.005053  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.35   1244.04
 DETECTOR ORIGIN (PIXELS) AT                     1265.24   1237.40
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.24
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     493 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   53     34     59     64     54     29     54     70     61
 SIGMAB (degree)       0.032  0.033  0.033  0.033  0.032  0.032  0.031  0.031  0.032
 SIGMAR (degree)       0.031  0.029  0.031  0.033  0.032  0.028  0.030  0.031  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    493 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  3  3  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  3 10 12 10  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  2  1  1  0  0  0     0  0  0  1  4 15 21 17  9  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  1  3  9 13 10  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  3  8 10  7  3  1  0  0  0     0  0  0  1  6 12 15  9  4  1  0  0  0     0  0  0  1  5  8  9  5  2  1  0  0  0
     0  0  0  1 12 34 42 31 15  6  2  0  0     0  0  0  4 24 57 63 41 19  8  3  1  0     0  0  0  4 20 38 37 22 10  4  1  0  0
     0  0  0  3 18 53 70 50 24 10  4  1  0     0  0  1  7 38 89100 66 30 13  5  1  0     0  0  1  7 31 58 58 36 16  7  2  0  0
     0  0  0  2 12 33 42 29 13  6  2  0  0     0  0  0  4 24 56 62 40 17  8  3  0  0     0  0  0  4 19 37 37 21  9  4  1  0  0
     0  0  0  1  3  8 10  6  3  1  0  0  0     0  0  0  1  6 13 16 10  4  2  0  0  0     0  0  0  1  5  9 10  6  2  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0 -1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 10  8  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  8 14 12  7  3  1  0  0  0     0  0  0  1  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  9  8  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    660 REFLECTION INTENSITIES INTEGRATED


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
   201   0  1.110  5167223    0    155      93     4  0.04285  0.09111
   202   0  1.109  5167724    0    147      84     3  0.04614  0.10773
   203   0  1.108  5167341    0    141      88     3  0.04190  0.10578
   204   0  1.110  5168464    0    131      76     5  0.04252  0.08845
   205   0  1.111  5168713    0    130      72     5  0.04434  0.10308
   206   0  1.111  5168192    0    114      80     2  0.04316  0.11844
   207   0  1.112  5168417    0    120      77     0  0.04375  0.10784
   208   0  1.108  5167790    0    134      83     3  0.04194  0.07708
   209   0  1.106  5166522    0    148      98     6  0.04051  0.07128
   210   0  1.107  5166949    0    151      95     2  0.04410  0.09773
   211   0  1.112  5166714    0    151      94     3  0.04409  0.09800
   212   0  1.108  5167096    0    139      86     1  0.04700  0.12478
   213   0  1.109  5167690    0    140      86     1  0.04610  0.12195
   214   0  1.112  5167423    0    139      81     5  0.04508  0.08717
   215   0  1.106  5167889    0    134      80     3  0.04380  0.09546
   216   0  1.106  5167460    0    134      88     3  0.04546  0.09618
   217   0  1.107  5167249    0    127      80     8  0.04160  0.06104
   218   0  1.107  5168194    0    119      78     8  0.04324  0.07334
   219   0  1.108  5168975    0    133      69     6  0.04322  0.07029
   220   0  1.105  5167924    0    139      83     6  0.04270  0.07690

   468 OUT OF    468 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.52
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.11
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.419    32.117    34.494  88.677 108.419 111.876
 REC. CELL PARAMETERS   0.041654  0.033747  0.030734  83.829  70.652  67.348
 COORDINATES OF UNIT CELL A-AXIS   -14.334   -23.350    -1.046
 COORDINATES OF UNIT CELL B-AXIS    31.409    -5.406     3.971
 COORDINATES OF UNIT CELL C-AXIS    -1.188    12.081    32.288
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.027    -0.028     0.010
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.040
 LAB COORDINATES OF ROTATION AXIS -0.999995 -0.002850 -0.001317
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.001141  0.005090  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.37   1244.08
 DETECTOR ORIGIN (PIXELS) AT                     1264.87   1237.39
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.23
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     468 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   49     29     54     66     55     31     58     62     51
 SIGMAB (degree)       0.033  0.033  0.033  0.032  0.032  0.033  0.031  0.032  0.032
 SIGMAR (degree)       0.033  0.029  0.032  0.034  0.032  0.028  0.029  0.032  0.031


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    468 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  0  3  9 15 12  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  2  1  0  0  0     0  0  0  1  4 15 23 19 10  4  2  0  0
     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  0  3  9 14 10  6  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  3  8 11  7  3  1  0  0  0     0  0  0  1  7 14 15  9  3  1  1  0  0     0  0  0  1  6 11  9  6  2  1  0  0  0
     0  0  0  2 15 38 50 36 16  7  3  1  0     0  0  1  4 28 62 68 45 20  9  3  1  0     0  0  1  4 22 45 42 26 11  4  2  0  0
     0  0  0  3 21 58 75 53 27 11  5  1  0     0  0  1  8 43 92100 66 33 13  5  1  0     0  0  1  7 34 63 62 39 19  8  3  1  0
     0  0  0  2 13 36 44 29 15  6  2  1  0     0  0  0  6 25 58 63 40 19  8  3  1  0     0  0  0  5 21 41 41 26 12  5  2  0  0
     0  0  0  0  3  8  9  7  3  1  0  0  0     0  0  0  1  7 15 14 10  4  2  1  0  0     0  0  0  1  6 12 11  8  3  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  3  3  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 14 11  6  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3 11 20 17 10  5  2  1  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  2  7 13 12  8  3  2  1  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0
     0  0  0  1  2  5  4  3  1  1  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    649 REFLECTION INTENSITIES INTEGRATED


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
   221   0  1.108  5167547    0    142      88     6  0.04297  0.08000
   222   0  1.112  5167211    0    152      95     5  0.04519  0.09697
   223   0  1.112  5167682    0    140      91     5  0.04064  0.09609
   224   0  1.111  5167482    0    146      89     9  0.03963  0.06861
   225   0  1.110  5168031    0    134      83     7  0.04281  0.07029
   226   0  1.114  5168091    0    134      85     2  0.04194  0.11898
   227   0  1.115  5168500    0    135      74     1  0.04341  0.10273
   228   0  1.116  5168110    0    139      78     5  0.04130  0.08698
   229   0  1.113  5167395    0    139      87     7  0.04076  0.07112
   230   0  1.112  5167582    0    138      87     9  0.04337  0.06495
   231   0  1.116  5168610    0    126      74     5  0.04321  0.08517
   232   0  1.117  5168526    0    129      74     5  0.04329  0.07682
   233   0  1.118  5168017    0    127      80     3  0.04073  0.08362
   234   0  1.117  5167913    0    129      81     6  0.03692  0.06077
   235   0  1.117  5168042    0    136      79     2  0.04116  0.08079
   236   0  1.122  5168165    0    129      77     1  0.04276  0.09184
   237   0  1.120  5168364    0    136      76     1  0.03989  0.08988
   238   0  1.118  5168771    0    129      70     0  0.03973  0.09669
   239   0  1.121  5169948    0    131      59     6  0.04269  0.07689
   240   0  1.123  5168928    0    143      75     2  0.04372  0.11033

   467 OUT OF    467 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.55
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.434    32.123    34.505  88.623 108.449 111.907
 REC. CELL PARAMETERS   0.041644  0.033746  0.030727  83.866  70.639  67.331
 COORDINATES OF UNIT CELL A-AXIS   -14.336   -23.366    -1.057
 COORDINATES OF UNIT CELL B-AXIS    31.415    -5.384     3.998
 COORDINATES OF UNIT CELL C-AXIS    -1.190    12.091    32.296
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.002    -0.030    -0.002
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.040
 LAB COORDINATES OF ROTATION AXIS -0.999995 -0.002850 -0.001317
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000520  0.005132  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.35   1244.07
 DETECTOR ORIGIN (PIXELS) AT                     1265.66   1237.33
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.31
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     467 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   45     30     53     64     55     28     53     66     59
 SIGMAB (degree)       0.033  0.033  0.033  0.033  0.033  0.033  0.031  0.032  0.032
 SIGMAR (degree)       0.031  0.029  0.032  0.034  0.032  0.028  0.029  0.031  0.030


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    467 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  3  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  0  2 10 13 11  6  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  2  1  0  0  0     0  0  0  1  4 14 20 17 10  4  2  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  0  3  8 12 10  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  3  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  1  1  1  0  0  0 -1  0     0  0  0  0  1  2  1  1  0  0  0 -1  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  3  8 10  8  3  1  0  0  0     0  0  0  1  6 14 14 10  4  2  0  0  0     0  0  0  2  5 10  7  5  2  1  0  0  0
     0  0  0  2 11 36 45 32 16  7  2  1  0     0  0  0  5 25 62 67 41 19  8  3  1  0     0  0  0  5 21 42 39 21  9  4  1  0  0
     0  0  0  3 19 52 68 49 26 12  4  1  0     0  0  1  7 39 91100 65 33 14  5  2  0     0  0  0  7 31 61 58 35 17  7  3  1  0
     0  0  0  2 13 31 41 30 15  7  2  1  0     0  0  0  4 25 54 61 41 20  8  3  1  0     0  0  0  4 20 38 37 22 10  4  1  0  0
     0  0  0  1  3 10 10  6  3  1  0  0  0     0  0  0  1  7 17 17 10  4  2  0  0  0     0  0  0  2  6 11 11  6  2  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0 -1 -1  0 -1  0  0  0     0  0  0  0  0  0 -1 -1  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 10  8  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  2  8 14 13  7  4  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  9  9  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  3  3  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    641 REFLECTION INTENSITIES INTEGRATED


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
   241   0  1.121  5167828    0    128      86     4  0.04105  0.09707
   242   0  1.123  5168062    0    143      77     4  0.04133  0.07141
   243   0  1.124  5167591    0    149      82     2  0.03977  0.07157
   244   0  1.120  5167233    0    137      89     1  0.03902  0.07805
   245   0  1.122  5168139    0    139      79     0  0.04494  0.08559
   246   0  1.122  5168667    0    131      74     2  0.04060  0.06929
   247   0  1.119  5168370    0    132      82     2  0.03858  0.07945
   248   0  1.124  5167694    0    143      82     1  0.03851  0.07451
   249   0  1.124  5167457    0    143      91     4  0.03839  0.07022
   250   0  1.126  5168609    0    146      76     3  0.04249  0.07295
   251   0  1.121  5168339    0    140      74     2  0.03973  0.07064
   252   0  1.121  5167807    0    137      79     4  0.04271  0.07420
   253   0  1.121  5167777    0    146      82     2  0.04023  0.07578
   254   0  1.126  5167389    0    144      88     0  0.04284  0.09222
   255   0  1.126  5167724    0    152      89     2  0.04057  0.10031
   256   0  1.128  5167340    0    149      91     3  0.04350  0.09548
   257   0  1.127  5167547    0    135      84     5  0.04222  0.08064
   258   0  1.126  5167902    0    140      81     2  0.04537  0.09940
   259   0  1.134  5168104    0    139      80     4  0.04113  0.07661
   260   0  1.132  5167774    0    135      83     0  0.04181  0.09101

   485 OUT OF    487 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.53
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.03
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.428    32.120    34.495  88.645 108.440 111.898
 REC. CELL PARAMETERS   0.041649  0.033748  0.030735  83.849  70.640  67.334
 COORDINATES OF UNIT CELL A-AXIS   -14.336   -23.360    -1.048
 COORDINATES OF UNIT CELL B-AXIS    31.413    -5.392     3.985
 COORDINATES OF UNIT CELL C-AXIS    -1.186    12.091    32.285
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.014    -0.021     0.003
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.040
 LAB COORDINATES OF ROTATION AXIS -0.999995 -0.002850 -0.001317
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000744  0.005092  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.35   1244.03
 DETECTOR ORIGIN (PIXELS) AT                     1265.37   1237.34
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.28
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     485 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   50     29     54     71     56     32     59     68     55
 SIGMAB (degree)       0.032  0.033  0.033  0.033  0.033  0.033  0.031  0.032  0.032
 SIGMAR (degree)       0.031  0.028  0.031  0.033  0.032  0.027  0.030  0.032  0.031


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    485 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  4  3  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  1  0  0  0     0  0  0  0  3 10 15 11  5  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  2  1  0  0  0     0  0  0  1  5 15 23 18 10  4  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  1  0  0  0     0  0  0  0  3  9 13 11  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  3  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0 -1 -1 -1  0  0  0     0  0  0  0  0  0  0  0 -1 -1  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  1  1  1  0  0 -1 -1  0     0  0  0  0  1  1  1  1  0  0  0 -1  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  3  8 10  7  4  1  0  0  0     0  0  0  1  6 13 15  9  4  2  0  0  0     0  0  0  1  5  8  7  4  2  1  0  0  0
     0  0  0  2 13 36 46 31 14  7  3  0  0     0  0  0  5 27 64 64 39 17  8  3  1  0     0  0  0  5 20 40 34 19  8  3  1  0  0
     0  0  0  3 18 53 71 51 26 11  4  1  0     0  0  0  7 40 92100 64 31 13  5  1  0     0  0  0  7 30 58 53 31 14  6  2  0  0
     0  0  0  2 12 34 42 30 16  7  2  1  0     0  0  0  5 24 58 63 39 19  8  3  1  0     0  0  0  5 19 38 35 20  9  4  2  1  0
     0  0  0  1  4 10 11  7  3  1  0  0  0     0  0  0  1  7 16 17 10  4  2  0  0  0     0  0  0  1  6 10 10  5  3  1  0  0  0
     0  0  0  0  1  1  1  1  0  0 -1 -1  0     0  0  0  0  1  2  2  1  0  0 -1 -1  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0 -1  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0 -1 -1  0  0  0  0     0  0  0  0  0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  9  7  3  1  1  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  8 14 12  6  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5 10  8  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  3  3  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

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
   261   0  1.129  5167555    0    131      83     0  0.04265  0.09199
   262   0  1.130  5168791    0    134      77     1  0.04286  0.11179
   263   0  1.129  5168808    0    136      76     0  0.03881  0.12403
   264   0  1.130  5168040    0    130      81     1  0.04477  0.11645
   265   0  1.132  5167803    0    146      85     3  0.04509  0.10757
   266   0  1.132  5167456    0    142      86     2  0.04645  0.11337
   267   0  1.132  5167680    0    131      81     2  0.04589  0.10829
   268   0  1.134  5168039    0    138      75     1  0.04575  0.10164
   269   0  1.135  5168573    0    133      68     0  0.04478  0.09618
   270   0  1.134  5168554    0    130      70     0  0.04288  0.10085
   271   0  1.136  5168025    0    126      74     2  0.04324  0.08304
   272   0  1.134  5167892    0    126      76     1  0.04260  0.08941
   273   0  1.134  5168506    0    127      70     0  0.04603  0.10438
   274   0  1.133  5168294    0    143      77     0  0.04378  0.11030
   275   0  1.133  5167564    0    139      85     1  0.03987  0.08291
   276   0  1.131  5167451    0    140      85     4  0.04189  0.08036
   277   0  1.131  5168295    0    152      74     3  0.04690  0.09078
   278   0  1.132  5168278    0    141      74     2  0.04350  0.07699
   279   0  1.133  5168164    0    129      79     3  0.03729  0.06615
   280   0  1.137  5168769    0    124      72     2  0.04094  0.07521

   457 OUT OF    460 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.49
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.04
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.432    32.122    34.503  88.636 108.440 111.900
 REC. CELL PARAMETERS   0.041642  0.033747  0.030728  83.858  70.644  67.334
 COORDINATES OF UNIT CELL A-AXIS   -14.337   -23.364    -1.048
 COORDINATES OF UNIT CELL B-AXIS    31.414    -5.388     3.991
 COORDINATES OF UNIT CELL C-AXIS    -1.188    12.095    32.292
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.007    -0.019     0.001
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.040
 LAB COORDINATES OF ROTATION AXIS -0.999995 -0.002850 -0.001317
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000585  0.005131  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.37   1244.08
 DETECTOR ORIGIN (PIXELS) AT                     1265.60   1237.34
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.34
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     457 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   43     32     61     67     49     23     50     66     58
 SIGMAB (degree)       0.033  0.033  0.033  0.033  0.033  0.033  0.031  0.032  0.032
 SIGMAR (degree)       0.031  0.028  0.031  0.034  0.032  0.028  0.030  0.031  0.031


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    457 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  3  4  3  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  3  2  1  1  0  0  0     0  0  0  0  3 10 16 13  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  4  3  2  1  0  0  0     0  0  0  1  5 15 23 20 11  5  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  1  0  0  0     0  0  0  1  4 11 14 12  7  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  4  3  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0 -1  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  3  9 12  8  4  1  0  0  0     0  0  0  1  7 16 17  9  4  1  0  0  0     0  0  0  1  5 10  9  5  2  1  0  0  0
     0  0  0  2 13 36 51 36 18  8  3  0  0     0  0  0  5 28 62 70 43 20  9  3  0  0     0  0  0  4 21 39 37 22 10  4  1  0  0
     0  0  0  3 21 57 74 55 28 12  5  1  0     0  0  1  8 42 95100 66 32 14  5  1  0     0  0  1  7 31 59 55 33 16  7  2  0  0
     0  0  0  3 15 40 46 33 17  7  3  1  0     0  0  0  5 28 63 64 41 20  8  3  1  0     0  0  0  4 20 38 35 21 10  4  1  0  0
     0  0  0  1  3  9 13  8  4  1  0  0  0     0  0  0  1  6 15 18 10  5  2  0  0  0     0  0  0  1  5  9 10  5  2  1  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0 -1  0 -1  0  0  0     0  0  0  0  0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  8  7  4  2  1  0  0  0     1  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  8 12 10  6  3  1  0  0  0     0  0  0  1  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  9  7  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0 -1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

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
   281   0  1.136  5169256    0    123      66     1  0.04354  0.08808
   282   0  1.137  5169483    0    122      61     1  0.03941  0.08001
   283   0  1.138  5168981    0    130      71     0  0.04363  0.09734
   284   0  1.137  5168196    0    145      82     0  0.04083  0.09437
   285   0  1.137  5167934    0    151      85     3  0.04077  0.07244
   286   0  1.137  5168212    0    142      80     4  0.03762  0.07383
   287   0  1.139  5169048    0    133      75     4  0.03975  0.07675
   288   0  1.141  5170308    0    130      58     2  0.04065  0.08012
   289   0  1.139  5169723    0    134      65     1  0.04077  0.08591
   290   0  1.134  5168817    0    140      76     1  0.03838  0.07336
   291   0  1.140  5167980    0    144      84     3  0.03891  0.06765
   292   0  1.141  5167800    0    140      86     1  0.04105  0.09061
   293   0  1.140  5168894    0    134      70     1  0.04216  0.08625
   294   0  1.141  5169089    0    141      70     2  0.04391  0.08019
   295   0  1.143  5169035    0    137      71     0  0.04502  0.10566
   296   0  1.140  5167541    0    128      84     2  0.04147  0.08030
   297   0  1.141  5167869    0    130      81     3  0.03980  0.06757
   298   0  1.140  5168489    0    137      77     2  0.04305  0.08530
   299   0  1.138  5168602    0    127      73     1  0.04309  0.09827
   300   0  1.136  5168538    0    125      73     1  0.04089  0.08493

   467 OUT OF    469 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.50
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.430    32.119    34.502  88.651 108.456 111.885
 REC. CELL PARAMETERS   0.041647  0.033748  0.030733  83.839  70.622  67.343
 COORDINATES OF UNIT CELL A-AXIS   -14.333   -23.364    -1.054
 COORDINATES OF UNIT CELL B-AXIS    31.411    -5.393     3.984
 COORDINATES OF UNIT CELL C-AXIS    -1.188    12.096    32.290
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.014    -0.025     0.001
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.040
 LAB COORDINATES OF ROTATION AXIS -0.999995 -0.002850 -0.001317
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000805  0.005037  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.34   1243.98
 DETECTOR ORIGIN (PIXELS) AT                     1265.29   1237.37
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.30
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     467 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   47     31     57     64     54     29     54     66     59
 SIGMAB (degree)       0.033  0.034  0.033  0.033  0.033  0.033  0.031  0.032  0.032
 SIGMAR (degree)       0.033  0.029  0.032  0.033  0.031  0.029  0.030  0.031  0.031


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    467 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  2  2  1  1  0  0  0     0  0  0  0  2  9 13 11  6  3  1  0  0
     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  1  2  3  3  2  1  0  0  0     0  0  0  1  4 12 19 17  9  4  2  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  2  2  2  1  0  0  0  0     0  0  0  0  3  9 13 10  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  0  1  1  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  0  0  1  0  0  0  0  0  0  0     0  0  0  1  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  3  8 11  7  3  1  0  0  0     0  0  0  1  6 14 16  9  4  1  0  0  0     0  0  0  2  6 10  9  5  2  1  0  0  0
     0  0  0  2 11 35 46 34 17  7  2  1  0     0  0  1  5 26 60 67 42 20  8  3  1  0     0  0  1  6 22 41 41 23 11  4  1  0  0
     0  0  0  3 19 52 68 53 26 12  4  1  0     0  0  1  9 42 94100 67 32 14  5  1  0     0  0  1  8 35 64 59 36 17  7  3  0  0
     0  0  0  2 14 35 46 31 15  7  3  1  0     0  0  0  6 29 63 66 41 19  8  3  1  0     0  0  0  5 22 43 38 23 10  4  1  0  0
     0  0  0  0  3  9 10  8  4  1  0  0  0     0  0  0  1  8 17 15 10  5  2  0  0  0     0  0  0  1  5 11  9  6  3  1  0  0  0
     0  0  0  0  0  1  1  1  0  0 -1 -1  0     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  1  3  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 11  9  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3 10 16 14  8  4  2  1  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 11 10  5  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  3  2  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    616 REFLECTION INTENSITIES INTEGRATED


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
   301   0  1.142  5169133    0    131      69     1  0.04498  0.09665
   302   0  1.142  5168410    0    132      78     2  0.04332  0.08592
   303   0  1.141  5167435    0    132      86     2  0.03923  0.06639
   304   0  1.139  5167422    0    146      85     2  0.04072  0.07309
   305   0  1.139  5167884    0    147      79     2  0.04115  0.07421
   306   0  1.141  5167683    0    140      82     2  0.03996  0.07007
   307   0  1.141  5168071    0    155      80     0  0.04258  0.09009
   308   0  1.142  5167812    0    155      86     3  0.04115  0.07917
   309   0  1.145  5166784    0    152      93     0  0.03914  0.06906
   310   0  1.144  5167292    0    153      88     1  0.03744  0.06717
   311   0  1.147  5167157    0    166      88     3  0.03892  0.06828
   312   0  1.147  5166773    0    159     100     3  0.03979  0.08382
   313   0  1.145  5166746    0    160      93     5  0.04158  0.07080
   314   0  1.145  5167214    0    163      87     5  0.04044  0.06631
   315   0  1.145  5166629    0    161      96     1  0.03851  0.07741
   316   0  1.147  5166975    0    154      95     1  0.04183  0.08722
   317   0  1.146  5167260    0    153      90     2  0.04541  0.07881
   318   0  1.145  5167867    0    145      78     4  0.03980  0.06726
   319   0  1.147  5168074    0    143      79     3  0.04141  0.07726
   320   0  1.146  5168821    0    144      76     2  0.04163  0.09054

   501 OUT OF    503 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.55
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.04
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.431    32.120    34.490  88.638 108.454 111.897
 REC. CELL PARAMETERS   0.041648  0.033748  0.030742  83.851  70.628  67.336
 COORDINATES OF UNIT CELL A-AXIS   -14.336   -23.363    -1.053
 COORDINATES OF UNIT CELL B-AXIS    31.412    -5.391     3.988
 COORDINATES OF UNIT CELL C-AXIS    -1.185    12.091    32.280
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.014    -0.025     0.002
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.040
 LAB COORDINATES OF ROTATION AXIS -0.999995 -0.002850 -0.001317
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000621  0.005056  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.37   1243.96
 DETECTOR ORIGIN (PIXELS) AT                     1265.56   1237.32
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.29
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     501 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   52     37     60     69     56     33     58     65     59
 SIGMAB (degree)       0.033  0.033  0.033  0.033  0.033  0.033  0.031  0.032  0.032
 SIGMAR (degree)       0.030  0.028  0.030  0.032  0.031  0.026  0.029  0.032  0.031


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    501 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  1  2  4  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  0  3 10 14 11  6  3  1  0  0
     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  0  2  3  2  2  1  0  0  0     0  0  0  1  5 14 21 17  9  4  2  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  1  0  0  0     0  0  0  0  3 10 14 11  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  3  3  3  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  3  9 11  7  4  1  0  0  0     0  0  0  2  7 15 16  9  5  2  0  0  0     0  0  0  1  5  8  8  4  2  1  0  0  0
     0  0  0  2 12 36 46 33 16  7  3  1  0     0  0  0  5 26 62 66 42 19  8  3  1  0     0  0  0  5 19 36 33 19  8  3  1  0  0
     0  0  0  3 20 54 71 51 25 11  4  1  0     0  0  0  8 43 92100 64 28 13  5  1  0     0  0  1  7 31 55 51 30 13  6  2  0  0
     0  0  0  2 13 36 43 30 15  7  2  1  0     0  0  0  5 27 60 61 38 18  8  3  1 -1     0  0  0  4 19 36 32 18  9  4  1  0 -1
     0  0  0  1  4 10 11  8  4  2  0  0  0     0  0  0  1  7 17 17 10  5  2  0  0  0     0  0  0  1  6 10  9  5  2  1  0  0  0
     0  0  0  0  1  1  2  1  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  1  0  0  0  0  0 -1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  5  7  6  3  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  8 12 10  5  2  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  5  8  7  3  2  1  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  1     1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
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
   321   0  1.149  5168686    0    135      84     3  0.04188  0.08878
   322   0  1.148  5169172    0    126      71     3  0.04001  0.08199
   323   0  1.147  5169473    0    139      65     2  0.04189  0.09177
   324   0  1.149  5168448    0    130      79     4  0.04020  0.07339
   325   0  1.146  5168101    0    137      84     2  0.03869  0.07687
   326   0  1.149  5168302    0    130      81     4  0.03892  0.07592
   327   0  1.149  5168538    0    121      76     3  0.04009  0.07188
   328   0  1.149  5168741    0    130      74     4  0.04003  0.07969
   329   0  1.148  5168152    0    137      81     4  0.04407  0.08085
   330   0  1.148  5167229    0    147      92     7  0.04083  0.07194
   331   0  1.145  5167678    0    151      79     8  0.04034  0.06788
   332   0  1.144  5167811    0    150      81     2  0.04287  0.09139
   333   0  1.146  5168307    0    137      76     3  0.04179  0.07271
   334   0  1.147  5168509    0    145      73     2  0.03969  0.07469
   335   0  1.149  5168042    0    140      79     1  0.03873  0.07595
   336   0  1.145  5167937    0    145      84     3  0.04318  0.07794
   337   0  1.146  5168318    0    147      78     4  0.04130  0.08255
   338   0  1.146  5168729    0    144      75     3  0.04300  0.09411
   339   0  1.146  5168725    0    144      76     7  0.03953  0.06966
   340   0  1.145  5168763    0    135      72     5  0.03986  0.06434

   466 OUT OF    467 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.57
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.419    32.115    34.507  88.689 108.422 111.876
 REC. CELL PARAMETERS   0.041656  0.033751  0.030724  83.814  70.643  67.343
 COORDINATES OF UNIT CELL A-AXIS   -14.332   -23.352    -1.035
 COORDINATES OF UNIT CELL B-AXIS    31.408    -5.402     3.964
 COORDINATES OF UNIT CELL C-AXIS    -1.187    12.102    32.294
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.025    -0.002     0.006
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.040
 LAB COORDINATES OF ROTATION AXIS -0.999995 -0.002850 -0.001317
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.001192  0.005211  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.41   1244.11
 DETECTOR ORIGIN (PIXELS) AT                     1264.84   1237.27
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.27
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     466 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   49     27     52     65     51     28     54     71     58
 SIGMAB (degree)       0.033  0.033  0.033  0.033  0.033  0.034  0.031  0.032  0.032
 SIGMAR (degree)       0.033  0.031  0.032  0.033  0.031  0.026  0.030  0.032  0.033


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    466 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  3  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  0  2  8 13 10  6  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  3  3  2  1  0  0  0     0  0  0  0  4 12 20 17 10  4  2  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  2  1  0  0  0  0     0  0  0  0  3  8 12 11  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  3  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0


     0  0  0  0  0  0 -1  0  0 -1  0  0  0     0  0  0  0  0  0 -1  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0 -1  0     0  0  0  0  1  1  1  0  0  0  0 -1  0     0  0  0  0  0  1  1  0  0  0  0 -1  0
     0  0  0  1  3  8 10  7  3  1  0  0  0     0  0  0  1  6 13 14  9  4  2  0  0  0     0  0  0  1  5 10  9  5  2  1  0  0  0
     0  0  0  2 12 33 44 31 17  7  3  0  0     0  0  0  6 26 59 64 41 20  9  3  1  0     0  0  1  5 23 42 41 25 11  4  2  0  0
     0  0  0  3 20 51 68 53 27 11  5  1  0     0  0  1  9 42 92100 66 33 13  5  1  0     0  0  1  8 35 64 62 37 18  7  3  1  0
     0  0  0  2 13 32 44 33 14  6  2  0  0     0  0  0  6 28 58 64 42 18  8  3  0  0     0  0  0  5 21 40 39 25 11  5  2  0  0
     0  0  0  1  3  9 10  7  3  1  0  0  0     0  0  0  2  7 16 15 10  4  2  0  0  0     0  0  0  1  5 12 11  7  3  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0 -1  0     0  0  0  0  1  3  2  1  0  0  0 -1  0     0  0  0  0  1  2  2  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  4  3  2  1  0  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0
     0  0  0  2  9 14 13  7  3  1  0  0  0     0  0  0  1  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3 13 20 19 10  5  2  0  0  0     0  0  0  1  2  3  2  1  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  2  8 13 11  7  3  1  0  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  4  3  2  1  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0 -1  0  0  0  0  0  0  0  0  0

    632 REFLECTION INTENSITIES INTEGRATED


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
   341   0  1.147  5168633    0    134      78     6  0.03836  0.07663
   342   0  1.145  5168701    0    131      77     6  0.03990  0.07570
   343   0  1.149  5168593    0    138      77     5  0.04463  0.07858
   344   0  1.148  5167847    0    157      86     2  0.04309  0.10502
   345   0  1.147  5166993    0    162      91     4  0.04087  0.09177
   346   0  1.147  5166716    0    159      94     6  0.04037  0.07190
   347   0  1.145  5167777    0    151      77     5  0.04221  0.06365
   348   0  1.147  5168275    0    141      79     4  0.04060  0.08773
   349   0  1.145  5168922    0    126      73     2  0.04239  0.11510
   350   0  1.143  5169300    0    125      68     1  0.04551  0.13542
   351   0  1.147  5169070    0    126      71     3  0.04531  0.11240
   352   0  1.146  5168324    0    118      84     4  0.04263  0.10850
   353   0  1.145  5168230    0    133      78     2  0.04328  0.09672
   354   0  1.144  5167484    0    135      83     9  0.03869  0.05991
   355   0  1.142  5167371    0    145      83     3  0.03906  0.08659
   356   0  1.142  5168122    0    159      77     4  0.04341  0.10494
   357   0  1.147  5167321    0    142      91     2  0.04195  0.10249
   358   0  1.146  5167619    0    141      85     3  0.04334  0.09686
   359   0  1.145  5169180    0    150      68     5  0.04448  0.09767
   360   0  1.148  5168640    0    140      70     5  0.04194  0.06765

   492 OUT OF    495 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.57
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.424    32.118    34.494  88.663 108.435 111.884
 REC. CELL PARAMETERS   0.041651  0.033749  0.030736  83.835  70.640  67.342
 COORDINATES OF UNIT CELL A-AXIS   -14.332   -23.358    -1.043
 COORDINATES OF UNIT CELL B-AXIS    31.410    -5.395     3.979
 COORDINATES OF UNIT CELL C-AXIS    -1.189    12.096    32.282
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.014    -0.014     0.002
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.040
 LAB COORDINATES OF ROTATION AXIS -0.999995 -0.002850 -0.001317
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000765  0.005089  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.39   1244.07
 DETECTOR ORIGIN (PIXELS) AT                     1265.39   1237.39
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.28
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     492 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   50     32     61     71     54     28     58     67     56
 SIGMAB (degree)       0.033  0.033  0.033  0.032  0.033  0.033  0.031  0.032  0.032
 SIGMAR (degree)       0.031  0.029  0.031  0.032  0.031  0.026  0.028  0.031  0.031


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    492 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  3  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  0  3  9 15 12  7  3  1  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  1  3  3  2  1  0  0  0     0  0  0  0  4 13 21 17 11  5  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  1  1  0  0  0  0     0  0  0  0  2  8 12  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  1  0  0  0 -1  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  1  4  8 10  7  3  1  0 -1  0     0  0  0  1  7 14 13  9  3  1  0  0  0     0  0  0  1  6  9  8  4  2  1  0  0  0
     0  0  0  2 14 37 50 35 18  7  3  1  0     0  0  0  6 29 61 66 43 20  8  3  1  0     0  0  0  5 22 39 36 22  9  4  1  0  0
     0  0  0  3 21 55 72 52 28 12  4  1  0     0  0  1  9 42 91100 64 31 13  5  1  0     0  0  1  8 31 56 55 33 15  6  2  1  0
     0  0  0  2 13 34 43 29 16  6  2  1  0     0  0  0  5 27 58 61 37 19  7  2  1  0     0  0  1  4 19 34 34 19  9  4  1  0  0
     0  0  0  0  3  8 11  8  4  1  0  0  0     0  0  0  1  6 14 16 11  5  2  0  0  0     0  0  0  1  5  9  9  6  2  1  0  0  0
     0  0  0  0  0  1  2  1  0  0 -1 -1  0     0  0  0  0  1  2  3  1  0  0 -1  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  1  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6  9  7  4  2  1  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  8 12 11  6  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  1  5  8  7  3  2  1  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0
     0  0  0  1  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    659 REFLECTION INTENSITIES INTEGRATED


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
   361   0  1.148  5168188    0    134      72     1  0.03948  0.07492
   362   0  1.148  5168833    0    126      69     4  0.04137  0.06867
   363   0  1.145  5169699    0    115      60     6  0.03984  0.06429
   364   0  1.145  5169716    0    119      61     7  0.04103  0.07333
   365   0  1.145  5169069    0    118      72     7  0.04237  0.09006
   366   0  1.147  5168425    0    128      74     4  0.04336  0.08292
   367   0  1.145  5167262    0    137      89     9  0.03951  0.06416
   368   0  1.144  5167493    0    139      84     9  0.03817  0.05913
   369   0  1.140  5168876    0    145      68     4  0.04220  0.10143
   370   0  1.139  5168197    0    138      77     6  0.04012  0.08541
   371   0  1.138  5167669    0    136      81     8  0.03851  0.06374
   372   0  1.138  5167160    0    160      88     8  0.04088  0.07095
   373   0  1.138  5167259    0    157      85     6  0.04195  0.06922
   374   0  1.139  5167118    0    147      92     6  0.04159  0.09115
   375   0  1.136  5167562    0    137      86     8  0.03946  0.08818
   376   0  1.134  5168861    0    137      69     7  0.04486  0.09685
   377   0  1.134  5168867    0    132      75     7  0.04694  0.09186
   378   0  1.133  5167611    0    129      81     8  0.04155  0.07298
   379   0  1.136  5168078    0    130      82     8  0.03914  0.09626
   380   0  1.136  5169058    0    132      74     6  0.04211  0.11576

   463 OUT OF    467 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.51
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.11
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.422    32.111    34.487  88.623 108.461 111.900
 REC. CELL PARAMETERS   0.041663  0.033758  0.030746  83.863  70.626  67.337
 COORDINATES OF UNIT CELL A-AXIS   -14.330   -23.356    -1.052
 COORDINATES OF UNIT CELL B-AXIS    31.403    -5.385     3.996
 COORDINATES OF UNIT CELL C-AXIS    -1.185    12.095    32.275
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.004    -0.017    -0.001
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.040
 LAB COORDINATES OF ROTATION AXIS -0.999995 -0.002850 -0.001317
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000531  0.005058  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.32   1244.08
 DETECTOR ORIGIN (PIXELS) AT                     1265.62   1237.44
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.18
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     463 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   48     31     55     66     52     27     48     66     56
 SIGMAB (degree)       0.033  0.033  0.033  0.033  0.033  0.034  0.032  0.032  0.032
 SIGMAR (degree)       0.030  0.029  0.031  0.031  0.030  0.027  0.029  0.031  0.031


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    463 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  2  7 10  8  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  0  0  0  0     0  0  0  0  3 10 17 14  8  4  2  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  2  7 10  9  4  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  2  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  1  1  1  0  0  0 -1  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  3  8 10  6  3  1  0  0  0     0  0  0  1  6 16 16  9  4  2  0  0  0     0  0  0  2  6 12 10  5  2  1  0  0  0
     0  0  0  3 12 31 38 27 14  6  2  1  0     0  0  0  7 27 61 65 39 19  8  3  1  0     0  0  1  7 23 45 42 23 11  5  2  0  0
     0  0  0  3 18 46 62 44 24 10  4  1  0     0  0  1  8 40 90100 65 32 14  6  2  0     0  0  1 10 37 67 64 39 18  7  3  1  0
     0  0  0  2 12 30 38 28 14  6  2  1  0     0  0  0  5 25 57 64 40 19  8  3  1  0     0  0  1  5 23 41 41 23 11  5  2  0  0
     0  0  0  1  3  8  9  7  3  1  0  0  0     0  0  0  1  6 14 15  9  5  2  1  0  0     0  0  0  1  5 10  9  5  3  1  0  0  0
     0  0  0  0  1  1  2  1  1  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  3  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  7 11  9  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3 10 16 14  8  4  2  1  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0  0  0  0  0
     0  0  0  2  6 11 10  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  1  0  0  0  0  0  0  0  0  0  0  0

    631 REFLECTION INTENSITIES INTEGRATED


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
   381   0  1.133  5169127    0    127      73     7  0.04499  0.11335
   382   0  1.137  5168366    0    129      72     6  0.04159  0.08929
   383   0  1.139  5168223    0    139      73     6  0.03892  0.06413
   384   0  1.137  5168166    0    153      74     8  0.04688  0.07679
   385   0  1.134  5166927    0    140      94     9  0.03955  0.07445
   386   0  1.137  5167002    0    140      95     6  0.04028  0.09284
   387   0  1.136  5167333    0    146      89     3  0.04357  0.12103
   388   0  1.137  5167291    0    133      93     2  0.04311  0.12533
   389   0  1.134  5167808    0    128      81     5  0.04243  0.10308
   390   0  1.134  5167153    0    130      90     2  0.04409  0.13914
   391   0  1.137  5167884    0    133      84     3  0.04293  0.13184
   392   0  1.132  5167666    0    147      87     6  0.04346  0.09937
   393   0  1.133  5167158    0    139      88     9  0.04258  0.07899
   394   0  1.132  5167028    0    134      94    10  0.04214  0.08991
   395   0  1.136  5168274    0    136      76     8  0.04439  0.10609
   396   0  1.138  5168458    0    136      81     9  0.04354  0.09997
   397   0  1.132  5168552    0    134      74     6  0.04032  0.08532
   398   0  1.133  5168231    0    129      82     9  0.03797  0.07174
   399   0  1.135  5168476    0    132      79     6  0.04304  0.09411
   400   0  1.133  5167915    0    127      82     7  0.04109  0.09172

   466 OUT OF    468 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.54
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.15
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.427    32.115    34.505  88.639 108.460 111.890
 REC. CELL PARAMETERS   0.041653  0.033752  0.030730  83.850  70.622  67.342
 COORDINATES OF UNIT CELL A-AXIS   -14.331   -23.361    -1.053
 COORDINATES OF UNIT CELL B-AXIS    31.407    -5.390     3.992
 COORDINATES OF UNIT CELL C-AXIS    -1.189    12.101    32.291
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.005    -0.018     0.000
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.040
 LAB COORDINATES OF ROTATION AXIS -0.999995 -0.002850 -0.001317
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000740  0.005044  1.023529
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.29   1244.06
 DETECTOR ORIGIN (PIXELS) AT                     1265.31   1237.43
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.27
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     466 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   46     31     54     62     52     28     52     68     58
 SIGMAB (degree)       0.033  0.033  0.033  0.033  0.033  0.034  0.031  0.032  0.032
 SIGMAR (degree)       0.031  0.032  0.031  0.032  0.030  0.026  0.029  0.032  0.033


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    466 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  2  7 11 10  6  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  2  1  0  0  0     0  0  0  0  3 11 17 16  9  4  2  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  2  1  1  0  0  0     0  0  0  0  2  6 11  9  6  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0 -1  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1 -1  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  2  7  9  6  3  1  0  0  0     0  0  0  1  7 14 15  8  4  1  0  0  0     0  0  0  2  6 10  9  5  2  1  0  0  0
     0  0  0  2 10 28 39 30 16  6  2  1  0     0  0  1  5 26 58 66 41 20  8  3  1  0     0  0  1  6 23 44 42 25 11  4  1  0  0
     0  0  0  3 17 47 61 48 25 11  4  1  0     0  0  1 10 44 94100 68 33 14  5  1  0     0  0  1 10 37 69 63 41 19  8  3  1  0
     0  0  0  2 10 29 38 27 15  6  3  1  0     0  0  0  6 26 59 63 40 19  8  3  1  0     0  0  0  6 23 43 41 24 11  5  1  0  0
     0  0  0  0  3  7 10  6  3  1  0  0  0     0  0  0  1  7 15 17  9  4  2  0  0  0     0  0  0  1  6 11 12  5  3  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0 -1  0     0  0  0  0  1  2  1  1  0  0  0 -1  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  3  3  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  8 12  9  6  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4 11 19 15 10  4  2  0  0  0     0  0  0  1  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 12 10  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  3  3  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    620 REFLECTION INTENSITIES INTEGRATED


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
   401   0  1.130  5167387    0    147      86     6  0.04061  0.07732
   402   0  1.135  5166785    0    148      92     5  0.03798  0.07663
   403   0  1.128  5166988    0    144      89     3  0.04219  0.08006
   404   0  1.131  5167937    0    144      84     9  0.04192  0.07557
   405   0  1.129  5167588    0    138      87     6  0.04136  0.08334
   406   0  1.131  5166780    0    148      95     7  0.04122  0.07052
   407   0  1.129  5166429    0    144     100    10  0.03811  0.06226
   408   0  1.127  5167337    0    148      92     4  0.04082  0.09898
   409   0  1.127  5167748    0    147      84     5  0.04257  0.08821
   410   0  1.127  5167424    0    140      94     4  0.04057  0.10673
   411   0  1.129  5167691    0    143      88     5  0.04199  0.10084
   412   0  1.125  5167804    0    141      85     5  0.04191  0.11861
   413   0  1.128  5167705    0    137      92     7  0.04210  0.10288
   414   0  1.125  5167590    0    142      88     5  0.04438  0.12501
   415   0  1.129  5168456    0    137      79     5  0.04046  0.09042
   416   0  1.124  5168479    0    133      79     6  0.03829  0.07529
   417   0  1.132  5168567    0    124      78     8  0.04053  0.08289
   418   0  1.136  5168955    0    112      77     7  0.04087  0.10204
   419   0  1.133  5168607    0    137      81     6  0.03973  0.11040
   420   0  1.129  5167695    0    130      88     6  0.04133  0.09000

   478 OUT OF    485 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.47
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.423    32.113    34.483  88.672 108.441 111.882
 REC. CELL PARAMETERS   0.041654  0.033754  0.030748  83.823  70.630  67.341
 COORDINATES OF UNIT CELL A-AXIS   -14.331   -23.358    -1.039
 COORDINATES OF UNIT CELL B-AXIS    31.406    -5.395     3.974
 COORDINATES OF UNIT CELL C-AXIS    -1.187    12.099    32.269
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.016    -0.009     0.002
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.040
 LAB COORDINATES OF ROTATION AXIS -0.999995 -0.002850 -0.001317
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000956  0.005096  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.36   1244.02
 DETECTOR ORIGIN (PIXELS) AT                     1265.10   1237.32
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.25
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     478 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   50     31     61     71     55     28     52     65     57
 SIGMAB (degree)       0.033  0.033  0.033  0.033  0.034  0.034  0.031  0.032  0.032
 SIGMAR (degree)       0.033  0.032  0.032  0.033  0.031  0.028  0.030  0.033  0.034


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    478 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  2  4  3  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  2  1  0  0  0     0  0  0  0  4 12 18 15  9  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  4  5  3  2  1  0  0     0  0  0  1  5 17 26 24 14  7  3  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  3  3  2  1  0  0  0     0  0  0  0  3 10 17 15  8  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  1  0  0  0  0     0  0  0  0  1  3  4  3  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0 -1  0     0  0  0  0  1  1  1  0  0  0  0 -1  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  1  3  9 14  9  4  1  0  0  0     0  0  0  2  7 14 16 10  4  1  0  0  0     0  0  0  1  5  8  8  5  2  1  0  0  0
     0  0  0  3 15 41 52 39 21  9  3  1  0     0  0  1  7 32 63 64 42 20  8  3  0  0     0  0  1  6 22 39 35 21  9  4  1  0  0
     0  0  0  4 25 62 79 63 32 15  6  2  0     0  0  1 11 51 97100 68 32 14  5  1  0     0  0  1  9 36 61 55 34 15  6  2  1  0
     0  0  0  3 15 39 51 40 19  9  3  1  0     0  0  1  7 30 63 67 45 20  9  3  1  0     0  0  1  5 20 39 37 23 10  4  1  0  0
     0  0  0  1  4  9 13  9  5  2  0  0  0     0  0  0  2  8 14 18 11  5  2  0  0  0     0  0  0  1  5  9 10  5  3  1  0  0  0
     0  0  0  0  0  1  2  1  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0 -1 -1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 11  8  4  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3 10 16 13  8  3  1  0  0 -1     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 10  9  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  2  3  1  1  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  1  0  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0 -1  0  0  0  0  0  0  0  0  0  0

    653 REFLECTION INTENSITIES INTEGRATED


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
   421   0  1.126  5167297    0    133      93    10  0.04272  0.09166
   422   0  1.127  5168459    0    133      82     9  0.04033  0.07533
   423   0  1.129  5169638    0    130      72     6  0.03908  0.11657
   424   0  1.128  5168743    0    126      79     5  0.04126  0.10386
   425   0  1.128  5168278    0    121      83     5  0.04019  0.09303
   426   0  1.129  5167923    0    130      84     9  0.04010  0.08629
   427   0  1.132  5168455    0    140      76     6  0.04159  0.08093
   428   0  1.131  5168009    0    132      83     6  0.04042  0.09809
   429   0  1.129  5167889    0    143      80     8  0.04003  0.07623
   430   0  1.129  5167410    0    129      87     9  0.03950  0.08409
   431   0  1.132  5167861    0    133      83     5  0.04039  0.09697
   432   0  1.132  5167667    0    146      91    10  0.04047  0.09298
   433   0  1.132  5166869    0    148      92     7  0.03977  0.08952
   434   0  1.129  5165923    0    151     111    11  0.03973  0.08367
   435   0  1.129  5167206    0    150      94     8  0.04154  0.09605
   436   0  1.135  5167539    0    147      92    10  0.04014  0.08738
   437   0  1.133  5167052    0    151      95    10  0.03980  0.07588
   438   0  1.131  5166889    0    150      97     9  0.03928  0.06909
   439   0  1.132  5167248    0    151      94     8  0.03819  0.07981
   440   0  1.130  5167725    0    146      91     5  0.03872  0.08562

   485 OUT OF    491 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.54
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.13
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.419    32.107    34.522  88.700 108.432 111.882
 REC. CELL PARAMETERS   0.041661  0.033761  0.030713  83.795  70.628  67.332
 COORDINATES OF UNIT CELL A-AXIS   -14.330   -23.353    -1.038
 COORDINATES OF UNIT CELL B-AXIS    31.402    -5.395     3.963
 COORDINATES OF UNIT CELL C-AXIS    -1.195    12.114    32.305
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.014     0.001     0.004
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.040
 LAB COORDINATES OF ROTATION AXIS -0.999995 -0.002850 -0.001317
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.001193  0.005112  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.39   1244.05
 DETECTOR ORIGIN (PIXELS) AT                     1264.83   1237.33
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.27
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     485 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   51     33     63     68     51     28     58     67     56
 SIGMAB (degree)       0.033  0.033  0.033  0.034  0.035  0.034  0.031  0.032  0.032
 SIGMAR (degree)       0.035  0.036  0.034  0.033  0.030  0.027  0.029  0.033  0.037


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    485 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  3  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  2  7 10  9  5  3  1  0  0
     0  0  0  0  0  0  0  0  1  0  0  0  0     0  0  0  0  0  1  2  2  1  1  0  0  0     0  0  0  0  3 10 15 14  8  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  1  0  0  0  0     0  0  0  0  2  6  9  9  5  2  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0 -1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0 -1  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  1  0  0  0 -1  0  0     0  0  0  0  0  1  1  0  0  0 -1 -1  0     0  0  0  0  0  1  1  0  0  0 -1 -1  0
     0  0  0  1  4  9 10  7  3  1  0  0  0     0  0  0  2  8 15 15  9  4  1  0  0  0     0  0  0  2  9 15 11  6  3  1  0  0  0
     0  0  0  1 12 30 41 31 17  8  3  1  0     0  0  0  5 28 57 63 41 21  9  3  1  0     0  0  1  7 29 49 46 28 12  5  2  0  0
     0  0  0  3 18 46 62 49 26 12  5  1  0     0  0  0  9 44 90100 71 34 15  6  1  0     0  0  1 11 44 74 74 47 21  9  3  1  0
     0  0  0  2 11 29 37 29 16  7  2  0  0     0  0  0  6 26 57 63 42 22  9  3  1  0     0  0  1  7 26 47 47 27 14  5  2  0  1
     0  0  0  0  2  7 10  7  4  1  0  0  0     0  0  0  1  6 14 17 11  6  2  0  0  0     0  0  0  1  6 12 12  8  4  1  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0 -1  0  0  0  0  0 -1 -1  0  0  0


     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  6  9  6  3  1  0  0  0  0     0  0  0  1  2  3  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  4 17 26 21 13  5  2  0  0  0     0  0  0  2  5  6  4  3  1  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0
     0  0  1  6 23 35 31 18  7  3  1  0  0     0  0  0  1  5  7  6  3  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  3 12 20 19 10  5  2  1  0  1     0  0  0  1  3  4  3  2  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  6  5  3  1  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
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
   441   0  1.131  5168581    0    134      81     7  0.04125  0.08638
   442   0  1.127  5169127    0    125      73     7  0.03897  0.09417
   443   0  1.125  5168541    0    123      78     8  0.04162  0.07827
   444   0  1.129  5167798    0    138      85     7  0.04350  0.09328
   445   0  1.132  5168099    0    139      82     4  0.04323  0.09152
   446   0  1.127  5168082    0    135      84     8  0.03920  0.09379
   447   0  1.125  5168033    0    131      87     3  0.03977  0.10930
   448   0  1.126  5167676    0    135      86     6  0.04179  0.08791
   449   0  1.124  5166961    0    136      95    10  0.04060  0.07528
   450   0  1.123  5167094    0    142      89     2  0.04116  0.09000
   451   0  1.120  5167480    0    143      86     2  0.04166  0.09177
   452   0  1.123  5168210    0    133      81     4  0.03989  0.09289
   453   0  1.120  5168677    0    129      72     4  0.04118  0.09965
   454   0  1.119  5168179    0    128      79     4  0.04440  0.09214
   455   0  1.118  5168179    0    132      83     5  0.03858  0.09319
   456   0  1.115  5166918    0    146      91     5  0.04198  0.08109
   457   0  1.117  5166592    0    153      97     9  0.03756  0.06746
   458   0  1.120  5166818    0    150      92     8  0.03777  0.06914
   459   0  1.116  5167874    0    148      83     7  0.04164  0.08092
   460   0  1.120  5167203    0    142      90     4  0.04030  0.10725

   475 OUT OF    480 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.55
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.16
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.421    32.107    34.490  88.662 108.451 111.891
 REC. CELL PARAMETERS   0.041662  0.033762  0.030744  83.827  70.622  67.334
 COORDINATES OF UNIT CELL A-AXIS   -14.330   -23.356    -1.040
 COORDINATES OF UNIT CELL B-AXIS    31.401    -5.388     3.983
 COORDINATES OF UNIT CELL C-AXIS    -1.192    12.110    32.272
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.001     0.000     0.001
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.040
 LAB COORDINATES OF ROTATION AXIS -0.999995 -0.002850 -0.001317
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000790  0.005160  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.34   1244.01
 DETECTOR ORIGIN (PIXELS) AT                     1265.30   1237.23
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.23
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     475 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   50     33     58     66     55     31     47     65     54
 SIGMAB (degree)       0.033  0.034  0.033  0.034  0.035  0.035  0.032  0.032  0.032
 SIGMAR (degree)       0.032  0.035  0.032  0.033  0.030  0.027  0.029  0.032  0.036


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    475 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  2  3  3  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  2  2  1  1  0  0  0     0  0  0  0  2  8 14 13  7  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  3  3  2  1  0  0  0     0  0  0  0  4 13 22 20 12  6  2  1  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  1  2  2  2  1  0  0  0     0  0  0  0  2  9 15 12  8  3  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  1  3  3  4  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0 -1  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0 -1 -1  0     0  0  0  0  1  1  1  1  0  0  0 -1  0
     0  0  0  1  4 10 10  8  4  1  0  0  0     0  0  0  2  8 16 14 10  4  1  0  0  0     0  0  0  2  7 12  9  5  2  1  0  0  0
     0  0  0  2 14 35 47 37 18  8  3  1  0     0  0  1  7 30 60 64 44 20  9  3  1  0     0  0  1  7 24 43 40 24 10  5  1  0 -1
     0  0  0  5 21 53 71 55 30 14  6  2  0     0  0  1 12 46 93100 67 34 15  6  2  0     0  0  1 11 36 66 60 37 17  7  2  0  0
     0  0  0  3 14 33 46 34 19  9  3  1  0     0  0  1  7 30 59 65 42 22 10  3  1  0     0  0  1  6 22 39 39 24 11  5  1  0  0
     0  0  0  1  4  9 11 10  5  2  0  0  0     0  0  0  2  7 15 16 13  6  2  0  0  0     0  0  0  1  5 10 10  8  3  1  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  1  1  2  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1 -1 -1  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  4  3  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  9 13 11  6  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  1  4 11 20 16  9  4  2  1  0  0     0  0  0  0  1  2  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  7 11 10  6  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  3  3  2  1  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    643 REFLECTION INTENSITIES INTEGRATED


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
   461   0  1.116  5167083    0    127      90     5  0.03957  0.08194
   462   0  1.114  5167830    0    129      79     7  0.04171  0.08095
   463   0  1.116  5168977    0    125      67     3  0.04265  0.09889
   464   0  1.111  5169082    0    127      67     2  0.04131  0.11119
   465   0  1.114  5168286    0    135      77     2  0.04178  0.10100
   466   0  1.112  5167238    0    139      90     8  0.04170  0.08107
   467   0  1.113  5167809    0    144      87     6  0.04126  0.09621
   468   0  1.113  5167740    0    136      90     7  0.04110  0.08856
   469   0  1.110  5167565    0    140      89     7  0.03961  0.07504
   470   0  1.110  5167001    0    138      94     9  0.04293  0.08903
   471   0  1.116  5167883    0    140      81     9  0.04110  0.08202
   472   0  1.120  5167808    0    129      83     6  0.04200  0.09305
   473   0  1.115  5168088    0    121      78     3  0.04247  0.09969
   474   0  1.112  5168164    0    125      84     3  0.04425  0.11943
   475   0  1.111  5168023    0    136      85     6  0.04529  0.10773
   476   0  1.107  5166345    0    128      98     8  0.04049  0.07957
   477   0  1.101  5166855    0    139      93    10  0.04318  0.07482
   478   0  1.096  5167538    0    138      88     7  0.04423  0.10324
   479   0  1.088  5167752    0    127      88     5  0.04554  0.10933
   480   0  1.092  5167420    0    125      88     4  0.04001  0.10676

   466 OUT OF    470 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.52
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.427    32.114    34.492  88.642 108.468 111.891
 REC. CELL PARAMETERS   0.041656  0.033754  0.030744  83.842  70.612  67.339
 COORDINATES OF UNIT CELL A-AXIS   -14.330   -23.362    -1.049
 COORDINATES OF UNIT CELL B-AXIS    31.406    -5.386     3.995
 COORDINATES OF UNIT CELL C-AXIS    -1.193    12.110    32.274
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.007    -0.008    -0.002
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.040
 LAB COORDINATES OF ROTATION AXIS -0.999995 -0.002850 -0.001317
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000493  0.005093  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.33   1244.01
 DETECTOR ORIGIN (PIXELS) AT                     1265.69   1237.32
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.29
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     466 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   48     31     55     66     51     27     55     64     53
 SIGMAB (degree)       0.033  0.033  0.033  0.034  0.035  0.035  0.032  0.032  0.032
 SIGMAR (degree)       0.032  0.034  0.033  0.032  0.030  0.027  0.029  0.033  0.035


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    466 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1 -1 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0     0  0  0  0  1  3  4  3  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  2  1  0  0  0     0  0  0  0  3 11 17 15  8  4  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  4  5  3  2  1  0  0     0  0  0  1  5 16 26 23 14  6  2  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  3  3  2  1  0  0  0     0  0  0  0  4 10 16 13  8  4  1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  1  0  0  0  0     0  0  0  0  1  3  4  4  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0 -1  0 -1 -1  0  0  0
     0  0  0  0  0  0  0  0 -1 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  1  0  0 -1 -1  0     0  0  0  0  1  2  1  1  0  0 -1 -1  0     0  0  0  0  1  1  1  0  0  0  0 -1  0
     0  0  0  1  4 10 12  8  4  1  0  0  0     0  0  0  2  9 15 14  9  4  1  0  0  0     0  0  0  2  6 10  8  4  2  0  0  0  0
     0  0  0  3 15 41 52 38 18  9  3  1  0     0  0  0  7 31 65 67 44 18  9  3  1  0    -1  0  1  6 23 40 36 22  9  4  1  0  0
     0  0  0  5 24 61 76 58 30 13  5  2  0     0  0  1 12 48 99100 66 32 14  5  1  0     0  0  1 11 36 63 56 35 16  7  2  0  0
     0  0  0  3 16 39 48 33 18  8  3  1  0     0  0  1  7 32 61 64 40 20  9  3  1  0     0  0  1  6 23 39 36 22 11  4  1  0  0
     0  0  0  1  4 11 12  9  4  2  0  0  0     0  0  0  2  7 17 16 11  5  2  0  0  0     0  0  0  2  5 11  9  6  3  1  0  0  0
     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1 -1  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0 -1 -1  0  0
     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0  0 -1  0  0  0  0


     0  0  0  0  0  0 -1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  0  0  0  0  0  0  0  0     0  0  0  0  1  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  2  2  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 10  7  4  1  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  9 15 12  7  3  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6  9  8  4  2  1  0  0  0     0  0  0  0  1  1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  2  3  2  1  1  0  0 -1  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0 -1  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    617 REFLECTION INTENSITIES INTEGRATED


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
   481   0  1.089  5167823    0    135      83     3  0.04254  0.10420
   482   0  1.086  5167635    0    143      86     2  0.04450  0.12237
   483   0  1.086  5166318    0    138      98     2  0.04457  0.12627
   484   0  1.083  5167339    0    142      91     4  0.04233  0.10074
   485   0  1.080  5166862    0    141      95     4  0.04283  0.11636
   486   0  1.083  5167142    0    139      97     6  0.04187  0.11807
   487   0  1.081  5167673    0    141      98     0  0.04301  0.17459
   488   0  1.081  5167884    0    135      91     5  0.04082  0.12273
   489   0  1.075  5167846    0    128      88     8  0.04012  0.08944
   490   0  1.073  5168098    0    138      84     9  0.04127  0.08302
   491   0  1.075  5168151    0    140      83     7  0.03957  0.08245
   492   0  1.073  5167932    0    132      87     8  0.03971  0.09181
   493   0  1.069  5168332    0    139      82     6  0.04267  0.10657
   494   0  1.070  5166840    0    138     101     3  0.03836  0.10822
   495   0  1.067  5166918    0    130      97     7  0.03884  0.07690
   496   0  1.068  5168333    0    140      78     3  0.03981  0.08168
   497   0  1.067  5168508    0    132      75     6  0.03998  0.06763
   498   0  1.065  5168066    0    137      81     8  0.03924  0.06810
   499   0  1.064  5167442    0    140      91     2  0.03826  0.10466
   500   0  1.062  5166760    0    137      95     2  0.04178  0.09887

   467 OUT OF    469 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   POSITION BEAM ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.50
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.14
 SPACE GROUP NUMBER      1
 UNIT CELL PARAMETERS     27.425    32.111    34.492  88.674 108.438 111.889
 REC. CELL PARAMETERS   0.041654  0.033758  0.030739  83.819  70.631  67.333
 COORDINATES OF UNIT CELL A-AXIS   -14.331   -23.359    -1.034
 COORDINATES OF UNIT CELL B-AXIS    31.405    -5.388     3.978
 COORDINATES OF UNIT CELL C-AXIS    -1.194    12.111    32.274
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.002     0.003    -0.000
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.040
 LAB COORDINATES OF ROTATION AXIS -0.999995 -0.002850 -0.001317
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.000788  0.005165  1.023528
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1266.39   1244.05
 DETECTOR ORIGIN (PIXELS) AT                     1265.35   1237.26
 CRYSTAL TO DETECTOR DISTANCE (mm)       231.30
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     467 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                   46     33     61     67     49     26     54     63     54
 SIGMAB (degree)       0.033  0.033  0.034  0.034  0.035  0.034  0.032  0.032  0.032
 SIGMAR (degree)       0.033  0.036  0.033  0.032  0.030  0.027  0.029  0.033  0.037


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    467 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1 -1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  3  4  4  1  1  0 -1  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  1  3  3  2  1  0  0  0     0  0  0  0  3 11 17 14  9  4  1  0  0
     0  0  0  0  0  0  0  1  0  0  0  0  0     0  0  0  0  0  2  5  5  3  2  1  0  0     0  0  0  1  5 17 27 24 14  7  3  1  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  2  3  2  1  0  0  0     0  0  0  1  4 11 16 15  9  4  2  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  1  0  0  0  0  0     0  0  0  0  1  3  4  4  2  1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0 -1  0 -1  0  0  0     0  0  0 -1  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0 -1  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0
     0  0  0  0  0  1  1  0  0  0 -1 -1  0     0  0  0  0  1  1  1  0  0  0 -1 -1  0     0  0  0  0  1  1  1  0  0  0  0  0  0
     0  0  0  1  5 11 13  9  4  2  0 -1  0     0  0  0  2  9 16 17  9  4  1  0 -1  0     0  0  0  2  7 11 11  5  2  1  0  0  0
     0  0  0  4 17 44 53 38 20  9  3  1  0     0  0  1  9 32 67 67 40 20  9  3  1  0     0  0  1  8 24 44 40 21 10  4  1  0  0
     0  0  1  5 26 64 80 60 32 14  6  1  0     0  0  1 13 49 95100 66 33 14  6  2  0     0  0  1 10 35 61 58 35 16  7  3  1  0
     0  0  0  3 16 41 48 38 20  9  3  1  0     0  0  1  7 30 62 62 42 21  9  3  1  0     0  0  1  5 21 38 35 21 10  4  2  0  0
     0  0  0  1  4  9 12  9  5  2  0  0  0     0  0  0  2  7 14 15 10  5  2  0  0  0     0  0  0  1  5  8  9  5  3  1  0  0  0
     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  1  2  2  1  0  0  0  0  0     0  0  0  0  0  1  1  1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0 -1  0  0    -1  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0 -1  0  0  0


     0  0  0  0  0  0  0 -1  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  1  3  3  3  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  8 12  9  5  2  1  0  0  0     0  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3 10 16 13  7  3  1  0  0  0     0  0  0  0  1  1  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 10  8  4  2  1  0  0  0     1  0  0  0  1  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  1  2  2  1  1  0  0  0  0     0  0  0  0  0  1  1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0

    631 REFLECTION INTENSITIES INTEGRATED


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
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  3  4  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  0  4 12 15 11  5  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  4  4  0  0  0  0  0     0  0  0  2  6 16 22 17  9  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  0  4 10 14 11  6  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  4  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  5 10 12  8  3  0  0  0  0     0  0  0  2  9 16 16 10  4  0  0  0  0     0  0  0  2  6 10 10  5  2  0  0  0  0
     0  0  0  4 17 41 49 33 15  6  0  0  0     0  0  0  7 32 65 68 42 18  7  0  0  0     0  0  0  6 23 42 40 23  9  3  0  0  0
     0  0  0  5 25 59 71 51 25  9  0  0  0     0  0  0 11 46 94100 66 30 11  0  0  0     0  0  0  9 34 60 58 36 16  6  0  0  0
     0  0  0  3 15 38 46 32 15  6  0  0  0     0  0  0  6 29 62 67 42 18  7  0  0  0     0  0  0  6 21 40 40 23 10  4  0  0  0
     0  0  0  0  4 10 13  8  4  0  0  0  0     0  0  0  2  7 17 19 12  5  0  0  0  0     0  0  0  0  6 11 12  7  3  0  0  0  0
     0  0  0  0  0  0  3  0  0  0  0  0  0     0  0  0  0  0  3  3  2  0  0  0  0  0     0  0  0  0  0  3  3  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  5  3  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4  7 12 10  5  3  0  0  0  0     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3 10 16 14  8  4  2  0  0  0     0  0  0  0  4  3  3  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  6 11 10  5  3  0  0  0  0     0  0  0  0  0  4  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  3  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
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
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  2  0  0  0  0  0     0  0  0  0  3  6  8  6  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  4  5  0  0  0  0  0     0  0  0  0  4  9 13 11  6  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  4  5  0  0  0  0  0     0  0  0  0  3  6 10  8  5  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  4  4  4  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  2  0  0  0  0  0  0  0     0  0  0  0  3  3  2  2  0  0  0  0  0     0  0  0  0  0  3  2  0  0  0  0  0  0
     0  0  0  3  4  7  8  5  3  0  0  0  0     0  0  0  3  9 18 18 11  5  0  0  0  0     0  0  0  3 10 17 15  8  4  0  0  0  0
     0  0  0  3 11 28 35 23 11  5  0  0  0     0  0  0  7 32 70 72 44 20  8  0  0  0     0  0  0  8 33 61 55 31 14  6  0  0  0
     0  0  0  3 16 40 51 37 19  8  0  0  0     0  0  0  8 43 94100 64 31 12  0  0  0     0  0  0  9 42 77 71 43 20  8  0  0  0
     0  0  0  2 10 27 34 24 12  5  0  0  0     0  0  0  5 26 58 63 40 19  8  0  0  0     0  0  0  5 22 43 42 25 12  5  0  0  0
     0  0  0  0  4 10 12  8  4  0  0  0  0     0  0  0  2  8 18 20 12  6  3  0  0  0     0  0  0  3  6 11 12  7  4  2  0  0  0
     0  0  0  0  0  0  4  0  0  0  0  0  0     0  0  0  0  0  4  3  2  0  0  0  0  0     0  0  0  0  0  3  2  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  5  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6  6  5  4  3  0  0  0  0     0  0  0  0  7 11  5  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  5 11 18 15  8  5  3  0  0  0     0  0  0  6  6  6  7 10  4  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4 13 22 19 11  5  3  0  0  0     0  0  0  4  5  7  5  4  4  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  6 12 11  7  4  3  0  0  0     0  0  0  0  4  5  4  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  4  4  3  2  0  0  0  0     0  0  0  0  0  2  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  3 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  3  8 10  6  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  0  0  0  0  0     0  0  0  0  5 13 16 11  5  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  0  0  0  0  0     0  0  0  0  4 10 13  9  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  2  4  5  4  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  0  0  0  0  0  0  0
     0  0  0  0  3  8  9  5  3  0  0  0  0     0  0  0  0  6 14 14  8  3  0  0  0  0     0  0  0  0  5 11 10  5  2  0  0  0  0
     0  0  0  3 12 30 34 22 10  4  0  0  0     0  0  0  4 23 55 58 36 17  7  0  0  0     0  0  0  4 19 42 42 25 12  5  0  0  0
     0  0  0  3 21 50 58 40 19  8  0  0  0     0  0  0  7 41 91100 66 32 13  0  0  0     0  0  0  7 34 69 72 47 22  9  0  0  0
     0  0  0  3 16 39 47 32 15  7  0  0  0     0  0  0  6 30 70 79 52 25 11  0  0  0     0  0  0  5 25 52 57 37 18  8  0  0  0
     0  0  0  0  6 15 18 12  6  3  0  0  0     0  0  0  2 11 25 29 19  9  4  0  0  0     0  0  0  2  9 18 21 14  7  3  0  0  0
     0  0  0  0  0  3  4  3  0  0  0  0  0     0  0  0  0  3  5  5  3  2  0  0  0  0     0  0  0  0  2  3  4  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  4  3  3  4  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  6 12 11  7  4  3  0  0  0     0  0  0  0  4  3  3  4  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3 11 22 21 14  6  3  0  0  0     0  0  0  0  4  4  4  4  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  8 17 18 11  5  4  0  0  0     0  0  0  0  3  4  4  3  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  6  7  5  3  0  0  0  0     0  0  0  0  0  3  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  4 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  5  7  5  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0     0  0  0  0  5 15 21 15  8  5  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  4  4  3  0  0  0  0     0  0  0  0  6 17 24 19 10  5  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  2  0  0  0  0  0     0  0  0  0  3  9 13 10  5  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  2  2  0  0  0  0  0  0     0  0  0  0  0  2  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  6 16 20 12  6  3  0  0  0     0  0  0  2  9 21 23 13  6  2  0  0  0     0  0  0  0  6 13 12  7  3  0  0  0  0
     0  0  0  3 17 48 59 41 20  9  0  0  0     0  0  0  5 27 67 75 48 23 10  0  0  0     0  0  0  4 19 42 43 26 12  5  0  0  0
     0  0  0  3 22 58 73 52 26 12  0  0  0     0  0  0  7 39 89100 66 32 14  6  0  0     0  0  0  6 30 61 62 39 18  8  0  0  0
     0  0  0  3 12 32 41 29 14  6  0  0  0     0  0  0  5 23 54 61 40 19  8  0  0  0     0  0  0  4 19 39 41 25 12  5  0  0  0
     0  0  0  0  3  8 10  7  4  0  0  0  0     0  0  0  2  6 14 16 11  5  3  0  0  0     0  0  0  0  6 12 12  8  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  0  0  0  4  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  4  3  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  6 12 11  6  3  2  0  0  0     0  0  0  0  0  2  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3 10 18 17 10  5  3  0  0  0     0  0  0  0  2  3  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  7 13 12  7  3  2  0  0  0     0  0  0  0  0  3  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  5  4  3  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  5 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  4  5  4  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  0  0  0  0  0     0  0  0  0  4 11 17 14  8  5  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  4  4  0  0  0  0     0  0  0  0  4 14 22 19 11  6  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  3  0  0  0  0  0     0  0  0  0  3  7 11  9  6  5  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  5 12 16 10  5  0  0  0  0     0  0  0  3  8 17 19 11  5  2  0  0  0     0  0  0  2  6 10 10  5  3  0  0  0  0
     0  0  0  3 15 42 55 41 20  9  0  0  0     0  0  2  6 29 65 72 47 22 10  4  0  0     0  0  0  6 21 40 39 23 11  4  0  0  0
     0  0  0  4 19 54 72 54 29 13  7  0  0     0  0  2  9 41 91100 67 33 15  6  0  0     0  0  2  9 32 59 57 35 17  7  0  0  0
     0  0  0  3 10 29 38 29 16  7  5  0  0     0  0  0  6 25 55 59 39 20  9  4  0  0     0  0  0  6 21 38 36 21 10  4  0  0  0
     0  0  0  0  3  6  8  6  4  2  0  0  0     0  0  0  2  7 13 14  9  4  3  0  0  0     0  0  0  0  6 11 10  5  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  2  3  3  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 10  9  5  3  0  0  0  0     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4  9 15 13  7  4  3  0  0  0     0  0  0  0  3  3  3  2  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  6 10  8  5  3  0  0  0  0     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  5  3  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  6 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  6  9  7  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  0  0  0  0  0  0     0  0  0  3  3  9 13 10  6  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  5  7  6  4  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  2  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  4  9 11  7  3  0  0  0  0     0  0  0  3  9 17 16  9  4  0  0  0  0     0  0  0  3  5  8  6  4  3  0  0  0  0
     0  0  0  3 12 34 41 27 13  6  0  0  0     0  0  3  8 35 68 65 38 17  7  0  0  0     0  0  0  6 22 35 29 16  7  3  0  0  0
     0  0  0  4 18 48 57 40 20  9  0  0  0     0  0  3 12 52100 95 58 27 12  4  0  0     0  0  2 11 35 55 46 26 12  5  0  0  0
     0  0  0  3 11 27 32 23 11  5  0  0  0     0  0  3  8 32 61 58 35 17  7  0  0  0     0  0  0  7 23 36 31 17  8  3  0  0  0
     0  0  0  0  3  6  7  5  3  0  0  0  0     0  0  0  2  8 16 15  9  4  3  0  0  0     0  0  0  2  7 10  9  5  3  0  0  0  0
     0  0  0  0  0  0  3  0  0  0  0  0  0     0  0  0  0  0  3  2  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  4  6  4  3  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  6  9  7  4  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  5  6  5  3  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  7 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  2  0  0  0  0  0  0     0  0  0  0  3  6 11  9  5  3  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  4  4  3  0  0  0  0     0  0  0  0  4 14 22 19 10  5  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  0  0  0  0  0     0  0  0  0  3 10 15 12  6  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  4  6  4  2  0  0  0  0     0  0  0  0  4  7  8  5  3  0  0  0  0     0  0  0  0  3  4  3  2  0  0  0  0  0
     0  0  0  0  8 28 39 28 14  6  0  0  0     0  0  0  4 20 48 53 33 15  7  0  0  0     0  0  0  3 13 25 24 14  6  3  0  0  0
     0  0  0  3 20 57 76 55 28 12  0  0  0     0  0  0  8 42 94100 65 31 13  5  0  0     0  0  0  6 27 51 47 28 13  5  0  0  0
     0  0  0  2 14 39 49 34 16  7  0  0  0     0  0  0  5 28 62 65 41 19  8  0  0  0     0  0  0  4 18 34 31 18  8  3  0  0  0
     0  0  0  0  3  8 10  6  3  0  0  0  0     0  0  0  0  5 13 13  8  3  0  0  0  0     0  0  0  0  4  7  6  4  2  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  4  4  3  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  2  6 10  8  5  3  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  4  7  6  4  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  3  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
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
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  0  0  0  0  0     0  0  0  0  3 10 15 14  8  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  5  5  4  0  0  0  0     0  0  0  0  7 21 32 27 16  7  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  4  4  4  3  0  0  0  0     0  0  0  0  5 16 23 18 10  5  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  3  0  0  0  0  0  0     0  0  0  0  2  4  5  4  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  3  5  7  5  3  0  0  0  0     0  0  0  0  4  8  8  6  3  0  0  0  0     0  0  0  0  4  5  5  3  0  0  0  0  0
     0  0  0  3 12 34 45 35 18  8  0  0  0     0  0  0  5 22 50 55 37 18  8  0  0  0     0  0  0  4 16 29 29 18  8  3  0  0  0
     0  0  0  5 26 68 86 65 34 15  6  0  0     0  0  0  9 44 94100 66 33 14  6  0  0     0  0  0  7 29 53 50 31 15  6  0  0  0
     0  0  0  4 19 48 57 40 20  9  0  0  0     0  0  0  6 30 63 64 41 19  8  0  0  0     0  0  0  5 19 35 32 19  8  4  0  0  0
     0  0  0  0  4 10 11  7  3  0  0  0  0     0  0  0  2  6 13 12  7  3  0  0  0  0     0  0  0  0  4  7  6  3  3  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  5  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4  6  8  6  4  3  0  0  0  0     0  0  0  0  0  3  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4  8 13 11  6  3  2  0  0  0     0  0  0  0  0  4  3  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  5  8  7  4  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  2  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  9 *****
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  3  3  4  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  3  3  0  0  0  0  0     0  0  0  0  4 11 17 14  8  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  5  5  4  0  0  0  0     0  0  0  0  6 17 26 24 14  7  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  4  4  0  0  0  0  0     0  0  0  0  4  9 14 13  7  4  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  3  3  4  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  2  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  5 10 11  8  4  0  0  0  0     0  0  0  3 10 17 17 10  5  2  0  0  0     0  0  0  3  8 13 11  6  3  0  0  0  0
     0  0  0  4 18 44 56 40 21  9  0  0  0     0  0  2  9 36 74 79 49 23 10  0  0  0     0  0  0  8 30 51 46 26 11  5  0  0  0
     0  0  0  4 23 59 77 59 31 14  0  0  0     0  0  0  9 44 94100 65 32 13  6  0  0     0  0  0  9 35 61 55 32 15  6  0  0  0
     0  0  0  3 10 28 36 27 14  6  0  0  0     0  0  0  3 18 41 44 28 13  5  0  0  0     0  0  0  3 13 25 23 13  6  3  0  0  0
     0  0  0  0  2  4  6  4  3  0  0  0  0     0  0  0  0  3  6  7  4  2  0  0  0  0     0  0  0  0  2  4  4  3  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4  7  6  4  4  3  0  0  0  0     0  0  0  0  0  8  4  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  5 11 15 11  6  4  3  0  0  0     0  0  0  5  9  6  5  6  2  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  4 11 16 13  7  4  3  0  0  0     0  0  0  3  4  5  5  5  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  3  5  6  5  4  3  0  0  0  0     0  0  0  0  3  4  2  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0  0


 STANDARD DEVIATION (DEGREES) OF BEAM DIVERGENCE
 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1    0.023  0.024  0.024  0.023  0.023  0.023  0.022  0.023  0.023
     2    0.023  0.024  0.024  0.024  0.023  0.023  0.022  0.023  0.023
     3    0.023  0.024  0.024  0.023  0.023  0.023  0.022  0.023  0.023
     4    0.023  0.024  0.024  0.024  0.023  0.023  0.023  0.023  0.023
     5    0.023  0.023  0.024  0.023  0.023  0.024  0.023  0.023  0.023
     6    0.023  0.023  0.024  0.023  0.023  0.023  0.022  0.023  0.023
     7    0.023  0.023  0.023  0.023  0.023  0.023  0.022  0.023  0.023
     8    0.023  0.024  0.024  0.023  0.023  0.023  0.022  0.023  0.023
     9    0.023  0.024  0.024  0.023  0.023  0.023  0.022  0.023  0.023
    10    0.023  0.023  0.024  0.023  0.023  0.023  0.022  0.022  0.023
    11    0.023  0.024  0.024  0.023  0.023  0.023  0.022  0.023  0.023
    12    0.023  0.023  0.024  0.023  0.023  0.023  0.022  0.023  0.023
    13    0.023  0.023  0.024  0.023  0.023  0.023  0.022  0.022  0.023
    14    0.023  0.024  0.024  0.023  0.023  0.023  0.022  0.023  0.023
    15    0.023  0.024  0.024  0.023  0.023  0.023  0.022  0.023  0.023
    16    0.023  0.023  0.024  0.023  0.023  0.024  0.022  0.023  0.023
    17    0.024  0.023  0.024  0.023  0.023  0.024  0.022  0.023  0.023
    18    0.023  0.024  0.024  0.023  0.023  0.023  0.022  0.023  0.023
    19    0.023  0.023  0.023  0.023  0.024  0.024  0.023  0.023  0.023
    20    0.023  0.024  0.024  0.023  0.024  0.024  0.022  0.023  0.023
    21    0.023  0.024  0.024  0.024  0.025  0.025  0.022  0.023  0.023
    22    0.024  0.024  0.024  0.024  0.025  0.024  0.022  0.023  0.024
    23    0.023  0.024  0.024  0.024  0.026  0.025  0.023  0.023  0.023
    24    0.024  0.023  0.024  0.024  0.025  0.025  0.023  0.023  0.023
    25    0.024  0.024  0.024  0.024  0.026  0.025  0.023  0.023  0.023


 STANDARD DEVIATION (DEGREES) OF REFLECTING RANGE
 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1    0.030  0.026  0.029  0.032  0.031  0.025  0.030  0.034  0.030
     2    0.032  0.028  0.031  0.033  0.031  0.025  0.029  0.034  0.031
     3    0.032  0.027  0.032  0.033  0.032  0.027  0.030  0.036  0.033
     4    0.032  0.027  0.030  0.033  0.031  0.026  0.031  0.038  0.034
     5    0.035  0.026  0.031  0.034  0.034  0.029  0.033  0.039  0.035
     6    0.033  0.028  0.031  0.033  0.031  0.026  0.031  0.036  0.031
     7    0.032  0.028  0.032  0.033  0.032  0.028  0.031  0.035  0.030
     8    0.033  0.028  0.030  0.034  0.032  0.026  0.032  0.037  0.035
     9    0.032  0.026  0.031  0.033  0.033  0.028  0.029  0.033  0.031
    10    0.030  0.030  0.032  0.033  0.031  0.026  0.029  0.031  0.030
    11    0.033  0.037  0.036  0.034  0.032  0.028  0.031  0.033  0.031
    12    0.030  0.028  0.031  0.034  0.032  0.027  0.029  0.032  0.030
    13    0.031  0.027  0.031  0.034  0.032  0.025  0.028  0.032  0.030
    14    0.030  0.027  0.031  0.034  0.032  0.026  0.029  0.032  0.030
    15    0.031  0.032  0.032  0.033  0.031  0.026  0.029  0.031  0.031
    16    0.030  0.026  0.030  0.033  0.030  0.025  0.028  0.032  0.031
    17    0.033  0.040  0.034  0.032  0.030  0.026  0.030  0.031  0.034
    18    0.030  0.032  0.031  0.033  0.031  0.025  0.028  0.031  0.031
    19    0.030  0.030  0.031  0.031  0.030  0.028  0.029  0.031  0.030
    20    0.030  0.035  0.032  0.033  0.029  0.024  0.027  0.031  0.033
    21    0.033  0.034  0.031  0.034  0.032  0.026  0.031  0.033  0.034
    22    0.035  0.057  0.041  0.033  0.030  0.026  0.029  0.033  0.046
    23    0.032  0.039  0.033  0.033  0.030  0.025  0.028  0.032  0.037
    24    0.031  0.035  0.033  0.034  0.030  0.026  0.028  0.033  0.036
    25    0.033  0.039  0.033  0.033  0.030  0.026  0.029  0.034  0.038


 ***** SUGGESTED VALUES FOR INPUT PARAMETERS *****
 BEAM_DIVERGENCE=   0.23257  BEAM_DIVERGENCE_E.S.D.=   0.02326
 REFLECTING_RANGE=  0.21900  REFLECTING_RANGE_E.S.D.=  0.03129



 NUMBER OF STRONG REFLECTIONS USED FOR LEARNING REFERENCE PROFILES

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1       55     33     61     66     52     29     56     64     56
     2       47     32     58     71     61     29     53     68     58
     3       52     27     48     62     51     30     57     77     58
     4       51     33     66     73     57     30     50     60     54
     5       49     32     49     65     55     28     53     67     57
     6       51     33     65     69     54     31     58     70     59
     7       43     29     51     62     51     22     49     58     56
     8       48     30     63     69     53     29     57     76     63
     9       47     33     58     71     55     31     53     61     54
    10       53     34     59     64     54     29     54     70     61
    11       49     29     54     66     55     31     58     62     51
    12       45     30     53     64     55     28     53     66     59
    13       50     29     54     71     56     32     59     68     55
    14       43     32     61     67     49     23     50     66     58
    15       47     31     57     64     54     29     54     66     59
    16       52     37     60     69     56     33     58     65     59
    17       49     27     52     65     51     28     54     71     58
    18       50     32     61     71     54     28     58     67     56
    19       48     31     55     66     52     27     48     66     56
    20       46     31     54     62     52     28     52     68     58
    21       50     31     61     71     55     28     52     65     57
    22       51     33     63     68     51     28     58     67     56
    23       50     33     58     66     55     31     47     65     54
    24       48     31     55     66     51     27     55     64     53
    25       46     33     61     67     49     26     54     63     54
 margin    1220    784   1438   1675   1337    717   1350   1661   1418



 NUMBER OF GRID POINTS IN REFERENCE PROFILES WITH  SIGNAL >   2.% OF MAXIMUM

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1      137    123    139    142    129    112    112    132    118
     2      144    128    150    144    127    112    112    131    127
     3      143    133    143    140    133    117    120    131    132
     4      148    145    144    151    126    112    122    137    132
     5      153    129    150    150    140    129    129    138    126
     6      138    125    142    142    131    110    118    136    123
     7      139    118    141    137    128    118    116    132    119
     8      142    134    141    149    130    114    116    132    127
     9      141    133    143    151    129    115    116    132    123
    10      131    122    147    143    128    112    112    122    118
    11      145    140    158    141    131    118    117    131    122
    12      131    124    148    150    124    111    109    125    118
    13      130    127    145    141    127    113    112    123    117
    14      129    126    145    148    131    111    115    128    124
    15      140    134    150    144    129    118    114    123    119
    16      133    123    142    138    125    109    111    126    121
    17      150    135    160    142    126    120    112    125    132
    18      129    131    144    132    129    110    109    124    124
    19      128    126    138    128    125    120    119    126    122
    20      128    139    147    134    123    113    111    127    128
    21      143    139    146    145    143    126    115    134    133
    22      157    159    163    144    133    127    112    124    145
    23      140    148    146    150    146    126    117    123    136
    24      142    148    147    147    142    130    114    132    138
    25      143    152    150    142    147    127    118    136    139
 margin     173    196    190    182    179    154    140    159    176

 MINIMUM INTEGRATION VOLUME IS  109 GRID POINTS



 CORRELATIONS BETWEEN REFERENCE PROFILES WITH CENTRAL PROFILE #1

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1     1.00   0.94   0.95   0.99   0.99   0.99   0.96   0.96   0.96
     2     1.00   0.95   0.92   0.97   0.98   0.98   0.94   0.93   0.95
     3     1.00   0.97   0.93   0.95   0.95   0.96   0.96   0.88   0.89
     4     1.00   0.98   0.97   0.98   0.98   0.96   0.96   0.90   0.90
     5     1.00   0.97   0.93   0.94   0.93   0.95   0.93   0.85   0.84
     6     1.00   0.96   0.92   0.96   0.96   0.97   0.95   0.91   0.94
     7     1.00   0.96   0.96   0.98   0.99   0.98   0.94   0.91   0.95
     8     1.00   0.98   0.96   0.98   0.98   0.96   0.96   0.95   0.90
     9     1.00   0.98   0.98   0.98   0.98   0.95   0.98   0.94   0.93
    10     1.00   0.92   0.95   0.98   0.98   0.98   0.98   0.98   0.96
    11     1.00   0.82   0.87   0.99   0.99   0.93   0.91   0.95   0.96
    12     1.00   0.97   0.95   0.97   0.97   0.98   0.95   0.94   0.96
    13     1.00   0.96   0.96   0.98   0.98   0.98   0.96   0.97   0.94
    14     1.00   0.97   0.97   0.98   0.98   0.98   0.97   0.96   0.96
    15     1.00   0.90   0.96   0.98   0.97   0.98   0.98   0.98   0.95
    16     1.00   0.97   0.98   0.98   0.98   0.97   0.97   0.97   0.94
    17     1.00   0.76   0.89   0.98   0.96   0.94   0.96   0.98   0.89
    18     1.00   0.90   0.96   0.98   0.98   0.97   0.99   0.99   0.93
    19     1.00   0.96   0.97   0.99   0.99   0.95   0.96   0.97   0.97
    20     1.00   0.89   0.95   0.98   0.99   0.96   0.98   0.97   0.96
    21     1.00   0.92   0.96   0.98   0.96   0.95   0.97   0.98   0.96
    22     1.00   0.68   0.88   0.98   0.95   0.93   0.97   0.98   0.81
    23     1.00   0.88   0.97   0.98   0.96   0.95   0.99   0.99   0.88
    24     1.00   0.96   0.98   0.97   0.97   0.96   0.96   0.97   0.96
    25     1.00   0.91   0.96   0.99   0.98   0.96   0.96   0.96   0.94
 margin    1.00   0.96   0.97   0.99   0.99   0.98   0.98   0.97   0.97



 CORRELATIONS BETWEEN REFERENCE PROFILES WITH RUN-AVERAGED PROFILES

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1     0.99   0.98   0.99   0.98   0.98   0.98   0.98   0.98   0.99
     2     1.00   0.99   0.99   0.99   0.99   0.99   0.99   0.99   0.99
     3     1.00   0.98   0.99   0.97   0.95   0.98   0.99   0.96   0.94
     4     1.00   0.93   0.98   0.99   0.99   0.99   0.99   0.97   0.93
     5     0.99   0.96   0.99   0.95   0.93   0.97   0.98   0.95   0.91
     6     1.00   0.99   0.98   0.98   0.97   0.99   0.99   0.98   0.97
     7     1.00   0.99   0.99   0.99   0.98   0.99   0.99   0.99   0.99
     8     0.98   0.89   0.97   0.99   0.99   0.99   0.98   0.97   0.90
     9     1.00   0.97   0.99   0.99   0.98   0.95   0.98   1.00   0.99
    10     1.00   0.99   1.00   1.00   0.99   1.00   1.00   0.99   0.99
    11     0.99   0.92   0.95   0.99   1.00   0.95   0.97   0.99   0.98
    12     1.00   0.99   1.00   0.99   0.99   0.99   0.99   1.00   0.99
    13     1.00   0.99   0.99   0.99   0.99   0.99   0.99   1.00   0.99
    14     1.00   0.98   0.98   0.99   0.99   0.99   0.99   1.00   1.00
    15     1.00   0.98   1.00   0.99   0.99   0.99   1.00   0.98   0.97
    16     1.00   0.97   0.97   0.99   1.00   0.99   0.99   0.99   0.99
    17     0.99   0.89   0.97   0.99   0.98   0.96   1.00   0.98   0.94
    18     1.00   0.99   0.99   0.99   0.98   0.99   0.99   0.99   0.99
    19     0.99   1.00   1.00   1.00   0.99   0.95   0.94   0.97   0.98
    20     0.99   0.96   0.99   0.99   0.99   0.99   0.99   0.99   0.96
    21     1.00   0.99   0.99   0.98   0.96   0.96   0.99   0.99   0.99
    22     0.98   0.72   0.90   0.98   0.99   0.99   0.99   0.93   0.79
    23     0.99   0.96   0.99   0.99   0.97   0.99   0.99   0.98   0.95
    24     0.99   0.98   0.96   0.97   0.98   0.99   0.98   0.99   0.99
    25     0.99   0.97   0.99   0.99   0.97   0.97   0.99   1.00   0.98
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

       -2...-1      -0.249     -0.000        -8.          1
       -1... 0      -0.188      0.071        -2.          5
        0... 1       0.256      0.106         5.         72
        1... 2       0.479      0.097        16.        113
        2... 3       0.615      0.085        33.        170
        3... 6       0.727      0.078        87.        694
        6... 9       0.774      0.073       203.        955
        9...12       0.785      0.074       371.       1163
       12...15       0.797      0.075       604.       1169
       15...18       0.799      0.074       883.       1053
       18...21       0.803      0.070      1248.        969
       21...24       0.797      0.074      1651.        797
       24...27       0.803      0.076      2156.        732
       27...30       0.804      0.075      2726.        640
       30...33       0.799      0.080      3512.        575
       33...36       0.791      0.078      4556.        559
       36...39       0.788      0.084      5898.        518
       39...42       0.791      0.082      8287.        560
       42...45       0.793      0.086     12091.        698
       45...48       0.787      0.091     22584.        898
       48...51       0.784      0.105     50537.        765

   13106 REFLECTIONS SAVED ON FILE "INTEGRATE.HKL"


 total elapsed wall-clock time for INTEGRATE      18.4 sec

 TASK     cpu time (sec)    elapsed wall-clock time (sec)
   1         288.2                    18.3
