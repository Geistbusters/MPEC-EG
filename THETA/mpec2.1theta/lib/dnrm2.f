      DOUBLE PRECISION FUNCTION DNRM2 ( N, DX, INCX)                    00000010
      INTEGER          NEXT                                             00000020
      DOUBLE PRECISION   DX(1), CUTLO, CUTHI, HITEST, SUM, XMAX,ZERO,ONE00000030
      DATA   ZERO, ONE /0.0D0, 1.0D0/                                   00000040
C                                                                       00000050
C     EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX() WITH STORAGE        00000060
C     INCREMENT INCX .                                                  00000070
C     IF    N .LE. 0 RETURN WITH RESULT = 0.                            00000080
C     IF N .GE. 1 THEN INCX MUST BE .GE. 1                              00000090
C                                                                       00000100
C           C.L.LAWSON, 1978 JAN 08                                     00000110
C                                                                       00000120
C     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE       00000130
C     HOPEFULLY APPLICABLE TO ALL MACHINES.                             00000140
C         CUTLO = MAXIMUM OF  DSQRT(U/EPS)  OVER ALL KNOWN MACHINES.    00000150
C         CUTHI = MINIMUM OF  DSQRT(V)      OVER ALL KNOWN MACHINES.    00000160
C     WHERE                                                             00000170
C         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.                 00000180
C         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)               00000190
C         V   = LARGEST  NO.            (OVERFLOW  LIMIT)               00000200
C                                                                       00000210
C     BRIEF OUTLINE OF ALGORITHM..                                      00000220
C                                                                       00000230
C     PHASE 1    SCANS ZERO COMPONENTS.                                 00000240
C     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO        00000250
C     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO                    00000260
C     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M                  00000270
C     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.                 00000280
C                                                                       00000290
C     VALUES FOR CUTLO AND CUTHI..                                      00000300
C     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER    00000310
C     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..                     00000320
C     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE00000330
C                   UNIVAC AND DEC AT 2**(-103)                         00000340
C                   THUS CUTLO = 2**(-51) = 4.44089E-16                 00000350
C     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.          00000360
C                   THUS CUTHI = 2**(63.5) = 1.30438E19                 00000370
C     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.             00000380
C                   THUS CUTLO = 2**(-33.5) = 8.23181D-11               00000390
C     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19                    00000400
C     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /                        00000410
C     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /                        00000420
      DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /                        00000430
C                                                                       00000440
      IF(N .GT. 0) GO TO 10                                             00000450
         DNRM2  = ZERO                                                  00000460
         GO TO 300                                                      00000470
C                                                                       00000480
   10 ASSIGN 30 TO NEXT                                                 00000490
      SUM = ZERO                                                        00000500
      NN = N * INCX                                                     00000510
C                                                 BEGIN MAIN LOOP       00000520
      I = 1                                                             00000530
   20    GO TO NEXT,(30, 50, 70, 110)                                   00000540
   30 IF( DABS(DX(I)) .GT. CUTLO) GO TO 85                              00000550
      ASSIGN 50 TO NEXT                                                 00000560
      XMAX = ZERO                                                       00000570
C                                                                       00000580
C                        PHASE 1.  SUM IS ZERO                          00000590
C                                                                       00000600
   50 IF( DX(I) .EQ. ZERO) GO TO 200                                    00000610
      IF( DABS(DX(I)) .GT. CUTLO) GO TO 85                              00000620
C                                                                       00000630
C                                PREPARE FOR PHASE 2.                   00000640
      ASSIGN 70 TO NEXT                                                 00000650
      GO TO 105                                                         00000660
C                                                                       00000670
C                                PREPARE FOR PHASE 4.                   00000680
C                                                                       00000690
  100 I = J                                                             00000700
      ASSIGN 110 TO NEXT                                                00000710
      SUM = (SUM / DX(I)) / DX(I)                                       00000720
  105 XMAX = DABS(DX(I))                                                00000730
      GO TO 115                                                         00000740
C                                                                       00000750
C                   PHASE 2.  SUM IS SMALL.                             00000760
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.     00000770
C                                                                       00000780
   70 IF( DABS(DX(I)) .GT. CUTLO ) GO TO 75                             00000790
C                                                                       00000800
C                     COMMON CODE FOR PHASES 2 AND 4.                   00000810
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.00000820
C                                                                       00000830
  110 IF( DABS(DX(I)) .LE. XMAX ) GO TO 115                             00000840
         SUM = ONE + SUM * (XMAX / DX(I))**2                            00000850
         XMAX = DABS(DX(I))                                             00000860
         GO TO 200                                                      00000870
C                                                                       00000880
  115 SUM = SUM + (DX(I)/XMAX)**2                                       00000890
      GO TO 200                                                         00000900
C                                                                       00000910
C                                                                       00000920
C                  PREPARE FOR PHASE 3.                                 00000930
C                                                                       00000940
   75 SUM = (SUM * XMAX) * XMAX                                         00000950
C                                                                       00000960
C                                                                       00000970
C     FOR REAL OR D.P. SET HITEST = CUTHI/N                             00000980
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)                         00000990
C                                                                       00001000
   85 HITEST = CUTHI/FLOAT( N )                                         00001010
C                                                                       00001020
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.            00001030
C                                                                       00001040
      DO 95 J =I,NN,INCX                                                00001050
      IF(DABS(DX(J)) .GE. HITEST) GO TO 100                             00001060
   95    SUM = SUM + DX(J)**2                                           00001070
      DNRM2 = DSQRT( SUM )                                              00001080
      GO TO 300                                                         00001090
C                                                                       00001100
  200 CONTINUE                                                          00001110
      I = I + INCX                                                      00001120
      IF ( I .LE. NN ) GO TO 20                                         00001130
C                                                                       00001140
C              END OF MAIN LOOP.                                        00001150
C                                                                       00001160
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.              00001170
C                                                                       00001180
      DNRM2 = XMAX * DSQRT(SUM)                                         00001190
  300 CONTINUE                                                          00001200
      RETURN                                                            00001210
      END                                                               00001220
