      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)                 00000010
C                                                                       00000020
C     FORMS THE DOT PRODUCT OF TWO VECTORS.                             00000030
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.                  00000040
C     JACK DONGARRA, LINPACK, 3/11/78.                                  00000050
C                                                                       00000060
      DOUBLE PRECISION DX(1),DY(1),DTEMP                                00000070
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N                                 00000080
C                                                                       00000090
      DDOT = 0.0D0                                                      00000100
      DTEMP = 0.0D0                                                     00000110
      IF(N.LE.0)RETURN                                                  00000120
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                               00000130
C                                                                       00000140
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS                00000150
C          NOT EQUAL TO 1                                               00000160
C                                                                       00000170
      IX = 1                                                            00000180
      IY = 1                                                            00000190
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 00000200
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 00000210
      DO 10 I = 1,N                                                     00000220
        DTEMP = DTEMP + DX(IX)*DY(IY)                                   00000230
        IX = IX + INCX                                                  00000240
        IY = IY + INCY                                                  00000250
   10 CONTINUE                                                          00000260
      DDOT = DTEMP                                                      00000270
      RETURN                                                            00000280
C                                                                       00000290
C        CODE FOR BOTH INCREMENTS EQUAL TO 1                            00000300
C                                                                       00000310
C                                                                       00000320
C        CLEAN-UP LOOP                                                  00000330
C                                                                       00000340
   20 M = MOD(N,5)                                                      00000350
      IF( M .EQ. 0 ) GO TO 40                                           00000360
      DO 30 I = 1,M                                                     00000370
        DTEMP = DTEMP + DX(I)*DY(I)                                     00000380
   30 CONTINUE                                                          00000390
      IF( N .LT. 5 ) GO TO 60                                           00000400
   40 MP1 = M + 1                                                       00000410
      DO 50 I = MP1,N,5                                                 00000420
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) +             00000430
     *   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)00000440
   50 CONTINUE                                                          00000450
   60 DDOT = DTEMP                                                      00000460
      RETURN                                                            00000470
      END                                                               00000480
