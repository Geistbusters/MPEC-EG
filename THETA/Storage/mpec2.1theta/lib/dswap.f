      SUBROUTINE  DSWAP (N,DX,INCX,DY,INCY)                             00000010
C                                                                       00000020
C     INTERCHANGES TWO VECTORS.                                         00000030
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL ONE.                     00000040
C     JACK DONGARRA, LINPACK, 3/11/78.                                  00000050
C                                                                       00000060
      DOUBLE PRECISION DX(1),DY(1),DTEMP                                00000070
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N                                 00000080
C                                                                       00000090
      IF(N.LE.0)RETURN                                                  00000100
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                               00000110
C                                                                       00000120
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL       00000130
C         TO 1                                                          00000140
C                                                                       00000150
      IX = 1                                                            00000160
      IY = 1                                                            00000170
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 00000180
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 00000190
      DO 10 I = 1,N                                                     00000200
        DTEMP = DX(IX)                                                  00000210
        DX(IX) = DY(IY)                                                 00000220
        DY(IY) = DTEMP                                                  00000230
        IX = IX + INCX                                                  00000240
        IY = IY + INCY                                                  00000250
   10 CONTINUE                                                          00000260
      RETURN                                                            00000270
C                                                                       00000280
C       CODE FOR BOTH INCREMENTS EQUAL TO 1                             00000290
C                                                                       00000300
C                                                                       00000310
C       CLEAN-UP LOOP                                                   00000320
C                                                                       00000330
   20 M = MOD(N,3)                                                      00000340
      IF( M .EQ. 0 ) GO TO 40                                           00000350
      DO 30 I = 1,M                                                     00000360
        DTEMP = DX(I)                                                   00000370
        DX(I) = DY(I)                                                   00000380
        DY(I) = DTEMP                                                   00000390
   30 CONTINUE                                                          00000400
      IF( N .LT. 3 ) RETURN                                             00000410
   40 MP1 = M + 1                                                       00000420
      DO 50 I = MP1,N,3                                                 00000430
        DTEMP = DX(I)                                                   00000440
        DX(I) = DY(I)                                                   00000450
        DY(I) = DTEMP                                                   00000460
        DTEMP = DX(I + 1)                                               00000470
        DX(I + 1) = DY(I + 1)                                           00000480
        DY(I + 1) = DTEMP                                               00000490
        DTEMP = DX(I + 2)                                               00000500
        DX(I + 2) = DY(I + 2)                                           00000510
        DY(I + 2) = DTEMP                                               00000520
   50 CONTINUE                                                          00000530
      RETURN                                                            00000540
      END                                                               00000550
