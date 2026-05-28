      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)                            00000010
C                                                                       00000020
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.                            00000030
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.                  00000040
C     JACK DONGARRA, LINPACK, 3/11/78.                                  00000050
C                                                                       00000060
      DOUBLE PRECISION DX(1),DY(1),DA                                   00000070
      INTEGER I,INCX,INCY,IXIY,M,MP1,N                                  00000080
C                                                                       00000090
      IF(N.LE.0)RETURN                                                  00000100
      IF (DA .EQ. 0.0D0) RETURN                                         00000110
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                               00000120
C                                                                       00000130
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS                00000140
C          NOT EQUAL TO 1                                               00000150
C                                                                       00000160
      IX = 1                                                            00000170
      IY = 1                                                            00000180
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 00000190
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 00000200
      DO 10 I = 1,N                                                     00000210
        DY(IY) = DY(IY) + DA*DX(IX)                                     00000220
        IX = IX + INCX                                                  00000230
        IY = IY + INCY                                                  00000240
   10 CONTINUE                                                          00000250
      RETURN                                                            00000260
C                                                                       00000270
C        CODE FOR BOTH INCREMENTS EQUAL TO 1                            00000280
C                                                                       00000290
C                                                                       00000300
C        CLEAN-UP LOOP                                                  00000310
C                                                                       00000320
   20 M = MOD(N,4)                                                      00000330
      IF( M .EQ. 0 ) GO TO 40                                           00000340
      DO 30 I = 1,M                                                     00000350
        DY(I) = DY(I) + DA*DX(I)                                        00000360
   30 CONTINUE                                                          00000370
      IF( N .LT. 4 ) RETURN                                             00000380
   40 MP1 = M + 1                                                       00000390
      DO 50 I = MP1,N,4                                                 00000400
        DY(I) = DY(I) + DA*DX(I)                                        00000410
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)                            00000420
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)                            00000430
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)                            00000440
   50 CONTINUE                                                          00000450
      RETURN                                                            00000460
      END                                                               00000470
