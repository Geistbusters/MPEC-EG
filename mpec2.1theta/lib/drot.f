      SUBROUTINE  DROT (N,DX,INCX,DY,INCY,C,S)                          00000010
C                                                                       00000020
C     APPLIES A PLANE ROTATION.                                         00000030
C     JACK DONGARRA, LINPACK, 3/11/78.                                  00000040
C                                                                       00000050
      DOUBLE PRECISION DX(1),DY(1),DTEMP,C,S                            00000060
      INTEGER I,INCX,INCY,IX,IY,N                                       00000070
C                                                                       00000080
      IF(N.LE.0)RETURN                                                  00000090
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                               00000100
C                                                                       00000110
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL       00000120
C         TO 1                                                          00000130
C                                                                       00000140
      IX = 1                                                            00000150
      IY = 1                                                            00000160
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 00000170
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 00000180
      DO 10 I = 1,N                                                     00000190
        DTEMP = C*DX(IX) + S*DY(IY)                                     00000200
        DY(IY) = C*DY(IY) - S*DX(IX)                                    00000210
        DX(IX) = DTEMP                                                  00000220
        IX = IX + INCX                                                  00000230
        IY = IY + INCY                                                  00000240
   10 CONTINUE                                                          00000250
      RETURN                                                            00000260
C                                                                       00000270
C       CODE FOR BOTH INCREMENTS EQUAL TO 1                             00000280
C                                                                       00000290
   20 DO 30 I = 1,N                                                     00000300
        DTEMP = C*DX(I) + S*DY(I)                                       00000310
        DY(I) = C*DY(I) - S*DX(I)                                       00000320
        DX(I) = DTEMP                                                   00000330
   30 CONTINUE                                                          00000340
      RETURN                                                            00000350
      END                                                               00000360
