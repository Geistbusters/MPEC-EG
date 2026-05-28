      SUBROUTINE  DSCAL(N,DA,DX,INCX)                                   00000010
C                                                                       00000020
C     SCALES A VECTOR BY A CONSTANT.                                    00000030
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.                   00000040
C     JACK DONGARRA, LINPACK, 3/11/78.                                  00000050
C                                                                       00000060
      DOUBLE PRECISION DA,DX(1)                                         00000070
      INTEGER I,INCX,M,MP1,N,NINCX                                      00000080
C                                                                       00000090
      IF(N.LE.0)RETURN                                                  00000100
      IF(INCX.EQ.1)GO TO 20                                             00000110
C                                                                       00000120
C        CODE FOR INCREMENT NOT EQUAL TO 1                              00000130
C                                                                       00000140
      NINCX = N*INCX                                                    00000150
      DO 10 I = 1,NINCX,INCX                                            00000160
        DX(I) = DA*DX(I)                                                00000170
   10 CONTINUE                                                          00000180
      RETURN                                                            00000190
C                                                                       00000200
C        CODE FOR INCREMENT EQUAL TO 1                                  00000210
C                                                                       00000220
C                                                                       00000230
C        CLEAN-UP LOOP                                                  00000240
C                                                                       00000250
   20 M = MOD(N,5)                                                      00000260
      IF( M .EQ. 0 ) GO TO 40                                           00000270
      DO 30 I = 1,M                                                     00000280
        DX(I) = DA*DX(I)                                                00000290
   30 CONTINUE                                                          00000300
      IF( N .LT. 5 ) RETURN                                             00000310
   40 MP1 = M + 1                                                       00000320
      DO 50 I = MP1,N,5                                                 00000330
        DX(I) = DA*DX(I)                                                00000340
        DX(I + 1) = DA*DX(I + 1)                                        00000350
        DX(I + 2) = DA*DX(I + 2)                                        00000360
        DX(I + 3) = DA*DX(I + 3)                                        00000370
        DX(I + 4) = DA*DX(I + 4)                                        00000380
   50 CONTINUE                                                          00000390
      RETURN                                                            00000400
      END                                                               00000410
