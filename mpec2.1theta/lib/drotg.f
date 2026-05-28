      SUBROUTINE DROTG(DA,DB,C,S)                                       00000010
C                                                                       00000020
C     CONSTRUCT GIVENS PLANE ROTATION.                                  00000030
C     JACK DONGARRA, LINPACK, 3/11/78.                                  00000040
C                                                                       00000050
      DOUBLE PRECISION DA,DB,C,S,ROE,SCALE,R,Z                          00000060
C                                                                       00000070
      ROE = DB                                                          00000080
      IF( DABS(DA) .GT. DABS(DB) ) ROE = DA                             00000090
      SCALE = DABS(DA) + DABS(DB)                                       00000100
      IF( SCALE .NE. 0.0D0 ) GO TO 10                                   00000110
         C = 1.0D0                                                      00000120
         S = 0.0D0                                                      00000130
         R = 0.0D0                                                      00000140
         GO TO 20                                                       00000150
   10 R = SCALE*DSQRT((DA/SCALE)**2 + (DB/SCALE)**2)                    00000160
      R = DSIGN(1.0D0,ROE)*R                                            00000170
      C = DA/R                                                          00000180
      S = DB/R                                                          00000190
   20 Z = 1.0D0                                                         00000200
      IF( DABS(DA) .GT. DABS(DB) ) Z = S                                00000210
      IF( DABS(DB) .GE. DABS(DA) .AND. C .NE. 0.0D0 ) Z = 1.0D0/C       00000220
      DA = R                                                            00000230
      DB = Z                                                            00000240
      RETURN                                                            00000250
      END                                                               00000260
