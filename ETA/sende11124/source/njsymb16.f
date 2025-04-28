c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      real*8 FUNCTION F3J(JD1,JD2,JD3,MD1,MD2,MD3,ITWO)
C FROM NBS TECHNICAL NOTE 409
C
C CALCULATE 3J SYMBOL. THIS FUNCTION WORKS WITH BOTH INTEGRAL AND HALF
C INTEGRAL ARGUMENTS. IF ALL INTEGRAL ARGUMENTS, USE ITWO=2 AND JD1,ETC
C EQUAL TO THE ARGUMENTS. IF SOME HALF INTEGRAL ARGUMENTS, USE ITWO=1
C AND JD1,ETC EQUAL TO TWICE THE ARGUMENTS.
C
      IMPLICIT REAL*16 (A-H,O-Z)
      COMMON/FACT16/FL(922),NCALL                                         9d23s99
      SAVE /FACT16/
      DIMENSION MTRI(9)
      DATA ZERO,HALF,ONE/0.0q+0,0.5q+0,1.0q+0/
      DATA EPS,EPS2,EPS3,EPS4,EPS5,EPS6/1.0q-10,1.0q+30,1.0q-30,8.0q+1,
     $1.0q+10,23.02585092994046q+0/
      J1=JD1*ITWO
      J2=JD2*ITWO
      J3=JD3*ITWO
      M1=MD1*ITWO
      M2=MD2*ITWO
      M3=MD3*ITWO
      IF(NCALL+1867)5,15,5
5     NCALL=-1867
      FL(1)=ZERO
      FL(2)=ZERO
      DO 50 N=3,922                                                     9d23s99
      FN=N-1
50    FL(N)=FL(N-1)+LOG(FN)
   15 continue                                                          11d25s98
      I=J1+J2-J3                                                        11d25s98
      I1=I/2
      IF(I-2*I1)1000,1010,1000
1010  MTRI(1)=I1
      I=J1-J2+J3
      I1=I/2
      IF(I-2*I1)1000,1020,1000
1020  MTRI(2)=I1
      I=-J1+J2+J3
      I1=I/2
      IF(I-2*I1)1000,1030,1000
1030  MTRI(3)=I1
      IF(M1+M2+M3)1000,1040,1000
1040  I=J1+M1
      I1=I/2
      IF(I-2*I1)1000,1050,1000
1050  MTRI(4)=I1
      MTRI(5)=(J1-M1)/2
      I=J2+M2
      I1=I/2
      IF(I-2*I1)1000,1060,1000
1060  MTRI(6)=I1
      MTRI(7)=(J2-M2)/2
      I=J3+M3
      I1=I/2
      IF(I-2*I1)1000,1070,1000
1070  MTRI(8)=I1
      MTRI(9)=(J3-M3)/2
      DO 30 N=1,9
      IF(MTRI(N))1000,30,30
30    CONTINUE
      IF(J3-J2+M1)40,45,45
40    KMIN=-J3+J2-M1
      GOTO 60
45    KMIN=0
60    IF(-J3+J1+M2-KMIN)80,80,70
70    KMIN=-J3+J1+M2
80    KMIN=KMIN/2
      IF(J2-J3+M1)90,100,100
90    KMAX=J1+J2-J3
      GOTO 110
100   KMAX=J1-M1
110   IF(J2+M2-KMAX)120,130,130
120   KMAX=J2+M2
130   KMAX=KMAX/2
      MIN1=MTRI(1)-KMIN+1
      MIN2=MTRI(5)-KMIN+1
      MIN3=MTRI(6)-KMIN+1
      MIN4=(J3-J2+M1)/2+KMIN
      MIN5=(J3-J1-M2)/2+KMIN
      UK=EPS
      S=UK
      NCUT=0
      KMAX=KMAX-KMIN
      IF(KMAX)165,165,155
155   DO 160 K=1,KMAX
       qval1=DFLOAT(MIN1-K)
       qval2=DFLOAT(MIN2-K)
       qval3=DFLOAT(MIN3-K)
       qval4=DFLOAT(KMIN+K)
       qval5=DFLOAT(MIN4+K)
       qval6=DFLOAT(MIN5+K)
      UK=-UK*qval1*qval2*qval3/(qval4*qval5*qval6)
      IF(ABS(UK)-EPS2)158,157,157
157   UK=EPS*UK
      S=EPS*S
      NCUT=NCUT+1
158   IF(ABS(UK)-EPS3)165,160,160
160   S=S+UK
165   DELOG=ZERO
      DO 170 N=1,9
      NUM=MTRI(N)
      if(num+1.gt.922)stop '3jf'                                        9d23s99
170   DELOG=DELOG+FL(NUM+1)
      NUM=(J1+J2+J3)/2+2
      if(num.gt.922)stop '3jf'                                          9d23s99
      DELOG=HALF*(DELOG-FL(NUM))
      if(max(KMIN+1,MIN1,MIN2,MIN3,MIN4+1,MIN5+1).gt.922)stop '3jf'     9d23s99
      ULOG=-FL(KMIN+1)-FL(MIN1)-FL(MIN2)-FL(MIN3)-FL(MIN4+1)-FL(MIN5+1)
      PLOG=DELOG+ULOG
      IF(PLOG+EPS4)172,171,171
171   IF(NCUT)175,175,172
172   SIG=SIGN(ONE,S)
      S=ABS(S)
      qval=DFLOAT(NCUT+1)
      SLOG=LOG(S)+qval*EPS6
      F3J=SIG*EXP(SLOG+PLOG)
      GOTO 178
175   S=S*EPS5
      P=EXP(PLOG)
      F3J=P*S
178   NUM=KMIN+(J1-J2-M3)/2
      IF(MOD(NUM,2))180,190,180
180   F3J=-F3J
190   CONTINUE
      GOTO 2000
1000  F3J=ZERO
2000  RETURN
      END
CC
      real*8 FUNCTION F6J(JD1,JD2,JD3,LD1,LD2,LD3,ITWO)
C  FROM NBS TECHNICAL NOTE 409
C
C  CALCULATE 6J SYMBOL. THIS FUNCTION WORKS WITH BOTH INTEGRAL AND HALF
C  INTEGRAL ARGUMENTS. IF ALL INTEGRAL ARGUMENTS, USE ITWO=2 AND JD1,ETC
C  EQUAL TO THE ARGUMENTS. IF SOME HALF INTEGRAL ARGUMENTS, USE ITWO=1
C  AND JD1,ETC EQUAL TO TWICE THE ARGUMENTS.
C
C  THIS CALLS FUNCTION S6J.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION MED(12)
      DATA ZERO/0.0d+0/
      J1=JD1*ITWO
      J2=JD2*ITWO
      J3=JD3*ITWO
      L1=LD1*ITWO
      L2=LD2*ITWO
      L3=LD3*ITWO
      I=-J1+J2+J3
      I1=I/2
      IF(I-2*I1)1000,1010,1000
1000  F6J=ZERO
      GOTO 100
1010  MED(1)=I1
      I=J1-J2+J3
      I1=I/2
      IF(I-2*I1)1000,1020,1000
1020  MED(2)=I1
      I=J1+J2-J3
      I1=I/2
      IF(I-2*I1)1000,1030,1000
1030  MED(3)=I1
      I=-J1+L2+L3
      I1=I/2
      IF(I-2*I1)1000,1040,1000
1040  MED(4)=I1
      I=J1-L2+L3
      I1=I/2
      IF(I-2*I1)1000,1050,1000
1050  MED(5)=I1
      I=J1+L2-L3
      I1=I/2
      IF(I-2*I1)1000,1060,1000
1060  MED(6)=I1
      I=-L1+J2+L3
      I1=I/2
      IF(I-2*I1)1000,1070,1000
1070  MED(7)=I1
      I=L1-J2+L3
      I1=I/2
      IF(I-2*I1)1000,1080,1000
1080  MED(8)=I1
      I=L1+J2-L3
      I1=I/2
      IF(I-2*I1)1000,1090,1000
1090  MED(9)=I1
      I=-L1+L2+J3
      I1=I/2
      IF(I-2*I1)1000,1100,1000
1100  MED(10)=I1
      I=L1-L2+J3
      I1=I/2
      IF(I-2*I1)1000,1110,1000
1110  MED(11)=I1
      I=L1+L2-J3
      I1=I/2
      IF(I-2*I1)1000,1120,1000
1120  MED(12)=I1
      DO 10 N=1,12
      IF(MED(N))1000,10,10
10    CONTINUE
      F6J=S6J(J1,J2,J3,L1,L2,L3)
100   RETURN
      END
CC
      real*8 FUNCTION S6J(JD1,JD2,JD3,LD1,LD2,LD3)
C FROM NBS TECHNICAL NOTE 409
      IMPLICIT REAL*16 (A-H,O-Z)
      COMMON/FACT/FL(922),NCALL                                         9d23s99
      SAVE /FACT/
      DIMENSION MA(4),MB(3),MED(12)
      DATA ZERO,ONE,HALF/0.0q+0,1.0q+0,0.5q+0/
      DATA EPS1,EPS2,EPS3,EPS4,EPS5/1.0q-15,1.0q-25,1.0q+15,64.0q+0,
     $  1.6038109389511792q-28/
      J1=JD1
      J2=JD2
      J3=JD3
      L1=LD1
      L2=LD2
      L3=LD3
      IF(NCALL+1867)5,15,5
5     NCALL=-1867
      FL(1)=ZERO
      FL(2)=ZERO
      DO 50 N=3,922                                                     9d23s99
      FN=DFLOAT(N-1)
50    FL(N)=FL(N-1)+LOG(FN)
   15 continue                                                          11d25s98
      MED(1)=(-J1+J2+J3)/2                                              11d25s98
      MED(2)=(J1-J2+J3)/2
      MED(3)=(J1+J2-J3)/2
      MED(4)=(-J1+L2+L3)/2
      MED(5)=(J1-L2+L3)/2
      MED(6)=(J1+L2-L3)/2
      MED(7)=(-L1+J2+L3)/2
      MED(8)=(L1-J2+L3)/2
      MED(9)=(L1+J2-L3)/2
      MED(10)=(-L1+L2+J3)/2
      MED(11)=(L1-L2+J3)/2
      MED(12)=(L1+L2-J3)/2
      MA(1)=MED(1)+MED(2)+MED(3)
      MA(2)=MED(4)+MED(5)+MED(6)
      MA(3)=MED(7)+MED(8)+MED(9)
      MA(4)=MED(10)+MED(11)+MED(12)
      MB(1)=MA(1)+MED(12)
      MB(2)=MA(1)+MED(4)
      MB(3)=MA(1)+MED(8)
      MAX=MA(1)
      DO 30 N=2,4
      IF(MAX-MA(N))20,30,30
20    MAX=MA(N)
30    CONTINUE
      MIN=MB(1)
      DO 51 N=2,3
      IF(MIN-MB(N))51,51,40
40    MIN=MB(N)
51    CONTINUE
      KMAX=MIN-MAX
      MINP1=MIN+1
      MIN1=MINP1-MA(1)
      MIN2=MINP1-MA(2)
      MIN3=MINP1-MA(3)
      MIN4=MINP1-MA(4)
      MIN5=MINP1+1
      MIN6=MB(1)-MIN
      MIN7=MB(2)-MIN
      MIN8=MB(3)-MIN
      UK=EPS1
      S=UK
      IF(KMAX)65,65,55
55    DO 60 K=1,KMAX
      UK=-UK*DFLOAT(MIN1-K)*DFLOAT(MIN2-K)*DFLOAT(MIN3-K)*DFLOAT(MIN4-K)
     1/(DFLOAT(MIN5-K)*DFLOAT(MIN6+K)*DFLOAT(MIN7+K)*DFLOAT(MIN8+K))
      IF(ABS(UK)-EPS2)65,65,60
60    S=S+UK
65    S=S*EPS3
      DELOG=ZERO
      DO 70 N=1,12
      NUM=MED(N)
      if(num+1.gt.922)stop 's6jf'                                       9d23s99
70    DELOG=DELOG+FL(NUM+1)
      NUM1=MA(1)+2
      NUM2=MA(2)+2
      NUM3=MA(3)+2
      NUM4=MA(4)+2
      if(max0(NUM1,NUM2,NUM3,NUM4).gt.922)stop 's6jf'                   9d23s99
      DELOG=DELOG-FL(NUM1)-FL(NUM2)-FL(NUM3)-FL(NUM4)
      DELOG=HALF*DELOG
      if(max0(MIN5,MIN1,MIN2,MIN3,MIN4,MIN6+1,MIN7+1,MIN8+1).gt.922)    9d23s99
     $     stop 's6jf'                                                  5d4s98
      ULOG=FL(MIN5)-FL(MIN1)-FL(MIN2)-FL(MIN3)-FL(MIN4)-FL(MIN6+1)-FL
     1(MIN7+1)-FL(MIN8+1)
      PLOG=DELOG+ULOG
      IF(PLOG+EPS4)72,75,75
72    Q=PLOG+EPS4
      Q=EXP(Q)
      S6J=Q*S
      IF(ABS(S6J)-ONE)73,73,74
73    S6J=ZERO
      GOTO 90
74    S6J=S6J*EPS5
      GOTO 78
75    P=EXP(PLOG)
      S6J=P*S
78    MIN2=MIN/2
      IF(MIN-2*MIN2)80,90,80
80    S6J=-S6J
90    CONTINUE
      RETURN
      END
CC
      real*8 FUNCTION F9J(JD1,JD2,JD3,JD4,JD5,JD6,JD7,JD8,JD9,ITWO)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     FROM NBS TECHNICAL NOTE 409.
C
C     CALCULATE 9J SYMBOL. THIS FUNCTION WORKS FOR BOTH INTEGRAL AND
C     HALF INTEGRAL ARGUMENTS. IF ALL INTEGRAL ARGUMENTS, USE ITWO=2 AND
C     JD1,ETC EQUAL TO THE ARGUMENTS. IF SAME HALF INTEGRAL ARGUMENTS,
C     USE ITWO=1 NAD JD1,ETC EQUAL TO TWICE THE ARGUMENTS.
C
C     THIS USES FUNCTION S6J.
C
      DIMENSION MTRIA(18),KN(6),KX(6),NN(6)
      DATA ZERO,TWO/0.0D+0,2.0D+0/
      J1=JD1*ITWO
      J2=JD2*ITWO
      J3=JD3*ITWO
      J4=JD4*ITWO
      J5=JD5*ITWO
      J6=JD6*ITWO
      J7=JD7*ITWO
      J8=JD8*ITWO
      J9=JD9*ITWO
C
      I=-J1+J2+J3
      I1=I/2
      IF(I-2*I1)1000,1010,1000
 1010 MTRIA(1)=I1
      I=J1-J2+J3
      I1=I/2
      IF(I-2*I1)1000,1020,1000
 1020 MTRIA(2)=I1
      I=J1+J2-J3
      I1=I/2
      IF(I-2*I1)1000,1030,1000
 1030 MTRIA(3)=I1
      I=-J4+J5+J6
      I1=I/2
      IF(I-2*I1)1000,1040,1000
 1040 MTRIA(4)=I1
      I=J4-J5+J6
      I1=I/2
      IF(I-2*I1)1000,1050,1000
 1050 MTRIA(5)=I1
      I=J4+J5-J6
      I1=I/2
      IF(I-2*I1)1000,1060,1000
 1060 MTRIA(6)=I1
      I=-J7+J8+J9
      I1=I/2
      IF(I-2*I1)1000,1070,1000
 1070 MTRIA(7)=I1
      I=J7-J8+J9
      I1=I/2
      IF(I-2*I1)1000,1080,1000
 1080 MTRIA(8)=I1
      I=J7+J8-J9
      I1=I/2
      IF(I-2*I1)1000,1090,1000
 1090 MTRIA(9)=I1
      I=-J1+J4+J7
      I1=I/2
      IF(I-2*I1)1000,1100,1000
 1100 MTRIA(10)=I1
      I=J1-J4+J7
      I1=I/2
      IF(I-2*I1)1000,1110,1000
 1110 MTRIA(11)=I1
      I=J1+J4-J7
      I1=I/2
      IF(I-2*I1)1000,1120,1000
 1120 MTRIA(12)=I1
      I=-J2+J5+J8
      I1=I/2
      IF(I-2*I1)1000,1130,1000
 1130 MTRIA(13)=I1
      I=J2-J5+J8
      I1=I/2
      IF(I-2*I1)1000,1140,1000
 1140 MTRIA(14)=I1
      I=J2+J5-J8
      I1=I/2
      IF(I-2*I1)1000,1150,1000
 1150 MTRIA(15)=I1
      I=-J3+J6+J9
      I1=I/2
      IF(I-2*I1)1000,1160,1000
 1160 MTRIA(16)=I1
      I=J3-J6+J9
      I1=I/2
      IF(I-2*I1)1000,1170,1000
 1170 MTRIA(17)=I1
      I=J3+J6-J9
      I1=I/2
      IF(I-2*I1)1000,1180,1000
 1180 MTRIA(18)=I1
C
      DO 30 N=1,18
      IF(MTRIA(N))1000,30,30
   30 CONTINUE
      KN(1)=MAX0(ABS(J2-J6),ABS(J1-J9),ABS(J4-J8))
      KN(2)=MAX0(ABS(J2-J7),ABS(J5-J9),ABS(J4-J3))
      KN(3)=MAX0(ABS(J6-J7),ABS(J5-J1),ABS(J8-J3))
      KN(4)=MAX0(ABS(J6-J1),ABS(J2-J9),ABS(J5-J7))
      KN(5)=MAX0(ABS(J2-J4),ABS(J3-J7),ABS(J6-J8))
      KN(6)=MAX0(ABS(J3-J5),ABS(J1-J8),ABS(J4-J9))
      KX(1)=MIN0(J2+J6,J1+J9,J4+J8)
      KX(2)=MIN0(J2+J7,J5+J9,J4+J3)
      KX(3)=MIN0(J6+J7,J5+J1,J8+J3)
      KX(4)=MIN0(J1+J6,J2+J9,J5+J7)
      KX(5)=MIN0(J2+J4,J3+J7,J6+J8)
      KX(6)=MIN0(J3+J5,J1+J8,J4+J9)
      DO 35 K=1,6
   35 NN(K)=KX(K)-KN(K)
      KSIGN=1
      I=MIN0(NN(1),NN(2),NN(3),NN(4),NN(5),NN(6))
      DO 40 K=1,6
      IF(I-NN(K))40,50,40
   40 CONTINUE
   50 KMIN=KN(K)+1
      KMAX=KX(K)+1
      GOTO (130,52,53,54,55,56),K
   52 J=J1
      J1=J5
      J5=J
      J=J3
      J3=J8
      J8=J
      J=J6
      J6=J7
      J7=J
      GOTO 130
   53 J=J2
      J2=J7
      J7=J
      J=J3
      J3=J4
      J4=J
      J=J5
      J5=J9
      J9=J
      GOTO 130
   54 J=J1
      J1=J2
      J2=J
      J=J4
      J4=J5
      J5=J
      J=J7
      J7=J8
      J8=J
      GOTO 120
   55 J=J1
      J1=J3
      J3=J
      J=J4
      J4=J6
      J6=J
      J=J7
      J7=J9
      J9=J
      GOTO 120
   56 J=J2
      J2=J3
      J3=J
      J=J5
      J5=J6
      J6=J
      J=J8
      J8=J9
      J9=J
  120 KSIGN=(1-MOD(J1+J2+J3+J4+J5+J6+J7+J8+J9,4))
C
  130 SUM=ZERO
      SIG=(-1)**(KMIN-1)*KSIGN
      FLK=dFLOAT(KMIN)
      DO 200 K=KMIN,KMAX,2
      TERM=FLK*S6J(J1,J4,J7,J8,J9,K-1)*S6J(J2,J5,J8,J4,K-1,J6)
     $*S6J(J3,J6,J9,K-1,J1,J2)
      FLK=FLK+TWO
  200 SUM=SUM+TERM
      F9J=SUM*SIG
      GOTO 2000
 1000 F9J=ZERO
 2000 RETURN
      END
