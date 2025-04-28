c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      SUBROUTINE LUSOLV(A,LDA,N,B,LDB,NRHS,IPVT,IERR,IFLAG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     SOLVE LINEAR EQUATIONS USING UNROLLED LOOPS AS DESCRIBED IN
C     1. J.J. DONGARRA AND S.C. EISENSTAT, ACM TRANS. MATH. SOFTWARE
C        VOL. 10, P. 221 (1984)
C     2. J.J. DONGARRA, L. KAUFMAN, AND S. HAMMARLING, ARGONNE NATIONAL
C        LABORATORY, MATHEMATICS AND COMPUTER SCIENCE DIVISION,
C        TECHNICAL MEMORANDUM NO. 46, JAN. 1985, UNPUBLISHED.
C     3. D.W. SCHWENKE AND D.G. TRUHLAR, APPENDIX OF 'CONVERGED
C        CALCULATIONS OF ROTATIONAL EXCITATION AND V-V ENERGY TRANSFER
C        IN THE COLLISION OF TWO MOLECULES', IN "SUPERCOMPUTER
C        SIMULATIONS IN CHEMISTRY", M. DUPUIS, ED., SPRINGER-VERLAG,
C        TO BE PUBLISHED. THIS IS ALSO AVAILABLE AS UNIVERSITY OF
C        MINNESOTA SUPERCOMPUTER INSTITUTE TECHNICAL REPORT UMSI 85/16.
C
C     SUBROUTINES LUX,SMXPYX,SXMPYX ARE FROM J.J. DONGARRA.
C     SXMPYXM IS A MODIFICATION OF SXMPYX, AND LUXS IS FROM REF. 3.
C
C     FOR THE PRESENT ROUTINES, ALL LOOPS ARE UNROLLED TO A LEVEL OF 16.
C
C     SOLVE   AX=B   FOR X.
C
C ON INPUT:
C
C     LDA IS THE FIRST DIMENSION OF A
C
C     N IS THE ORDER OF THE MATRIX A
C
C     LDB IS THE FIRST DIMENSION OF B
C
C     NRHS IS THE NUMBER OF COLUMNS OF B AND HENCE THE NUMBER OF X
C     VECTORS TO BE SOLVED FOR.
C
C     IPVT IS A SCRATCH INTEGER ARRAY OF LENGTH N
C
C     IERR INDICATES AN ERROR CONDITION. IF IERR = 0, NO ERROR HAS
C     OCCURRED. IF IERR NE 0, IERR IS THE NUMBER OF THE PIVOT
C     ELEMENT WHICH IS ZERO.  IF IERR NE 0, THE CONTENTS OF A
C     ARE GARBAGE, WHILE THE CONTENTS OF B ARE UNALTERED.
C
C     IFLAG CONTROLS THE PATH THROUGH ROUTINES. IF IFLAG = 1, ONLY
C     THE LU DECOMPOSITION OF A IS CALCULATED. IF IFLAG = 2,
C     THE LU DECOMPOSITION OF A FROM A PREVIOUS CALL TO LUSOLV
C     IS USED TO DETERMINE THE X MATRIX. IF IFLAG = 3, THE
C     LU DECOMPOSITION OF A IS CALCULATED AND THE MATRIX X IS FOUND.
C     IF IFLAG = 1, LDB,NRHS,B CAN BE DUMMIES.
C     IF IFLAG = 2, A,LDA,N,IPVT MUST BE UNCHANGED FROM THE CALL TO
C     LUSOLV WHICH CALCULATED THE LU DECOMPOSTION OF A.
C
C     FOR OPTIMUM PERFORMANCE, LDA AND LDB SHOULD NOT BE MULTIPLES OF
C     2*(CP/ MCT)*NMB, WHERE CP IS THE CPU CLOCK PERIOD, MCT IS THE
C     MEMORY BANK CYCLE TIME, AND NMB IS THE NUMBER OF MEMORY BANKS.
C     ON THE CRAY-1 WITH NMB = 16, CP/ MCT = 4, THE FIRST DIMENSIONS
C     SHOULD NOT BE MULTIPLES OF 8.
C
C ON RETURN:
C
C     A IS OVERWRITTEN WITH INFORMATION CONCERNING THE LU DECOMPOSITION
C     OF THE ORIGINAL CONTENTS OF A
C
C     IPVT CONTAINS THE PIVOT INDICIES
C
C     B IS OVERWRITTEN WITH THE CALCULATED X MATRIX
C
      DIMENSION A(LDA,N),B(LDB,NRHS),IPVT(N)
      IF(IFLAG.NE.2)CALL LUX(A,LDA,N,IPVT,IERR)
      IF((IERR.NE.0.and.iflag.eq.3).OR.IFLAG.EQ.1)RETURN
      CALL LUXS(A,LDA,N,IPVT,B,LDB,NRHS)
      RETURN
      END
C
C****
C
      SUBROUTINE LUX(A,LDA,N,IPVT,INFO)
C ** LU.F -- LU DECOMPOSITION
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LDA,N,IPVT(N),INFO
      DOUBLE PRECISION A(LDA,1),T
C
      INFO = 0
      DO 40 J = 1, N
C
C       FORM J-TH COLUMN OF L
C
         CALL SMXPYX(N-J+1,A(J,J),J-1,LDA,A(1,J),A(J,1))
C
C       SEARCH FOR PIVOT
C
         T = ABS(A(J,J))
         K = J
         DO 10 I = J+1, N
            IF (ABS(A(I,J)) .GT. T) THEN
               T = ABS(A(I,J))
               K = I
            END IF
   10    CONTINUE
         IPVT(J) = K
C
C       TEST FOR ZERO PIVOT
C
         IF (T .EQ. 0) THEN
            INFO = J
            GO TO 50
         ENDIF
C
C       INTERCHANGE ROWS
C
         DO 20 I = 1, N
            T = A(J,I)
            A(J,I) = A(K,I)
            A(K,I) = T
   20    CONTINUE
C
C       FORM J-TH ROW OF U
C
         A(J,J) = 1/A(J,J)
         CALL SXMPYX(N-J,LDA,A(J,J+1),J-1,LDA,A(J,1),LDA,A(1,J+1))
         T = -A(J,J)
         DO 30 I = J+1, N
            A(J,I) = T*A(J,I)
   30    CONTINUE
   40 CONTINUE
   50 RETURN
      END
C
C***
C
      SUBROUTINE LUXS(A,LDA,N,IPVT,B,LDB,NRHS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     SOLVE AX=B GIVEN LU DECOMPOSITION OF A.
C     SEE REF. 2 IN MAIN PROGRAM (LUSOLV) FOR A DESCRIPTION OF THE
C     ALGORITHM.
C
      DIMENSION A(LDA,N),B(LDB,NRHS),IPVT(N)
C
C      SWAP ROWS
C
      DO 20 K=1,N
      L=IPVT(K)
      IF(L.EQ.K)GOTO 20
      DO 10 J=1,NRHS
      T=B(L,J)
      B(L,J)=B(K,J)
   10 B(K,J)=T
   20 CONTINUE
C
C     FORWARD SUBSTITUTION
C
      DO 40 J=1,N
      CALL SXMPYXM(NRHS,LDB,B(J,1),J-1,LDA,A(J,1),LDB,B)
      DO 120 K=1,NRHS
  120 B(J,K)=B(J,K)*A(J,J)
   40 CONTINUE
C
C     BACKWARD ELIMINATION
C
      DO 60 I=N-1,1,-1                                                  1d11s23
      CALL SXMPYX(NRHS,LDB,B(I,1),N-I,LDA,A(I,I+1),LDB,B(I+1,1))
   60 CONTINUE
C
      RETURN
      END
C
C***
C
      SUBROUTINE SMXPYX(N1,Y,N2,LDM,X,M)
C ** SMXPY16.F -- SMXPY UNROLLED TO A DEPTH OF 16
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LDM,N1,N2
      DOUBLE PRECISION Y(*),X(*),M(LDM,*)
C
C   CLEANUP ODD VECTOR
C
      J = MOD(N2,2)
      IF (J .GE. 1) THEN
         DO 10 I = 1, N1
            Y(I) = (Y(I)) + X(J)*M(I,J)
   10    CONTINUE
      ENDIF
C
C   CLEANUP ODD GROUP OF TWO VECTORS
C
      J = MOD(N2,4)
      IF (J .GE. 2) THEN
         DO 20 I = 1, N1
            Y(I) = ( (Y(I))
     $             + X(J-1)*M(I,J-1)) + X(J)*M(I,J)
   20    CONTINUE
      ENDIF
C
C   CLEANUP ODD GROUP OF FOUR VECTORS
C
      J = MOD(N2,8)
      IF (J .GE. 4) THEN
         DO 30 I = 1, N1
            Y(I) = ((( (Y(I))
     $             + X(J-3)*M(I,J-3)) + X(J-2)*M(I,J-2))
     $             + X(J-1)*M(I,J-1)) + X(J)  *M(I,J)
   30    CONTINUE
      ENDIF
C
C   CLEANUP ODD GROUP OF EIGHT VECTORS
C
      J = MOD(N2,16)
      IF (J .GE. 8) THEN
         DO 40 I = 1, N1
            Y(I) = ((((((( (Y(I))
     $             + X(J-7)*M(I,J-7)) + X(J-6)*M(I,J-6))
     $             + X(J-5)*M(I,J-5)) + X(J-4)*M(I,J-4))
     $             + X(J-3)*M(I,J-3)) + X(J-2)*M(I,J-2))
     $             + X(J-1)*M(I,J-1)) + X(J)  *M(I,J)
   40    CONTINUE
      ENDIF
C
C   MAIN LOOP - GROUPS OF SIXTEEN VECTORS
C
      JMIN = J+16
      DO 60 J = JMIN, N2, 16
         DO 50 I = 1, N1
            Y(I) = ((((((((((((((( (Y(I))
     $             + X(J-15)*M(I,J-15)) + X(J-14)*M(I,J-14))
     $             + X(J-13)*M(I,J-13)) + X(J-12)*M(I,J-12))
     $             + X(J-11)*M(I,J-11)) + X(J-10)*M(I,J-10))
     $             + X(J- 9)*M(I,J- 9)) + X(J- 8)*M(I,J- 8))
     $             + X(J- 7)*M(I,J- 7)) + X(J- 6)*M(I,J- 6))
     $             + X(J- 5)*M(I,J- 5)) + X(J- 4)*M(I,J- 4))
     $             + X(J- 3)*M(I,J- 3)) + X(J- 2)*M(I,J- 2))
     $             + X(J- 1)*M(I,J- 1)) + X(J)   *M(I,J)
   50    CONTINUE
   60 CONTINUE
      RETURN
      END
C
C***
C
      SUBROUTINE SXMPYX(N1,LDY,Y,N2,LDX,X,LDM,M)
C ** SXMPY16.F -- SXMPY UNROLLED TO A DEPTH OF 16
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N1,LDY,N2,LDX,LDM
      DOUBLE PRECISION Y(LDY,*),X(LDX,*),M(LDM,*)
C
C   CLEANUP ODD VECTOR
C
      J = MOD(N2,2)
      IF (J .GE. 1) THEN
         DO 10 I = 1, N1
            Y(1,I) = (Y(1,I)) + X(1,J)*M(J,I)
   10    CONTINUE
      ENDIF
C
C   CLEANUP ODD GROUP OF TWO VECTORS
C
      J = MOD(N2,4)
      IF (J .GE. 2) THEN
         DO 20 I = 1, N1
            Y(1,I) = ( (Y(1,I))
     $               + X(1,J-1)*M(J-1,I)) + X(1,J)*M(J,I)
   20    CONTINUE
      ENDIF
C
C   CLEANUP ODD GROUP OF FOUR VECTORS
C
      J = MOD(N2,8)
      IF (J .GE. 4) THEN
         DO 30 I = 1, N1
            Y(1,I) = ((( (Y(1,I))
     $               + X(1,J-3)*M(J-3,I)) + X(1,J-2)*M(J-2,I))
     $               + X(1,J-1)*M(J-1,I)) + X(1,J)  *M(J,I)
   30    CONTINUE
      ENDIF
C
C   CLEANUP ODD GROUP OF EIGHT VECTORS
C
      J = MOD(N2,16)
      IF (J .GE. 8) THEN
         DO 40 I = 1, N1
            Y(1,I) = ((((((( (Y(1,I))
     $               + X(1,J-7)*M(J-7,I)) + X(1,J-6)*M(J-6,I))
     $               + X(1,J-5)*M(J-5,I)) + X(1,J-4)*M(J-4,I))
     $               + X(1,J-3)*M(J-3,I)) + X(1,J-2)*M(J-2,I))
     $               + X(1,J-1)*M(J-1,I)) + X(1,J)  *M(J,I)
   40    CONTINUE
      ENDIF
C
C   MAIN LOOP - GROUPS OF SIXTEEN VECTORS
C
      JMIN = J+16
      DO 60 J = JMIN, N2, 16
         DO 50 I = 1, N1
            Y(1,I) = ((((((((((((((( (Y(1,I))
     $               + X(1,J-15)*M(J-15,I)) + X(1,J-14)*M(J-14,I))
     $               + X(1,J-13)*M(J-13,I)) + X(1,J-12)*M(J-12,I))
     $               + X(1,J-11)*M(J-11,I)) + X(1,J-10)*M(J-10,I))
     $               + X(1,J- 9)*M(J- 9,I)) + X(1,J- 8)*M(J- 8,I))
     $               + X(1,J- 7)*M(J- 7,I)) + X(1,J- 6)*M(J- 6,I))
     $               + X(1,J- 5)*M(J- 5,I)) + X(1,J- 4)*M(J- 4,I))
     $               + X(1,J- 3)*M(J- 3,I)) + X(1,J- 2)*M(J- 2,I))
     $               + X(1,J- 1)*M(J- 1,I)) + X(1,J)   *M(J,I)
   50    CONTINUE
   60 CONTINUE
      RETURN
      END
C
C***
C
      SUBROUTINE SXMPYXM(N1,LDY,Y,N2,LDX,X,LDM,M)
C ** SXMPYM16.F -- SXMPYM UNROLLED TO A DEPTH OF 16
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N1,LDY,N2,LDX,LDM
      DOUBLE PRECISION Y(LDY,*),X(LDX,*),M(LDM,*)
C
C   CLEANUP ODD VECTOR
C
      J = MOD(N2,2)
      IF (J .GE. 1) THEN
         DO 10 I = 1, N1
            Y(1,I) = (Y(1,I)) - X(1,J)*M(J,I)
   10    CONTINUE
      ENDIF
C
C   CLEANUP ODD GROUP OF TWO VECTORS
C
      J = MOD(N2,4)
      IF (J .GE. 2) THEN
         DO 20 I = 1, N1
            Y(1,I) = ( (Y(1,I))
     $               - X(1,J-1)*M(J-1,I)) - X(1,J)*M(J,I)
   20    CONTINUE
      ENDIF
C
C   CLEANUP ODD GROUP OF FOUR VECTORS
C
      J = MOD(N2,8)
      IF (J .GE. 4) THEN
         DO 30 I = 1, N1
            Y(1,I) = ((( (Y(1,I))
     $               - X(1,J-3)*M(J-3,I)) - X(1,J-2)*M(J-2,I))
     $               - X(1,J-1)*M(J-1,I)) - X(1,J)  *M(J,I)
   30    CONTINUE
      ENDIF
C
C   CLEANUP ODD GROUP OF EIGHT VECTORS
C
      J = MOD(N2,16)
      IF (J .GE. 8) THEN
         DO 40 I = 1, N1
            Y(1,I) = ((((((( (Y(1,I))
     $               - X(1,J-7)*M(J-7,I)) - X(1,J-6)*M(J-6,I))
     $               - X(1,J-5)*M(J-5,I)) - X(1,J-4)*M(J-4,I))
     $               - X(1,J-3)*M(J-3,I)) - X(1,J-2)*M(J-2,I))
     $               - X(1,J-1)*M(J-1,I)) - X(1,J)  *M(J,I)
   40    CONTINUE
      ENDIF
C
C   MAIN LOOP - GROUPS OF SIXTEEN VECTORS
C
      JMIN = J+16
      DO 60 J = JMIN, N2, 16
         DO 50 I = 1, N1
            Y(1,I) = ((((((((((((((( (Y(1,I))
     $               - X(1,J-15)*M(J-15,I)) - X(1,J-14)*M(J-14,I))
     $               - X(1,J-13)*M(J-13,I)) - X(1,J-12)*M(J-12,I))
     $               - X(1,J-11)*M(J-11,I)) - X(1,J-10)*M(J-10,I))
     $               - X(1,J- 9)*M(J- 9,I)) - X(1,J- 8)*M(J- 8,I))
     $               - X(1,J- 7)*M(J- 7,I)) - X(1,J- 6)*M(J- 6,I))
     $               - X(1,J- 5)*M(J- 5,I)) - X(1,J- 4)*M(J- 4,I))
     $               - X(1,J- 3)*M(J- 3,I)) - X(1,J- 2)*M(J- 2,I))
     $               - X(1,J- 1)*M(J- 1,I)) - X(1,J)   *M(J,I)
   50    CONTINUE
   60 CONTINUE
      RETURN
      END
