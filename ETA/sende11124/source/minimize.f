c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      SUBROUTINE MINIMIZE (X,FX,STEP,IDONE,thrs,thrsv,lpowell,ldeb)     3d23s23
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     MINIMIZE A FUNCTION OF 1 VARIABLE.
C     ON INPUT:
C     X IS A GUESS FOR THE MINIMUM
C     FX IS THE FUNCTION TO BE MINIMIZED EVALUATED AT X
C     STEP IS A STEPSIZE FOR INITIAL SEARCHED FOR THE MINIMUM.
C     IDONE IS 1 ON FIRST CALL TO MINIMIZE.
C     ON RETURN:
C     X IS A NEW ESTIMATE FOR THE MINIMUM
C     IDONE INDICATES CONVERGENCE:
C           IDONE=0 MEANS CONTINUE
C           IDONE=1 MEANS SEARCH HAS COMPLETED.
C
      DIMENSION XS(3),FS(3),A(3,3),B(3),IPVT(3)
      SAVE ISTEP,IDATA,SIG,XS,FS
      logical ldeb,lpowell                                              3d23s23
C      PRINT*,'X,FX ',X,FX
      IF(IDONE.EQ.1)THEN
C
C     FIRST CALL. SAVE X AND FX AND TAKE A STEP
C
       if(ldeb)write(6,*)('}set up in minimize ...'),lpowell                    3d14s23
       XS(1)=X
       FS(1)=FX
       SIG=1D0
       X=X+STEP
       ISTEP=1
       IDATA=1
       IDONE=0
       RETURN
      END IF
      IDONE=0
C
C     NOT FIRST CALL. ARE WE PERFORMING QUADRATIC INTERPOLATION?
C
      IF(ISTEP.EQ.2)THEN
       FMIN=DMIN1(FS(1),FS(2),FS(3),FX)
       IF(X.GT.XS(2))THEN
        IF(FMIN.EQ.FS(2))THEN
         XS(3)=X
         FS(3)=FX
        ELSE
         XS(1)=XS(2)
         FS(1)=FS(2)
         XS(2)=X
         FS(2)=FX
        END IF
       ELSE
        IF(FMIN.EQ.FX)THEN
         XS(3)=XS(2)
         FS(3)=FS(2)
         XS(2)=X
         FS(2)=FX
        ELSE
         XS(1)=X
         FS(1)=FX
        END IF
       END IF
      END IF
    4 CONTINUE
      IF(ISTEP.EQ.2)THEN
C       PRINT*,'XS ',XS
C       PRINT*,'FS ',FS
       avg=(xs(1)+xs(2)+xs(3))/3d0                                      3d13s23
       DO 2 I=1,3
        A(I,1)=1D0
        delta=xs(i)-avg                                                 3d13s23
        A(I,2)=delta
        A(I,3)=delta**2
        B(I)=FS(I)
    2  CONTINUE
       if(ldeb)then                                                     3d14s23
        write(6,*)('}fitting ')
        call prntm2(fs,1,3,1)
        write(6,*)('}at points ')
        call prntm2(xs,1,3,1)
       end if                                                           3d14s23
       CALL LUSOLV(A,3,3,B,3,1,IPVT,IERR,3)
       if(ldeb)then                                                     3d14s23
        write(6,*)('}coefficients ')
        call prntm2(b,1,3,1)
       end if                                                           3d14s23
       IF(IERR.NE.0)THEN
        WRITE(6,3)IERR
    3   FORMAT(/1X,41HIN MINIMIZE, ON RETURN FROM LUSOLV, IERR=,I5)
        STOP
       END IF
       X=avg-0.5D0*B(2)/B(3)                                            3d13s23
       val=b(1)+(x-avg)*(b(2)+(x-avg)*b(3))                             3d14s23
       if(ldeb)then                                                     3d14s23
        write(6,*)('}predicted minimum: '),x
        write(6,*)('predicted value: '),val                             3d23s23
        nnegp=0                                                         3d23s23
        do i=1,3                                                        3d23s23
         diffp=fs(i)-val                                                3d23s23
         write(6,*)('fs('),i,(')-pre = '),diffp                         3d23s23
         if(i.eq.2.and.diffp.lt.0d0)nnegp=nnegp+1                       3d23s23
        end do                                                          3d23s23
        write(6,*)('nnegp = '),nnegp
       end if                                                           3d14s23
       IF(abs(X).LT.1D-5)THEN
        BOT=1D0
       ELSE
        BOT=X
       END IF
C       WRITE(6,*)('INTERPOLATED VALUE '),X,BOT
       if(ldeb)then                                                     3d14s23
        write(6,*)('test3: '),abs((xs(3)-x)/bot),abs(fs(1)-val)
        write(6,*)('test2: '),abs((xs(2)-x)/bot),abs(fs(2)-val)
        write(6,*)('test1: '),abs((xs(1)-x)/bot),abs(fs(3)-val)
       end if
       IF(ABS((XS(3)-X)/BOT).LT.thrs.or.abs(fs(3)-val).lt.thrsv)IDONE=1 3d14s23
       IF(ABS((XS(2)-X)/BOT).LT.thrs.or.abs(fs(2)-val).lt.thrsv)IDONE=1 3d14s23
       IF(ABS((XS(1)-X)/BOT).LT.thrs.or.abs(fs(1)-val).lt.thrsv)IDONE=1 3d14s23
       RETURN
      END IF
C
C     WE ARE STEPPING. IF THIS IS THE FIRST STEP, ARE WE GOING IN
C     THE CORRECT DIRECTION?
C
      IF(IDATA.EQ.1)THEN
       if(ldeb)write(6,*)('}stepping ')                                 3d14s23
        IF(FX.GT.FS(1))THEN
         SIG=-1D0
         XS(2)=X
         FS(2)=FX
         IDATA=2
         X=XS(1)-STEP
         RETURN
        END IF
       if(ldeb)write(6,*)('energy went down, so continue marching')
       XS(2)=X
       FS(2)=FX
       IDATA=2
       X=X+STEP
       RETURN
      END IF
C
C     PUT NEW POINT INTO ARRAY
C
      IDATA=IDATA+1
      IF(SIG.EQ.1D0)THEN
        IF(IDATA.GT.3)THEN
         DO 10 I=1,2
          XS(I)=XS(I+1)
          FS(I)=FS(I+1)
   10    CONTINUE
        END IF
        XS(3)=X
        FS(3)=FX
        X=X+STEP
      ELSE
       DO 11 I=1,2
        XS(4-I)=XS(3-I)
        FS(4-I)=FS(3-I)
   11  CONTINUE
       XS(1)=X
       FS(1)=FX
       X=X-STEP
      END IF
      IF(IDATA.GE.3)THEN
       step=step*1.4d0
C      PRINT*,'X ',XS
C      PRINT*,'F ',FS
C
C     CHECK TO SEE IF WE CAN SWITCH OVER TO QUADRATIC INTERPOLATION
C
       if(ldeb)then                                                     3d14s23
        write(6,*)('}data so far ')
        call prntm2(fs,1,3,1)
       end if                                                           3d14s23
       IF((FS(1).GT.FS(3).AND.FS(3).GT.FS(2)).OR.
     $    (FS(3).GT.FS(1).AND.FS(1).GT.FS(2)))THEN
        if(ldeb)write(6,*)('}we''ve bracketed the minimum! ')           3d14s23
        ISTEP=2
        GO TO 4
       END IF
      END IF
      RETURN
      END
