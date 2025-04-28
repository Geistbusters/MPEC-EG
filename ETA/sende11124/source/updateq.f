c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine updateq(vecq,gq,eig,nroot,ncsfq,hdig,ipointq,nfcnq,    3d25s21
     $     iqsbas,ncsf,mdon,isto,savqv,ngot,iprint)                     3d25s21
      implicit real*8 (a-h,o-z)
      dimension vecq(ncsfq,*),gq(ncsfq,*),eig(*),hdig(*),ipointq(*),    3d25s21
     $     iqsbas(3,*),ncsf(*),savqv(ncsfq,*)                           3d25s21
      data iseedx/-15382567/                                            6d20s24
      save                                                              6d20s24
      iqs=1                                                             3d25s21
      if(iprint.ne.0)then                                               3d25s21
       write(6,*)('Hi, my name is updateq')
       write(6,*)('input vecq ')
       call prntm2(vecq,ncsfq,nroot,ncsfq)
       write(6,*)('input g')
       call prntm2(gq,ncsfq,nroot,ncsfq)
       write(6,*)('input eig ')
       call prntm2(eig,1,nroot,1)
      end if                                                            3d25s21
      do if=1,nfcnq                                                     3d25s21
       ihf=ipointq(if)                                                  3d25s21
       huse=hdig(ihf)                                                   3d25s21
       nclo=iqsbas(1,if)                                                3d25s21
       nclop=nclo+1                                                     3d25s21
       iarg=nclop-mdon                                                  3d25s21
       do ir=1,nroot                                                    3d25s21
        bot=huse-eig(ir)                                                3d25s21
        bot=bot+1d-14                                                   3d25s21
        boti=1d0/bot                                                    3d25s21
        do i=0,ncsf(iarg)-1                                              3d25s21
         resid=gq(iqs+i,ir)-eig(ir)*vecq(iqs+i,ir)                      3d25s21
         update=-resid*boti                                             3d25s21
         vecq(iqs+i,ir)=update                                          3d25s21
        end do                                                          3d25s21
       end do                                                           3d25s21
       iqs=iqs+ncsf(iarg)                                               3d25s21
      end do                                                            3d25s21
      if(iprint.ne.0)then                                               3d25s21
       write(6,*)('unnormalized trial')
       call prntm2(vecq,ncsfq,nroot,ncsfq)
      end if                                                            3d25s21
      isto=0                                                            3d25s21
      do i=1,nroot                                                      3d25s21
       itryx=0                                                          6d20s24
 1215  continue                                                         6d20s24
       do j=1,ngot                                                      3d25s21
        dot=0d0                                                         3d25s21
        do k=1,ncsfq                                                    3d25s21
         dot=dot+vecq(k,i)*savqv(k,j)                                   3d25s21
        end do                                                          3d25s21
        do k=1,ncsfq                                                    3d25s21
         vecq(k,i)=vecq(k,i)-dot*savqv(k,j)                             3d25s21
        end do                                                          3d25s21
       end do                                                           3d25s21
       im=i-1                                                           3d25s21
       do j=1,im                                                        3d25s21
        dot=0d0                                                         3d25s21
        do k=1,ncsfq                                                    3d25s21
         dot=dot+vecq(k,i)*vecq(k,j)                                    3d25s21
        end do                                                          3d25s21
        do k=1,ncsfq                                                    3d25s21
         vecq(k,i)=vecq(k,i)-dot*vecq(k,j)                              3d25s21
        end do                                                          3d25s21
       end do                                                           3d25s21
       dot=0d0                                                          3d25s21
       do k=1,ncsfq                                                     3d25s21
        dot=dot+vecq(k,i)**2                                             3d25s21
       end do                                                           3d25s21
c     if this threshold does not reproduce the same number of fcns
c     as in grest, we have a problem!
       if(dot.gt.1d-12)then                                             2d28s22
        isto=isto+1                                                     3d25s21
        doti=1d0/sqrt(dot)                                               3d25s21
        do k=1,ncsfq                                                     3d25s21
         vecq(k,isto)=vecq(k,i)*doti                                    3d25s21
        end do                                                           3d25s21
       end if                                                           3d25s21
      end do                                                            3d25s21
      xsto=dfloat(isto)                                                 4d15s21
      call dws_bcast(xsto,1)                                            4d15s21
      isto=nint(xsto)                                                   4d15s21
      call dws_bcast(vecq,ncsfq*isto)                                   4d15s21
      return                                                            3d25s21
      end                                                               3d25s21
