c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcds12c(nrootu,nll,ip,idd,njhere,nfd,icol,ncol,itmpdc1,9d20s23
     $     itmpdc2,i11s,nnt,ncsf2,ncsf,prod,ldcont,iad2,sumdc,ff,bc,ibc,7d15s24
     $     igqqq,lab)                                                   7d15s24
      implicit real*8 (a-h,o-z)
      integer*8 ip(*)
      character*(*) lab                                                 7d15s24
      include "common.store"
      dimension prod(*)
      do ir=0,nrootu-1
       do iii=1,nll
        iiim=iii-1                                                      9d20s23
        ipp=ip(iii)-1                                                   9d20s23
        iad=idd+njhere*(ipp+nfd*(icol+ncol*ir))                         9d20s23
        jtmpdc1=itmpdc1+iiim+nll*(ir+nrootu*(i11s-1))                   9d20s23
        do j=0,njhere-1                                                 9d20s23
         bc(jtmpdc1)=bc(iad+j)                                          9d20s23
         jtmpdc1=jtmpdc1+nnt                                            9d20s23
        end do                                                          9d20s23
       end do                                                           9d20s23
      end do                                                            9d20s23
      call dgemm('n','n',nnt,ncsf2,ncsf,ff,bc(itmpdc1),nnt,             9d20s23
     $     prod,ncsf,0d0,bc(itmpdc2),nnt,'hcds12c')                     9d20s23
      if(itmpdc2.le.22759755.and.itmpdc2+nnt*ncsf2.gt.22759755)then
      idelta=          22759755-itmpdc2
      ii2=idelta/nnt
      ii1=idelta-nnt*ii2
      zum=0d0
      do j=0,ncsf-1
       iadx1=itmpdc1+ii1+nnt*j
       iadx2=1+j+ncsf*ii2
       term=bc(iadx1)*prod(iadx2)*ff
       zum=zum+term
       if(abs(term).gt.1d-10)write(6,*)j,bc(iadx1),prod(iadx2),zum,
     $      iadx1
      end do
      end if
      iadl=ldcont                                                       9d20s23
      do ir=0,nrootu-1                                                  9d20s23
       do iii=0,nll-1                                                   9d20s23
        do j=0,ncsf2-1                                                  9d20s23
         iadc2=itmpdc2+iii+nll*(ir+nrootu*j)                            9d20s23
         bc(iadl+j)=bc(iadl+j)+bc(iadc2)                                9d20s23
        end do                                                          9d20s23
        iadl=iadl+ncsf2                                                 9d20s23
       end do                                                           9d20s23
      end do                                                            9d20s23
      return
      end
