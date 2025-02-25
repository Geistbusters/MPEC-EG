c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine run(bc,cfile,ehf,i3x,iapairg,ibc,ibdat,ibstor,idenergy,2d14s23
     $     ieigr,ih0ao,ionex,ioooo,ispt1,isstor,isymdws,ivecr,jmats,    2d14s23
     $     kmats,maxdws,medws,mgmat,ncfile,potdws)                      2d14s23
      implicit real*8 (a-h,o-z)                                         2d14s23
      include "common.basis"                                            2d14s23
      include "common.input"                                            2d14s23
      include "common.mrci"
      include "common.store"                                            2d14s23
      write(6,*)('Hi, my name is run ')
      potn=0d0
      do ia=1,natom
       do ib=ia+1,natom                                                 2d14s23
        dist=sqrt((xcart(1,ia)-xcart(1,ib))**2                          2d14s23
     $           +(xcart(2,ia)-xcart(2,ib))**2                          2d14s23
     $           +(xcart(3,ia)-xcart(3,ib))**2)                         2d14s23
        potn=potn+((atnum(1,ia)*atnum(1,ib))/dist)                      2d14s23
       end do                                                           2d14s23
      end do
      write(6,*)('computed potn: '),potn,potdws-potn
      ngaus2=ngaus*2
      ngaus3=ngaus2+ngaus
      ngaus4=ngaus3+ngaus
      ngaus5=ngaus4+ngaus
      ngaus6=ngaus5+ngaus
      ngaus7=ngaus6+ngaus
      jbdat=ibdat-1
      rms=0d0
      do ig=1,ngaus
       iad=jbdat+ig                                                     2d14s23
       l=ibc(iad)
       iad=iad+ngaus
       xp=bc(iad)
       write(6,*)('l,xp: '),l,xp
       iad=iad+ngaus
       xnorm=bc(iad)
       iad=iad+ngaus
       ioff=ibc(iad)
       iad=iad+ngaus
       iadx=iad+ngaus3
       nhere=ibc(iadx)
       write(6,*)('atom no.: '),nhere
       bc(iad)=xcart(1,nhere)
       iad=iad+ngaus
       bc(iad)=xcart(2,nhere)
       iad=iad+ngaus
       bc(iad)=xcart(3,nhere)
      end do
      return
      end
