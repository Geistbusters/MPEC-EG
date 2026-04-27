c mpec2.1 version theta. For copyright and Disclaimers, see start.f
      subroutine bodc1fd(iorbb,iorbk,iwavb,iwavk,joder1,wcoef,          2d3s25
     $     natom,ngaus,jbdatb,jbdatk,ibdat,nbasdws,nbasp,multh,isym,    2d14s25
     $     iapair,ibstor,isstor,idorel,ascale,nsymb,bc,ibc)             2d14s25
      implicit real*8 (a-h,o-z)                                         2d3s25
      integer*8 ibcoffo,ixmt,koder1,itransp
      include "common.store"                                            2d3s25
      character*9 name                                                  2d3s25
      dimension iwavb(*),iwavk(*),iorbb(*),iorbk(*),idata(7),nbasdws(*),2d3s25
     $     nbasp(*),multh(8,8),ixmt(8),ibdat(*)                         2d14s25
      write(6,*)('Hi, my name is bodc1fd ')                             2d3s25
      write(6,*)('nbasp: '),(nbasp(isb),isb=1,nsymb)
      write(6,*)('fd weight: '),wcoef
      ibdatb=ibdat(jbdatb)                                              2d14s25
      ibdatk=ibdat(jbdatk)                                              2d14s25
      ibcoffo=ibcoff                                                    2d3s25
      ipt=1
      npt=1
      iosym=1
      name='overlap  '
      i2e=0
      iprt=1
      do i=1,7                                                          5d27s21
       idata(i)=0
      end do
      write(6,*)('stamp1'),idorel,nsymb
      data=+1d0
      nbb=0                                                             2d3s25
      if(idorel.eq.0)then                                               2d3s25
       ncomp=1                                                          2d3s25
      else                                                              2d3s25
       ncomp=2                                                          2d3s25
      end if                                                            2d3s25
      do isb=1,nsymb                                                    2d3s25
       ixmt(isb)=ibcoff                                                 2d3s25
       ibcoff=ixmt(isb)+(nbasp(isb)*ncomp)**2                           2d3s25
       nbb=nbb+(nbasp(isb)*ncomp)**2                                    2d3s25
      end do                                                            2d3s25
      call enough('bodc1fd.xmt',bc,ibc)                                 2d3s25
      do iz=ixmt(1),ibcoff-1                                            2d13s25
       bc(iz)=0d0                                                       2d13s25
      end do                                                            2d13s25
      write(6,*)('stamp2')
      write(6,*)('ibdatb,ibdatk '),ibdatb,ibdatk
      call parap2(natom,ngaus,ibdatb,ibdatk,ixmt,isym,iapair,ibstor,    2d3s25
     $     isstor,iso,0000,idorel,ascale,ipt,npt,data,idata,iosym,1,    2d13s25
     $     multh,nbb,nbasp,nbasdws,iorbb,iorbk,name,1,bc,ibc)           2d14s25
      write(6,*)('stamp3')
      koder1=joder1                                                     2d3s25
      do isb=1,nsymb                                                    2d3s25
       if(nbasdws(isb).gt.0)then                                        2d3s25
        write(6,*)('orbital overlap for sym '),isb
        call prntm2(bc(ixmt(isb)),nbasdws(isb),nbasdws(isb),            2d3s25
     $       nbasdws(isb))                                              2d3s25
        nn=nbasdws(isb)*nbasdws(isb)                                    2d3s25
        do i=0,nn-1                                                     2d3s25
         bc(koder1+i)=bc(koder1+i)+bc(ixmt(isb)+i)*wcoef                2d13s25
        end do                                                          2d3s25
        write(6,*)('fd orb matrix so far '),koder1,loc(bc(koder1))
        call prntm2(bc(koder1),nbasdws(isb),nbasdws(isb),               2d3s25
     $       nbasdws(isb))                                              2d3s25
        itransp=ibcoff
        ibcoff=itransp+nbasdws(isb)*nbasdws(isb)
        do i=0,nbasdws(isb)-1
         do j=0,nbasdws(isb)-1
          ji=koder1+j+nbasdws(isb)*i
          ij=itransp+i+nbasdws(isb)*j
          bc(ij)=bc(ji)
         end do
        end do
        write(6,*)('transposed ')
        call prntm2(bc(itransp),nbasdws(isb),nbasdws(isb),nbasdws(isb))
        ibcoff=itransp
        koder1=koder1+nn                                                2d3s25
       end if                                                           2d3s25
      end do                                                            2d3s25
      ibcoff=ibcoffo                                                    2d3s25
      return                                                            2d3s25
      end                                                               2d3s25
