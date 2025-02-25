c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine moldenf(ibdat,iorbn,nsymb,ngaus,ibstor,isstor,noocc,   11d22s19
     $     idoub,iacto,nbasisp,bc,ibc)                                  11d9s22
      implicit real*8 (a-h,o-z)                                         11d22s19
      include "common.store"                                            11d22s19
      dimension iorbn(nsymb),noocc(nsymb),idoub(nsymb),iacto(nsymb),    11d22s19
     $     nbasisp(nsymb)                                               11d22s19
      write(3,*)('[MO]')                                                11d22s19
      write(3,*)('[5D]')                                                11d22s19
      write(3,*)('[9G]')                                                11d22s19
      do isb=1,nsymb                                                    11d22s19
       tail=dfloat(isb)*0.1d0                                           11d22s19
       iad=iorbn(isb)-1                                                 11d22s19
       do im=1,idoub(isb)                                               11d22s19
        sym=dfloat(im)+tail                                             11d22s19
        write(3,1)sym                                                   11d22s19
    1   format('Sym= ',f5.1)                                            11d22s19
        write(3,*)('Ene= 0.')
        write(3,*)('Spin= Alpha')                                       11d22s19
        write(3,*)('Occup= '),2d0                                       11d22s19
        do ip=1,nbasisp(isb)                                            11d22s19
         write(3,*)ip,bc(iad+ip)                                        11d22s19
        end do                                                          11d22s19
        iad=iad+nbasisp(isb)                                            11d22s19
       end do
       do im=1,iacto(isb)                                               11d22s19
        sym=dfloat(im+idoub(isb))+tail                                  11d22s19
        write(3,1)sym                                                   11d22s19
        write(3,*)('Ene= 0.')
        write(3,*)('Spin= Alpha')                                       11d22s19
        write(3,*)('Occup= '),bc(noocc(isb)+im-1)                       11d22s19
        do ip=1,nbasisp(isb)                                            11d22s19
         write(3,*)ip,bc(iad+ip)                                        11d22s19
        end do                                                          11d22s19
        iad=iad+nbasisp(isb)                                            11d22s19
       end do
      end do                                                            11d22s19
      return                                                            11d22s19
      end                                                               11d22s19
