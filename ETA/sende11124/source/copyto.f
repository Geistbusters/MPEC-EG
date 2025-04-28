c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine copyto(iorb,iorbao,ifrom,ito,bc,ibc)                   11d9s22
      implicit real*8 (a-h,o-z)                                         2d9s20
      include "common.store"                                            2d9s20
      include "common.basis"                                            2d9s20
      include "common.input"                                            2d9s20
c
c     copy orbital file from unit ifrom to unit ito, except replace     2d9s20
c     orbitals with what is under iorb and iorbao.                      2d9s20
c
      logical ldebug                                                    4d15s24
      dimension iorb(*),iorbao(*),isinfo(11),                           2d9s20
     $     nbasdws(8),nbasisc(8),nbasisp(8),idoub(8),iacto(8)           2d9s20
      include "common.print"                                            4d15s24
      ldebug=iprtr(32).ne.0                                             4d15s24
      read(ifrom)nsymb,idorel,ngaus,natom,nwcont,numeminus,lmax,        2d9s20
     $     nbasallp,multh,ascale,potdws,ipropsym,isinfo,nextradata      2d9s20
      write(ito)nsymb,idorel,ngaus,natom,nwcont,numeminus,lmax,         2d9s20
     $     nbasallp,multh,ascale,potdws,ipropsym,isinfo,nextradata      2d9s20
      read(ifrom)(nbasdws(isb),isb=1,nsymb)                             2d9s20
      write(ito)(nbasdws(isb),isb=1,nsymb)                              2d9s20
      read(ifrom)(nbasisc(isb),isb=1,nsymb)                             2d9s20
      write(ito)(nbasisc(isb),isb=1,nsymb)                              2d9s20
      read(ifrom)(nbasisp(isb),isb=1,nsymb)                             2d9s20
      write(ito)(nbasisp(isb),isb=1,nsymb)                              2d9s20
      nneed=natom*6+nextradata+1                                        2d9s20
      read(ifrom)isym,(bc(ibcoff+i),i=0,nneed-1),(idoub(isb),iacto(isb),2d9s20
     $     isb=1,nsymb),nlzz                                            2d9s20
      write(6,*)('last word of extra data '),bc(ibcoff+nneed-1),
     $     nneed
      write(ito)isym,(bc(ibcoff+i),i=0,nneed-1),(idoub(isb),iacto(isb), 2d9s20
     $     isb=1,nsymb),nlzz                                            2d9s20
      read(ifrom)((iapair(j,i),j=1,3),i=1,natom)                        2d9s20
      write(ito)((iapair(j,i),j=1,3),i=1,natom)                         2d9s20
      nbdat=ngaus*9+nwcont                                              2d9s20
      read(ifrom)(bc(ibcoff+i),i=0,nbdat-1)                             2d9s20
      write(ito)(bc(ibcoff+i),i=0,nbdat-1)                              2d9s20
      ibstor=ibcoff-1                                                   2d9s20
      isstor=ibstor+nbasallp*3                                          2d9s20
      read(ifrom)(ibc(ibstor+i),ibc(isstor+i),i=1,nbasallp*3)           2d9s20
      write(ito)(ibc(ibstor+i),ibc(isstor+i),i=1,nbasallp*3)            2d9s20
      do isb=1,nsymb                                                    2d9s20
       read(ifrom)(ibc(ibcoff+i),i=0,nbasdws(isb)-1)                    2d9s20
       write(ito)(ibc(ibcoff+i),i=0,nbasdws(isb)-1)                     2d9s20
      end do                                                            2d9s20
      if(idorel.eq.0)then                                               2d9s20
       ncomp=1                                                          2d9s20
      else                                                              2d9s20
       ncomp=2                                                          2d9s20
      end if                                                            2d9s20
      do isb=1,nsymb                                                    2d9s20
       nrow=nbasisp(isb)*ncomp                                          2d9s20
       if(nrow.gt.0)then                                                2d9s20
        read(ifrom)(bc(ibcoff+i),i=0,nbasdws(isb)*nrow-1)               2d9s20
        if(ldebug)then                                                  4d15s24
         write(6,*)('saving ao orbitals for symmetry '),isb             4d15s24
         call prntm2(bc(iorbao(isb)),nrow,nbasisp(isb),nrow)            4d15s24
         write(6,*)('in orthogonalized basis ')                         4d15s24
         call prntm2(bc(iorb(isb)),nbasc(isb),nbasdws(isb),             4d15s24
     $        nbasisc(isb))                                             4d15s24
        end if                                                          4d15s24
        write(ito)(bc(iorbao(isb)+i),i=0,nbasdws(isb)*nrow-1)           2d9s20
        read(ifrom)(bc(ibcoff+i),i=0,nbasdws(isb)*nbasc(isb)-1)         2d9s20
        write(ito)(bc(iorb(isb)+i),i=0,nbasdws(isb)*nbasc(isb)-1)       2d9s20
       else                                                             2d9s20
        read(ifrom)                                                     2d9s20
        write(ito)                                                      2d9s20
        read(ifrom)                                                     2d9s20
        write(ito)                                                      2d9s20
       end if                                                           2d9s20
      end do                                                            2d9s20
      read(ifrom)nstate                                                 2d9s20
      write(ito)nstate                                                  12d2s22
      do i=1,nstate                                                     2d9s20
       read(ifrom)nhere                                                 2d9s20
       write(ito)nhere                                                  2d9s20
       read(ifrom)(bc(ibcoff+j),j=0,nhere-1)                            2d9s20
       write(ito)(bc(ibcoff+j),j=0,nhere-1)                             2d9s20
      end do                                                            2d9s20
      close(unit=ifrom)                                                 2d9s20
      close(unit=ito)                                                   2d9s20
      return                                                            2d9s20
      end                                                               2d9s20
