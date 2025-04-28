c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine ntrans(ivec,wavef,bc,ibc)                              4d25s23
      implicit real*8 (a-h,o-z)                                         4d25s23
c
c     take vectors in ivec, read in vector file information from head   4d25s23
c     of wavef file, and create new file norbs with old orbs transformed4d25s23
c     to new orbs.                                                      4d25s23
c
      character*(*) wavef                                               4d25s23
      dimension ivec(*),nameci(6),isinfo(11),                           4d25s23
     $     nbasdws(8),nlorcont(8),nbasisc(8),nbasisp(8),                4d25s23
     $     idoub(8),iacto(8),icanog(8)                                  5d3s23
      include "common.store"                                            4d25s23
      include "common.basis"                                            4d25s23
      include "common.input"                                            4d25s23
      ibcoffo=ibcoff                                                    4d25s23
      write(6,*)('hi, my name is ntrans for wavef file '),wavef         4d25s23
      open(unit=1,file=wavef,form='unformatted')                        4d25s23
      open(unit=7,file='norbs',form='unformatted')                      4d25s23
      read(1)ismult,isymmrci,nroot,norb,nsb,nameci                      12d27s19
      read(1)idum
      read(1)idum
      read(1)nsymb,idorel,ngaus,natom,nwcont,numeminus,lmax,            4d25s23
     $        nbasallp,multh,ascale,potdws,ipropsym,(isinfo(j),j=1,11)  4d25s23
     $        ,nextradata                                               4d25s23
      write(7)nsymb,idorel,ngaus,natom,nwcont,numeminus,lmax,            4d25s23
     $        nbasallp,multh,ascale,potdws,ipropsym,(isinfo(j),j=1,11)  4d25s23
     $        ,nextradata                                               4d25s23
      read(1)(nbasdws(isb),isb=1,nsymb)
      write(7)(nbasdws(isb),isb=1,nsymb)
      read(1)(nbasisc(isb),isb=1,nsymb)                                 4d25s23
      write(7)(nbasisc(isb),isb=1,nsymb)                                 4d25s23
      read(1)(nbasisp(isb),isb=1,nsymb)                                 4d25s23
      write(7)(nbasisp(isb),isb=1,nsymb)                                 4d25s23
      iext=ibcoff                                                       4d25s23
      ibcoff=iext+nextradata                                            4d25s23
      call enough('ntrans.iext',bc,ibc)                                 4d25s23
      jext=iext-1                                                       4d25s23
      read(1)isym,(xcart(i,1),atnum(i,1),i=1,natom*3),                      4d25s23
     $        (bc(jext+i),i=1,nextradata),ehfbsofarl,                   4d25s23
     $        (idoub(isb),iacto(isb),isb=1,nsymb),nlzz                  1d2s20
      write(7)isym,(xcart(i,1),atnum(i,1),i=1,natom*3),                      4d25s23
     $        (bc(jext+i),i=1,nextradata),ehfbsofarl,                   4d25s23
     $        (idoub(isb),iacto(isb),isb=1,nsymb),nlzz                  1d2s20
      ibcoff=iext
      read(1)((iapair(j,i),j=1,3),i=1,natom)                            4d25s23
      write(7)((iapair(j,i),j=1,3),i=1,natom)                           4d25s23
      nbdat=ngaus*9+nwcont
      ibdat=ibcoff                                                      4d25s23
      ibcoff=ibdat+nbdat                                                4d25s23
      call enough('ntrans.ibdat',bc,ibc)                                4d25s23
      read(1)(bc(ibdat+i),i=0,nbdat-1)                                  4d25s23
      write(7)(bc(ibdat+i),i=0,nbdat-1)                                 4d25s23
      ibcoff=ibdat                                                      4d25s23
      ibstor=ibcoff                                                     4d25s23
      isstor=ibstor+nbasallp*3                                          4d25s23
      ibcoff=isstor+nbasallp*3                                          4d25s23
      call enough('ntrans.ibstor',bc,ibc)                               4d25s23
      jbstor=ibstor-1                                                   4d25s23
      jsstor=isstor-1                                                   4d25s23
      read(1)(ibc(jbstor+i),ibc(jsstor+i),i=1,nbasallp*3)               4d25s23
      write(7)(ibc(jbstor+i),ibc(jsstor+i),i=1,nbasallp*3)              4d25s23
      ibcoff=ibstor                                                     4d25s23
      do isb=1,nsymb                                                    4d25s23
       ibcode=ibcoff                                                    4d25s23
       ibcoff=ibcode+nbasdws(isb)                                       4d25s23
       call enough('ntrans.bcode',bc,ibc)                               4d25s23
       write(7)(ibc(ibcode+i),i=0,nbasdws(isb)-1)                        4d25s23
       ibcoff=ibcode                                                    4d25s23
      end do                                                            4d25s23
      ncomp=1                                                           4d25s23
      if(idorel.ne.0)ncomp=2                                            4d25s23
      do isb=1,nsymb                                                    4d25s23
       nrow=nbasisp(isb)*ncomp                                          4d25s23
       if(nrow.gt.0)then                                                4d25s23
        iorbo=ibcoff                                                    4d25s23
        iorbn=iorbo+nbasdws(isb)*nrow                                   4d25s23
        ibcoff=iorbn+nbasdws(isb)*nrow                                  4d25s23
        call enough('ntrans.orbo',bc,ibc)                               4d25s23
        read(1)(bc(iorbo+i),i=0,nbasdws(isb)*nrow-1)                    4d25s23
        call dgemm('n','n',nrow,nbasdws(isb),nbasdws(isb),1d0,bc(iorbo),4d25s23
     $       nrow,bc(ivec(isb)),nbasdws(isb),0d0,bc(iorbn),nrow,        4d25s23
     $       'ntrans.ao')                                               4d25s23
        write(7)(bc(iorbn+i),i=0,nbasdws(isb)*nrow-1)                   4d25s23
        ibcoff=iorbo                                                    4d25s23
       else                                                             4d25s23
        read(1)                                                         4d25s23
        write(7)                                                        4d25s23
       end if                                                           4d25s23
       if(nbasdws(isb).gt.0)then                                        4d25s23
        morbo=ibcoff                                                    4d25s23
        morbn=morbo+nbasisc(isb)*nbasdws(isb)                           4d25s23
        ibcoff=morbn+nbasisc(isb)*nbasdws(isb)                          4d25s23
        call enough('ntrans.morbo',bc,ibc)                              4d25s23
        read(1)(bc(morbo+i),i=0,nbasisc(isb)*nbasdws(isb)-1)            4d25s23
        call dgemm('n','n',nbasisc(isb),nbasdws(isb),nbasdws(isb),1d0,  4d25s23
     $       bc(morbo),nbasisc(isb),bc(ivec(isb)),nbasdws(isb),0d0,     4d25s23
     $       bc(morbn),nbasisc(isb),'ntrans.morbn')                     4d25s23
        write(7)(bc(morbn+i),i=0,nbasisc(isb)*nbasdws(isb)-1)           4d25s23
        ibcoff=morbo                                                    4d25s23
       else                                                             4d25s23
        read(1)                                                         4d25s23
        write(7)                                                        4d25s23
       end if                                                           4d25s23
      end do                                                            4d25s23
      close(unit=1)                                                     4d25s23
      close(unit=7)                                                     4d25s23
      ibcoff=ibcoffo                                                    4d25s23
      return                                                            4d25s23
      end                                                               4d25s23
