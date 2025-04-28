c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine graborb(iorb,fname,bc,ibc)                             11d10s22
      implicit real*8 (a-h,o-z)                                         7d26s22
      character*(*) fname                                               7d26s22
      dimension iorb(*),idumx18(18),nbasdws(8),nbasc(8),nbasp(8),       7d27s22
     $     multhx(8,8),isymx(3,8)                                       7d27s22
      include "common.store"                                            7d26s22
      include "common.basis"                                            7d27s22
      write(6,*)('in graborb for '),fname                               7d26s22
       open(unit=1,file=fname,form='unformatted')
       nwavrec=1                                                         5d4s21
       read(1)ismult,isymmrci,nroot,norb,nsb,nameci                      12d27s19
       nwavrec=nwavrec+1                                                 5d4s21
       read(1)idum
       if(norb.gt.0)then                                                 5d4s21
        nwavrec=nwavrec+1                                                 5d4s21
        read(1)idum                                                      5d4s21
       end if                                                            5d4s21
        write(6,*)('spin multiplicity '),ismult
        write(6,*)('symmetry of wavefunction '),isymmrci
        write(6,*)('number of roots = '),nroot
       read(1)nsymb,idorelx,ngaus,natom,nwcont,numeminus,lmax,nbasallp, 8d30s21
     $     multhx,ascale,potdws,idumx18
       write(6,*)('nsymb,idorelx '),nsymb,idorelx
       nextradatatag=idumx18(18)+1                                      2d14s22
       nvec=natom-1                                                      5d26s21
       nexpdata=0                                                        5d26s21
       if(nvec.eq.1)then                                                 5d26s21
        nexpdata=1                                                       5d26s21
       else if(nvec.gt.1)then                                            5d26s21
        nexpdata=3*nvec-3                                                5d26s21
       end if                                                            5d26s21
       read(1)(nbasdws(isb),isb=1,nsymb)                                 11d20s19
       write(6,*)('nbasdws: '),(nbasdws(isb),isb=1,nsymb)
       read(1)(nbasc(isb),isb=1,nsymb)                                   11d20s19
       write(6,*)('nbasc: '),(nbasc(isb),isb=1,nsymb)
       read(1)(nbasp(isb),isb=1,nsymb)                                   11d20s19
       write(6,*)('nbasp: '),(nbasp(isb),isb=1,nsymb)
       itmp1=ibcoff
       itmp2=itmp1+natom*3
       iextrad=itmp2+natom*3
       ibcoff=iextrad+nextradatatag                                      5d27s21
       call enough('graborb.  1',bc,ibc)
       nwavrec=nwavrec+1                                                 5d4s21
       read(1)isymx,(bc(itmp1+i),bc(itmp2+i),i=0,natom*3-1),             5d27s21
     $     (bc(iextrad+i),i=0,nextradatatag-1)                          5d27s21
       if(nextradatatag-nexpdata.ge.6+natom)then                        7d27s22
        write(6,*)
     $     ('it looks like we have cmx and iagrp data tagging along')
        jextrad=iextrad+nexpdata                                         5d27s21
        jextrad=ibcoff-7-natom
        do i=1,2                                                         5d27s21
         do j=1,3                                                        5d27s21
          cmx(j,i)=bc(jextrad)                                           5d27s21
          jextrad=jextrad+1                                              5d27s21
         end do                                                          5d27s21
        end do                                                           5d27s21
        do i=1,natom                                                     5d27s21
         iagrp(i)=ibc(jextrad)                                           5d27s21
         jextrad=jextrad+1                                               5d27s21
        end do                                                           5d27s21
       end if                                                            5d26s21
       ibcoff=iextrad                                                    5d27s21
       read(1)dum
       read(1)dum
       read(1)idum
       if(idorelx.eq.0)then                                             7d26s22
        ncomp=1                                                         7d26s22
       else                                                             7d26s22
        ncomp=2                                                         7d26s22
       end if                                                           7d26s22
       do isb=1,nsymb                                                   7d26s22
        need=nbasp(isb)*nbasdws(isb)*ncomp                               1d11s20
        iorb(isb)=ibcoff                                                7d26s22
        ibcoff=iorb(isb)+need                                           7d26s22
        if(need.gt.0)then
         read(1)(bc(iorb(isb)+i),i=0,need-1)                            7d26s22
         read(1)dum
        end if
       end do                                                           7d26s22
      close(unit=1)                                                     7d26s22
      return                                                            7d26s22
      end                                                               7d26s22
