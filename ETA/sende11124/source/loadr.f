c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine loadr(potdws,ibdat,ibstor,isstor,nsymb,istinfo,        1d2s20
     $     idoub,iacto,nlzzq,iextradata,nextradata,iflag,nbasdws,       2d24s21
     $     element,nlorcont,nec,mdon,mdoop,iflagr,bc,ibc,icanog)        5d5s23
      implicit real*8 (a-h,o-z)
c
c     we are getting our basic input from forbs file rather than from   11d20s19
c     input file ...                                                    11d20s19
c
      character*3 element(*)                                            2d24s21
      dimension nbasdws(8),istinfo(11),idoub(8),iacto(8),nlorcont(*),   5d3s23
     $     icanog(8)                                                    5d3s23
      integer*2 ipack2(2)                                               5d4s23
      equivalence (ipack2,ipack4)                                       5d4s23
      include "common.basis"                                            11d20s19
      include "common.input"                                            11d20s19
      include "common.store"                                            11d20s19
      read(1)nsymb,idorel,ngaus,natom,nwcont,numeminus,lmax,nbasallp,   11d20s19
     $     multh,ascale,potdws,ipropsym,istinfo,nextradata              5d25s21
      nsymbb=nsymb                                                      11d20s19
      write(6,*)('nbasallp from readin: '),nbasallp
      write(6,*)('nsymb = '),nsymb
      write(6,*)('idorel = '),idorel,(' ascale '),ascale
      write(6,*)('nuclear repulsion potential: '),potdws
      write(6,*)('numeminus = '),numeminus
      write(6,*)('istinfo: '),istinfo
      write(6,*)('ngaus = '),ngaus
      write(6,*)('natom = '),natom
      write(6,*)('lmax = '),lmax
      if(iflag.eq.0)then                                                6d3s21
       read(1)(nbasdws(isb),isb=1,nsymb)                                 11d20s19
      else                                                              6d3s21
       read(1)(nbasdws(isb),isb=1,nsymb),                               6d3s21
     $      (nlorcont(isb),isb=1,nsymb)                                 6d3s21
      end if                                                            6d3s21
      read(1)(nbasc(isb),isb=1,nsymb)                                   11d20s19
      read(1)(nbasb(isb),isb=1,nsymb)                                   11d20s19
c
c     non rel:
c     nbasdws=nbasc, and if uncontracted, nbasc=nbasb?
c     rel:
c     nbasdws=number of electronic states
c     nbasc number of electronic and positronic states
c     nbasb number of unique gaussians.
c     if uncontracted, nbasb=nbasdws=nbasc/2
c
      write(6,*)('nbasdws: '),(nbasdws(isb),isb=1,nsymb)
      write(6,*)('nbasc: '),(nbasc(isb),isb=1,nsymb)
      write(6,*)('nbasb: '),(nbasb(isb),isb=1,nsymb)
      if(iflag.ne.0)then                                                2d24s21
       write(6,*)('cas input')                                          2d24s21
       write(6,*)('set nbasdws to nbasb: ')
       do isb=1,nsymb                                                   2d24s21
        nbasdws(isb)=nbasb(isb)                                         2d24s21
       end do                                                           2d24s21
      end if                                                            2d24s21
      nextradata=nextradata+1                                           1d10s19
      iextradata=ibcoff                                                 1d10s19
      ibcoff=iextradata+nextradata                                      1d10s19
      nexpect=1+2*natom*3+nextradata+nsymb*2+1
      read(1)isym,((xcart(i,j),atnum(i,j),i=1,3),j=1,natom),            1d19s23
     $     (bc(iextradata+i),i=0,nextradata-1),                         1d10s19
     $     (idoub(isb),iacto(isb),isb=1,nsymb),nlzzq                    1d2s20
      do isb=1,nsymb                                                    5d4s23
       ipack4=idoub(isb)                                                5d4s23
       idoub(isb)=ipack2(1)                                             5d4s23
       icanog(isb)=ipack2(2)                                            5d4s23
      end do                                                            5d4s23
      write(6,*)('atoms: ')                                             1d22s10
      if(idorel.eq.0)then                                               8d17s15
       write(6,*)(' #   Z             x              y              z'),
     $     (' (a.u.)')
      else
       write(6,*)('Nuclei will have finite radii ')
       write(6,310)
  310  format(2x,'#',3x,'Z',13x,'x',14x,'y',14x,'z',7x,                 9d12s17
     $      'nuc weigt (amu) nuc rad   nuc gaus (a.u.)')                9d12s17
      end if                                                            8d17s15
      do i=1,natom                                                      2d24s21
       nuc=nint(atnum(1,i))                                             2d24s21
       atom(i)=element(nuc)                                             4d30s21
       if(idorel.eq.0)then                                              2d24s21
        write(6,10)i,atnum(1,i),element(nuc),                            2d24s21
     $       (xcart(j,i),j=1,3)
       else                                                             2d24s21
        write(6,10)i,atnum(1,i),element(nuc),                            2d24s21
     $       (xcart(j,i),j=1,3),atnum(2,i)                              2d24s21
       end if                                                           2d24s21
   10  format('>',i3,f5.0,x,a3,2x,3f15.8,23x,es11.3)                    7d18s24
      end do                                                            2d24s21
      write(6,*)(' ')                                                   2d24s21
      read(1)((iapair(j,i),j=1,3),i=1,natom)                            11d20s19
      nbdat=ngaus*9+nwcont                                              11d20s19
      ibdat=ibcoff                                                      11d20s19
      ibcoff=ibdat+nbdat                                                11d20s19
      read(1)(bc(ibdat+i),i=0,nbdat-1)                                  11d20s19
      ibstor=ibcoff                                                     11d20s19
      isstor=ibstor+nbasallp*3                                          11d21s19
      ibcoff=isstor+nbasallp*3                                          11d21s19
      nbasis=nbasallp                                                   11d20s19
      read(1)(ibc(ibstor+i),ibc(isstor+i),i=0,nbasallp*3-1)             11d21s19
      if(iflagr.eq.1)return                                             8d12s22
      do isb=1,nsymb                                                    7d26s22
       if(nbasdws(isb).gt.0)then                                        7d26s22
        read(1)dum                                                      7d26s22
        read(1)dum                                                      7d26s22
       end if                                                           7d26s22
      end do                                                            7d26s22
      read(1)nec,jff0,mdon,mdoop,nct                                    7d26s22
      return                                                            11d20s19
      end                                                               11d20s19
