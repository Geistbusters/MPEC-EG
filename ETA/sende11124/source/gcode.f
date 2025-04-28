c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine gcode(idorbb,nclob,isorbb,nopenb,idorbk,nclok,isorbk,  7d15s19
     $     nopenk,icode,imap,nnot,nx,iflag,nab)                         7d16s19
      implicit real*8 (a-h,o-z)                                         7d15s19
      integer*1 idorbb(*),isorbb(*),idorbk(*),isorbk(*),icode(*),       7d15s19
     $     iorb(64),imap(*),nab(2)                                      7d16s19
      data icall/0/                                                     9d9s19
      save icall
      icall=icall+1
      if(iflag.ne.0)write(6,*)('in gcode: '),nclob,nclok
      nnot=0                                                            7d15s19
      if(iabs(nclob-nclok).gt.2)return                                  7d15s19
      if(iflag.ne.0)write(6,*)('cb: '),(idorbb(j),j=1,nclob)
      if(iflag.ne.0)write(6,*)('ck: '),(idorbk(j),j=1,nclok)
      do i=1,64                                                         7d15s19
       iorb(i)=0                                                        7d15s19
      end do                                                            7d15s19
      nb=1
      nk=1
      idifc=0
      nmin=64                                                           7d15s19
      nmax=0                                                            7d15s19
c
c     code = 1 means o bra
c     code = 2 means o ket
c     code = 3 means o bra and ket
c     code = 4 means c bra
c     code = 5 means c bra o ket
c     code = 6 means c ket                                              6d24s19
c     code = 7 means c ket o bra
c
      if(nclob.eq.0.and.nclok.gt.0)then                                 7d16s19
       do i=1,nclok                                                     7d16s19
        idifc=idifc+1                                                   7d16s19
        ival=idorbk(i)                                                  7d16s19
        nmin=min(nmin,ival)                                             7d16s19
        nmax=max(nmax,ival)                                             7d16s19
        iorb(ival)=6                                                    7d16s19
       end do                                                           7d16s19
      else if(nclob.gt.0.and.nclok.eq.0)then                            7d16s19
       do i=1,nclob                                                     7d16s19
        idifc=idifc+1                                                   7d16s19
        ival=idorbb(i)                                                  7d16s19
        nmin=min(nmin,ival)                                             7d16s19
        nmax=max(nmax,ival)                                             7d16s19
        iorb(ival)=4                                                    7d16s19
       end do                                                           7d16s19
      else if(nclob.gt.0.and.nclok.gt.0)then
    1  continue
       if(iflag.ne.0)
     $      write(6,*)('comparing '),nb,nk,idorbb(nb),idorbk(nk)
        if(idorbb(nb).eq.idorbk(nk))then                                 7d15s19
         iorb(idorbk(nk))=0                                             7d16s19
         nb=nb+1                                                         7d15s19
         nk=nk+1                                                         7d15s19
        else if(idorbb(nb).gt.idorbk(nk))then
         idifc=idifc+1
         if(idifc.gt.4)then                                             8d8s19
          nnot=0                                                         7d15s19
          return                                                         7d15s19
         end if                                                          7d15s19
         iorb(idorbk(nk))=6
         ival=idorbk(nk)                                                 7d15s19
         nmin=min(nmin,ival)                                             7d15s19
         nmax=max(nmax,ival)                                             7d15s19
         iorb(idorbb(nb))=4
         ival=idorbb(nb)                                                 7d15s19
         nmin=min(nmin,ival)                                             7d15s19
         nmax=max(nmax,ival)                                             7d15s19
         nk=nk+1                                                         7d15s19
        else
         idifc=idifc+1                                                     7d15s19
         if(idifc.gt.4)then                                             8d8s19
          nnot=0                                                         7d15s19
          return                                                         7d15s19
         end if                                                          7d15s19
         iorb(idorbb(nb))=4
         ival=idorbb(nb)                                                 7d15s19
         nmin=min(nmin,ival)                                             7d15s19
         nmax=max(nmax,ival)                                             7d15s19
         iorb(idorbk(nk))=6
         ival=idorbk(nk)                                                 7d15s19
         nmin=min(nmin,ival)                                             7d15s19
         nmax=max(nmax,ival)                                             7d15s19
         nb=nb+1                                                         7d15s19
        end if                                                           7d15s19
       if(nb.le.nclob.and.nk.le.nclok)go to 1                           7d15s19
       if(nb.le.nclob)then                                              7d16s19
        do i=nb,nclob                                                   7d16s19
         idifc=idifc+1                                                  7d16s19
         ival=idorbb(i)                                                 7d16s19
         iorb(ival)=4                                                   7d16s19
         nmin=min(nmin,ival)                                            7d16s19
         nmax=max(nmax,ival)                                            7d16s19
        end do                                                          7d16s19
       else if(nk.le.nclok)then                                         7d16s19
        do i=nk,nclok                                                   7d16s19
         idifc=idifc+1                                                  7d16s19
         ival=idorbk(i)                                                 7d16s19
         iorb(ival)=6                                                   7d16s19
         nmin=min(nmin,ival)                                            7d16s19
         nmax=max(nmax,ival)                                            7d16s19
        end do                                                          7d16s19
       end if                                                           7d16s19
      end if                                                            7d16s19
      if(iflag.ne.0)write(6,*)('in gcode, idifc: '),idifc
      if(iflag.ne.0)write(6,*)('orb: '),(iorb(i),i=1,nmax)
      if(iflag.ne.0)write(6,*)('ob: '),(isorbb(j),j=1,nopenb)
      if(iflag.ne.0)write(6,*)('ok: '),(isorbk(j),j=1,nopenk)
      nb=1                                                              7d15s19
      nk=1                                                              7d15s19
      idifo=0                                                           7d15s19
      if(nopenb.eq.0.and.nopenk.gt.0)then                                 7d16s19
       do i=1,nopenk                                                     7d16s19
        idifo=idifo+1                                                   7d16s19
        ival=isorbk(i)                                                  7d16s19
        nmin=min(nmin,ival)                                             7d16s19
        nmax=max(nmax,ival)                                             7d16s19
        if(iorb(ival).eq.4)then                                         7d16s19
         iorb(ival)=5                                                    7d16s19
        else                                                            7d16s19
         iorb(ival)=2                                                   7d16s19
        end if
       end do                                                           7d16s19
      else if(nopenb.gt.0.and.nopenk.eq.0)then                            7d16s19
       do i=1,nopenb                                                     7d16s19
        idifo=idifo+1                                                   7d16s19
        ival=isorbb(i)                                                  7d16s19
        nmin=min(nmin,ival)                                             7d16s19
        nmax=max(nmax,ival)                                             7d16s19
        if(iorb(ival).eq.6)then                                         7d16s19
         iorb(ival)=7                                                   7d16s19
        else                                                            7d16s19
         iorb(ival)=1                                                   7d16s19
        end if                                                          7d16s19
       end do                                                           7d16s19
      else if(nopenb.gt.0.and.nopenk.gt.0)then
    2  continue                                                         7d15s19
        if(iflag.ne.0)write(6,*)('comparing '),nb,nk,isorbb(nb),
     $      isorbk(nk),iorb(isorbb(nb)),iorb(isorbk(nk))
        if(isorbb(nb).eq.isorbk(nk))then                                7d15s19
         ival=isorbb(nb)                                                7d15s19
         nmin=min(nmin,ival)                                             7d15s19
         nmax=max(nmax,ival)                                             7d15s19
         iorb(ival)=3                                                   7d15s19
         nb=nb+1                                                        7d15s19
         nk=nk+1                                                        7d15s19
        else if(isorbb(nb).gt.isorbk(nk))then                           7d15s19
         idifo=idifo+1                                                  7d15s19
         if(idifo.gt.4)then                                             7d15s19
          nnot=0                                                        7d15s19
          return                                                        7d15s19
         end if                                                         7d15s19
         ival=isorbb(nb)                                                7d15s19
         nmin=min(nmin,ival)                                             7d15s19
         nmax=max(nmax,ival)                                             7d15s19
         if(iorb(ival).ge.6)then                                        7d15s19
          iorb(ival)=7                                                  7d15s19
         else                                                           7d15s19
          iorb(ival)=1                                                   7d15s19
         end if                                                         7d15s19
         ival=isorbk(nk)                                                7d15s19
         nmin=min(nmin,ival)                                             7d15s19
         nmax=max(nmax,ival)                                             7d15s19
         nk=nk+1                                                        7d16s19
         if(iorb(ival).ge.4)then                                        7d15s19
          iorb(ival)=5                                                  7d15s19
         else                                                           7d15s19
          iorb(ival)=2                                                   7d15s19
         end if                                                         7d15s19
        else                                                            7d15s19
         idifo=idifo+1                                                  7d15s19
         if(idifo.gt.4)then                                             7d15s19
          nnot=0                                                        7d15s19
          return                                                        7d15s19
         end if                                                         7d15s19
         ival=isorbk(nk)                                                7d15s19
         nmin=min(nmin,ival)                                             7d15s19
         nmax=max(nmax,ival)                                             7d15s19
         if(iorb(ival).ge.4)then                                        7d15s19
          iorb(ival)=5                                                  7d15s19
         else                                                           7d15s19
          iorb(ival)=2                                                   7d15s19
         end if                                                         7d15s19
         ival=isorbb(nb)                                                7d15s19
         nmin=min(nmin,ival)                                             7d15s19
         nmax=max(nmax,ival)                                             7d15s19
         if(iorb(ival).ge.6)then                                        7d15s19
          iorb(ival)=7                                                  7d15s19
         else                                                           7d15s19
          iorb(ival)=1                                                   7d15s19
         end if                                                         7d15s19
         nb=nb+1                                                        7d16s19
        end if                                                          7d15s19
       if(nb.le.nopenb.and.nk.le.nopenk)go to 2                         7d15s19
       if(iflag.ne.0)write(6,*)('we are at bottom of loop with '),
     $      nb,nopenb,nk,nopenk
       if(nb.le.nopenb)then                                             7d16s19
        do i=nb,nopenb                                                  7d16s19
         ival=isorbb(i)                                                 7d16s19
         nmin=min(nmin,ival)                                            7d16s19
         nmax=max(nmax,ival)                                            7d16s19
         if(iorb(ival).ge.6)then                                        7d16s19
          iorb(ival)=7                                                  7d16s19
         else                                                           7d16s19
          iorb(ival)=1                                                  7d16s19
         end if                                                         7d16s19
        end do                                                          7d16s19
       else if(nk.le.nopenk)then                                        7d16s19
        do i=nk,nopenk                                                  7d16s19
         ival=isorbk(i)                                                 7d16s19
         nmin=min(nmin,ival)                                            7d16s19
         nmax=max(nmax,ival)                                            7d16s19
         if(iorb(ival).ge.4)then                                        7d16s19
          iorb(ival)=5                                                  7d16s19
         else                                                           7d16s19
          iorb(ival)=2                                                  7d16s19
         end if                                                         7d16s19
        end do                                                          7d16s19
       end if                                                           7d16s19
      end if                                                            7d16s19
      if(iflag.ne.0)write(6,*)('in gcode, idifo: '),idifo,nmin,nmax
      if(iflag.ne.0)write(6,*)('orb: '),(iorb(i),i=1,nmax)
      ii=0                                                              7d15s19
      nnot=0                                                            7d15s19
      do i=nmin,nmax                                                    7d15s19
       if(iorb(i).ne.0)then                                             7d15s19
        ii=ii+1                                                         7d15s19
        icode(ii)=iorb(i)                                               7d15s19
        imap(ii)=i                                                      7d15s19
        if(iorb(i).ne.3)then                                            7d16s19
         nnot=nnot+1                                                    7d16s19
         nab(min(2,nnot))=i                                             7d16s19
        end if                                                          7d16s19
       end if                                                           7d15s19
      end do                                                            7d15s19
      if(iflag.ne.0)write(6,*)('code: '),(icode(i),i=1,ii)
c
c     if nnot=4 and we get 4 and 6, then this is more than a double,
c     so can it.
c
      if(nnot.eq.4)then                                                 10d2s19
       i4=0                                                             10d2s19
       i6=0                                                             10d2s19
       do i=1,ii                                                        10d2s19
        if(icode(i).eq.4)i4=i4+1                                        10d2s19
        if(icode(i).eq.6)i6=i6+1                                        10d2s19
       end do                                                           10d2s19
       if(i4*i6.gt.0)nnot=50                                            10d2s19
      end if                                                            10d2s19
      if(nnot.eq.0.and.(ii.gt.0.or.nclob.gt.0))then                     4d6s20
       nnot=1                                                           7d16s19
      end if                                                            7d16s19
      if(nnot.gt.4)then                                                 7d15s19
       nnot=0                                                           7d15s19
       return                                                           7d15s19
      end if                                                            7d15s19
      nx=ii                                                             7d15s19
      return                                                            7d15s19
      end                                                               7d15s19
