c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine sospher3(numj,ijjs,ijjm,jjmmm,iwavedat,nspc,elowest,   6d4s21
     $     esave,iqnsave,tm,ndim,tm2,tm3,ndim3,lprint,bc,ibc)           11d9s22
      implicit real*8 (a-h,o-z)                                         5d18s21
      integer*8 ijjs(*),ijjm(*),ipack8                                  5d19s21
      integer*2 ipack2(4),iqnsave(4,*)                                  6d4s21
      character*6 lab                                                   5d19s21
      logical lprint                                                    3d2s22
      equivalence (ipack8,ipack2)                                       5d19s21
      dimension iwavedat(nspc,*),esave(*),tm(ndim,ndim,2),              6d4s21
     $     tm2(ndim,*),tm3(ndim3,*)                                     6d4s21
      include "common.store"
      if(ndim.eq.0)then                                                 3d2s22
       return                                                           3d2s22
      end if                                                            3d2s22
      if(lprint)then                                                    3d8s22
       write(6,*)(' ')                                                  3d8s22
       write(6,*)('in sospher3, coupling electronic states')            3d8s22
       write(6,*)(' ')                                                  3d8s22
       write(6,1215)elowest                                             3d8s22
 1215  format('term',5x,'rt  J',3x,'Etot',15x,'Erel(1/cm) wrt',f14.8)   3d8s22
      end if                                                            3d8s22
      iimmm=jjmmm                                                       5d18s21
      ntot=0                                                            5d19s21
      do j=1,numj                                                       6d3s21
       nn=ijjs(j)                                                       6d3s21
       ntot=ntot+nn                                                     6d3s21
      end do                                                            6d3s21
      ieiga=ibcoff                                                      6d3s21
      ieign=ieiga+ntot                                                  6d3s21
      ieigr=ieign+ntot                                                  6d3s21
      ieigj=ieigr+ntot                                                  6d3s21
      ibcoff=ieigj+ntot                                                 6d3s21
      call enough('sospher3.  1',bc,ibc)
      jeiga=ieiga                                                       6d3s21
      jeign=ieign                                                       6d3s21
      jeigr=ieigr                                                       6d3s21
      jeigj=ieigj                                                       6d3s21
      ircoff=1                                                          6d4s21
      do j=1,numj
       nn=ijjs(j)                                                       5d19s21
       ll=ijjm(j)+nn*nn                                                 5d19s21
       ieig=ibcoff                                                      5d19s21
       ivec=ieig+nn                                                     5d19s21
       isym=ivec+nn*nn                                                  5d19s21
       ibcoff=isym+nn                                                   5d19s21
       call enough('sospher3.  2',bc,ibc)
       call diagx(nn,bc(ijjm(j)),bc(ieig),bc(ivec),ibc(isym),bc,ibc)    11d14s22
       ivect=ibcoff                                                     6d4s21
       itmp=ivect+nn*nn                                                 6d4s21
       ibcoff=itmp+ndim*nn                                               6d4s21
       call enough('sospher3.  3',bc,ibc)
       do ii=0,nn-1                                                     6d4s21
        do jj=0,nn-1                                                    6d4s21
         ji=ivec+jj+nn*ii                                               6d4s21
         ij=ivect+ii+nn*jj                                              6d4s21
         bc(ij)=bc(ji)                                                  6d4s21
        end do                                                          6d4s21
       end do                                                           6d4s21
       do ip=1,2                                                        6d4s21
        call dgemm('n','n',ndim,nn,nn,1d0,tm(1,ircoff,ip),ndim,bc(ivec),
     $      nn,0d0,bc(itmp),ndim,                                       6d4s21
     d' sospher3.  1')
        do ii=0,nn-1                                                     6d4s21
         jtmp=itmp+ndim*ii                                               6d4s21
         do jj=0,ndim-1                                                  6d4s21
          tm(1+jj,ircoff+ii,ip)=bc(jtmp+jj)
         end do                                                         6d4s21
        end do                                                          6d4s21
        call dgemm('n','n',nn,ndim,nn,1d0,bc(ivect),nn,tm(ircoff,1,ip),
     $       ndim,0d0,bc(itmp),nn,                                      6d4s21
     d' sospher3.  2')
        do ii=0,ndim-1                                                   6d4s21
         jtmp=itmp+nn*ii                                                 6d4s21
         do jj=0,nn-1                                                    6d4s21
          tm(ircoff+jj,ii+1,ip)=bc(jtmp+jj)                                6d4s21
         end do                                                         6d4s21
        end do                                                          6d4s21
       end do                                                           6d4s21
       ibcoff=itmp+ndim3*nn                                             6d4s21
       call enough('sospher3.  4',bc,ibc)
       if(ndim3.gt.0)then                                               3d2s22
        call dgemm('n','n',ndim3,nn,nn,1d0,tm3(1,ircoff),ndim3,bc(ivec), 6d4s21
     $      nn,0d0,bc(itmp),ndim3,                                      6d4s21
     d' sospher3.  3')
        do ii=0,nn-1                                                     6d4s21
         jtmp=itmp+ndim3*ii                                              6d4s21
         do jj=0,ndim3-1                                                 6d4s21
          tm3(1+jj,ircoff+ii)=bc(jtmp+jj)                                6d4s21
         end do                                                          6d4s21
        end do                                                           6d4s21
        call dgemm('n','n',nn,ndim3,nn,1d0,bc(ivect),nn,tm2(ircoff,1),   6d4s21
     $       ndim,0d0,bc(itmp),nn,                                      6d4s21
     d' sospher3.  4')
        do ii=0,ndim3-1                                                  6d4s21
         jtmp=itmp+nn*ii                                                 6d4s21
         do jj=0,nn-1                                                    6d4s21
          tm2(ircoff+jj,ii+1)=bc(jtmp+jj)                                6d4s21
         end do                                                          6d4s21
        end do                                                           6d4s21
       end if                                                           3d2s22
       ircoff=ircoff+nn                                                 6d4s21
       ibcoff=itmp                                                      6d4s21
       do i=0,nn-1                                                      6d3s21
        bc(jeiga+i)=bc(ieig+i)                                          6d3s21
        ibc(jeigj+i)=iimmm                                              6d3s21
        vx=0d0                                                          6d3s21
        jvec=ivec+nn*i                                                  6d3s21
        do ii=0,nn-1                                                    6d3s21
         if(abs(bc(jvec+ii)).gt.vx)then                                 6d3s21
          vx=abs(bc(jvec+ii))                                           6d3s21
          ivx=ii                                                        6d3s21
         end if                                                         6d3s21
         ibc(jeigr+i)=ibc(ll+ivx)                                        6d3s21
         ibc(jeign+i)=ibc(ll+nn+ivx)                                    6d3s21
        end do                                                          6d3s21
       end do                                                           6d3s21
       jeiga=jeiga+nn                                                   6d3s21
       jeigr=jeigr+nn                                                   6d3s21
       jeign=jeign+nn                                                   6d3s21
       jeigj=jeigj+nn                                                   6d3s21
       if(j.eq.1)then                                                   5d19s21
        e0=bc(ieig)                                                     5d19s21
       else                                                             5d19s21
        e0=min(e0,bc(ieig))                                             5d19s21
       end if                                                           5d19s21
       do i=0,nn-1                                                      5d19s21
        bc(ijjm(j)+i)=bc(ieig+i)                                        5d19s21
       end do                                                           5d19s21
       ibcoff=ieig                                                      5d19s21
       iimmm=iimmm+2                                                    5d18s21
      end do
      e0=elowest                                                        6d1s21
      isort=ibcoff                                                      5d19s21
      isort2=isort+ntot                                                 5d19s21
      isort3=isort2+ntot                                                5d19s21
      ibcoff=isort3+ntot                                                5d19s21
      call enough('sospher3.  5',bc,ibc)
      call dsortdws(bc(ieiga),ibc(isort),ntot)                          1d18s23
      itmp=ibcoff                                                       6d4s21
      ibcoff=itmp+ndim*ndim                                             6d4s21
      do ipass=1,2                                                         6d4s21
       do i=0,ndim-1                                                    6d4s21
        ii=ibc(isort+i)                                                 6d4s21
        jtmp=itmp-1+ndim*i                                              6d4s21
        do j=1,ndim                                                     6d4s21
         bc(jtmp+j)=tm(j,ii,ipass)                                      6d4s21
        end do                                                          6d4s21
       end do                                                           6d4s21
       do i=1,ndim                                                      6d4s21
        jtmp=itmp-1+ndim*(i-1)
        do j=1,ndim
         tm(j,i,ipass)=bc(jtmp+j)                                       6d4s21
        end do                                                          6d4s21
       end do                                                           6d4s21
       do i=1,ndim
        jtmp=itmp-1+ndim*(i-1)
        do j=1,ndim
         jj=ibc(isort+j-1)
         bc(jtmp+j)=tm(jj,i,ipass)
        end do
       end do
       do i=1,ndim                                                      6d4s21
        jtmp=itmp-1+ndim*(i-1)
        do j=1,ndim
         tm(j,i,ipass)=bc(jtmp+j)                                       6d4s21
        end do                                                          6d4s21
       end do                                                           6d4s21
      end do                                                            6d4s21
      ibcoff=itmp+ndim*ndim3                                            6d4s21
      call enough('sospher3.  6',bc,ibc)
      do i=0,ndim-1                                                     6d4s21
       ii=ibc(isort+i)                                                  6d4s21
       jtmp=itmp-1+ndim3*i                                              6d4s21
       do j=1,ndim3                                                     6d4s21
        bc(jtmp+j)=tm3(j,ii)                                             6d4s21
       end do                                                           6d4s21
      end do                                                            6d4s21
      do i=1,ndim                                                       6d4s21
       jtmp=itmp-1+ndim3*(i-1)                                          6d4s21
       do j=1,ndim3                                                     6d4s21
        tm3(j,i)=bc(jtmp+j)                                             6d4s21
       end do                                                           6d4s21
      end do                                                            6d4s21
      do i=1,ndim3                                                      6d4s21
       jtmp=itmp+ndim*(i-1)                                             6d4s21
       do j=0,ndim-1                                                    6d4s21
        jj=ibc(isort+j)                                                 6d4s21
        bc(jtmp+j)=tm2(jj,i)                                            6d4s21
       end do                                                           6d4s21
      end do                                                            6d4s21
      do i=1,ndim3                                                      6d4s21
       jtmp=itmp-1+ndim*(i-1)                                           6d4s21
       do j=1,ndim                                                      6d4s21
        tm2(j,i)=bc(jtmp+j)                                             6d4s21
       end do                                                           6d4s21
      end do                                                            6d4s21
      last=-1                                                           6d3s21
      lastr=-1
      do i=0,ntot-1                                                     6d3s21
       ii=ibc(isort+i)-1                                                6d3s21
       if(last.ne.ibc(ieign+ii).or.lastr.ne.ibc(ieigr+ii))then          6d3s21
        if(last.ne.ibc(ieign+ii))then                                   6d3s21
         do l=1,6                                                          10d27s20
          if(iwavedat(13+l,ibc(ieign+ii)).ne.0)then                           5d19s21
           lab(l:l)=char(iwavedat(13+l,ibc(ieign+ii)))                    6d3s21
          else                                                             10d27s20
           lab(l:l)=' '                                                   5d19s21
          end if                                                           10d27s20
         end do                                                           5d19s21
        end if                                                          6d3s21
       end if
       jj=ibc(ieigj+ii)
       ip=i+1                                                           6d4s21
       esave(ip)=bc(ieiga+i)-e0                                         6d4s21
       iqnsave(1,ip)=ibc(ieign+ii)                                      6d4s21
       iqnsave(2,ip)=ibc(ieigr+ii)                                      6d4s21
       iqnsave(3,ip)=jj                                                 6d4s21
       if(lprint)then                                                   3d2s22
        if(mod(jj,2).eq.0)then
         jj=jj/2
         write(6,22)iwavedat(1,ibc(ieign+ii)),lab,ibc(ieigr+ii),
     $      jj,
     $      bc(ieiga+i),(bc(ieiga+i)-e0)*219474.635d0
        else
         write(6,23)iwavedat(1,ibc(ieign+ii)),lab,ibc(ieigr+ii),
     $      jj,
     $      bc(ieiga+i),(bc(ieiga+i)-e0)*219474.635d0
        end if
       end if                                                           3d2s22
       last=ibc(ieign+ii)                                               6d3s21
       lastr=ibc(ieigr+ii)                                              6d3s21
   22  format(1x,i1,a6,i3,i3,f14.8,f14.2)
   23  format(1x,i1,a6,i3,i3,'/2',f14.8,f14.2)
      end do                                                            6d3s21
      jsort=isort                                                       5d19s21
      jsort3=isort3                                                     5d19s21
      iimmm=jjmmm                                                       5d19s21
      do j=1,numj                                                       5d19s21
       nn=ijjs(j)                                                       5d19s21
       ll=ijjm(j)+nn*nn                                                 5d19s21
       do i=0,nn-1                                                      5d19s21
        bc(jsort+i)=bc(ijjm(j)+i)                                       5d19s21
        ipack2(1)=iimmm                                                 5d19s21
        ipack2(2)=ibc(ll+i)                                             5d19s21
        ipack2(3)=ibc(ll+i+nn)                                          5d19s21
        ibc(jsort3+i)=ipack8                                             5d19s21
       end do                                                           5d19s21
       jsort=jsort+nn                                                   5d19s21
       jsort3=jsort3+nn                                                 5d19s21
       iimmm=iimmm+2                                                    5d19s21
      end do                                                            5d19s21
      call dsortdws(bc(isort),ibc(isort2),ntot)                         1d18s23
      do i=0,ntot-1                                                     5d19s21
       j=ibc(isort2+i)-1                                                5d19s21
       ipack8=ibc(isort3+j)                                             5d19s21
       do l=1,6                                                          10d27s20
        if(iwavedat(13+l,ipack2(3)).ne.0)then                           5d19s21
         lab(l:l)=char(iwavedat(13+l,ipack2(3)))                        5d19s21
        else                                                             10d27s20
         lab(l:l)=' '                                                   5d19s21
        end if                                                           10d27s20
       end do                                                           5d19s21
       itry=ipack2(1)                                                   6d1s21
      end do                                                            5d19s21
      return                                                            5d18s21
      end                                                               5d18s21
