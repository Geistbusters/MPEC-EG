c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine sortoutsym(nlzz,iorbsym,orbh0,vecob,                   4d23s21
     $     nbasisc,nbasdws,iflag,ixsb,nxsb,mxsb,natom,ngaus,            12d4s22
     $     ibdat,isym,iapair,ibstor,isstor,idorel,ascale,multh,nbasisp, 12d4s22
     $     isb,nsymb,bc,ibc)                                            12d4s22
c
c     iflag=0: vecob is in orthogonal basis                             5d5s21
c     iflag=1: vecob is in ao basis                                     5d5s21
c
      implicit real*8 (a-h,o-z)                                         4d22s21
      logical ldebug
      integer*8 ixsb(mxsb,*)                                            6d7s21
      dimension orbh0(nbasisc,*),vecob(*),iptr(2),nbasisp(*),           1d19s23
     $     iorb(8),noc(8),ixmtr(8,2),islz(3),iorbq(8),iorbqz(8)         12d4s22
      include "common.store"                                            4d22s21
      include "common.print"                                            4d25s21
      if(iprtr(21).ne.0)then                                            4d25s21
       ldebug=.true.
      else                                                              4d25s21
       ldebug=.false.
      end if                                                            4d25s21
      if(ldebug)write(6,*)('sortoutsym: iflag = '),iflag,mxsb,loc(mxsb)
      ipt=ibcoff                                                        4d22s21
      inuq=ipt+nbasisc                                                  4d22s21
      ibcoff=inuq+nbasisc                                               4d22s21
      call enough('sortoutsym.  1',bc,ibc)
      jnuq=inuq-1                                                       4d22s21
      nuq=0                                                             4d22s21
      if(ldebug)write(6,88)(ibc(iorbsym+i),i=0,nbasisc-1)               4d23s21
      do i=0,nbasisc-1                                                  4d22s21
       iqu=ibc(iorbsym+i)                                               4d23s21
       if(ldebug)write(6,*)i+1,ibc(iorbsym+i)                           4d23s21
       do j=1,nuq                                                       4d22s21
        if(iqu.eq.ibc(jnuq+j))then                                      4d22s21
         ibc(ipt+i)=j-1                                                 4d22s21
         go to 1                                                        4d22s21
        end if                                                          4d22s21
       end do                                                           4d22s21
       ibc(ipt+i)=nuq                                                   4d22s21
       nuq=nuq+1                                                        4d22s21
       ibc(jnuq+nuq)=iqu                                                4d22s21
    1  continue                                                         4d22s21
      end do
      if(ldebug)then
       write(6,*)('unique orbital symmetries ')
       do i=1,nuq
        write(6,*)i,ibc(jnuq+i)
       end do
      end if
      if(ldebug)then
       write(6,*)('what we have for orbh0 ')
       call prntm2(orbh0,nbasisc,nbasdws,nbasisc)
      end if
      isz=inuq+nuq                                                      4d22s21
      iqn=isz+nuq                                                       4d22s21
      inbr=iqn+nbasdws                                                  4d22s21
      ibcoff=inbr+nuq                                                   4d22s21
      call enough('sortoutsym.  2',bc,ibc)
      jqn=iqn-1                                                         4d22s21
      jpt=ipt-1                                                         4d22s21
      jsz=isz-1                                                         4d22s21
      jnbr=inbr-1                                                       4d22s21
      do i=1,nuq                                                        4d22s21
       ibc(jnbr+i)=0                                                    4d22s21
      end do                                                            4d22s21
      do i=1,nbasdws                                                    4d22s21
       do j=0,nuq-1                                                     4d22s21
        bc(isz+j)=0d0                                                   4d22s21
       end do                                                           4d22s21
       do j=1,nbasisc                                                   4d22s21
        bc(isz+ibc(jpt+j))=bc(isz+ibc(jpt+j))+orbh0(j,i)**2             4d22s21
       end do                                                           4d22s21
    2  format(2i5,10es10.3)
       xmx=0d0                                                          4d22s21
       do j=1,nuq                                                       4d22s21
        if(bc(jsz+j).gt.xmx)then                                        4d22s21
         xmx=bc(jsz+j)
         ibc(jqn+i)=j
        end if
       end do
       if(ldebug)write(6,2)i,ibc(jqn+i),(bc(isz+j),j=0,nuq-1)
       ibc(jnbr+ibc(jqn+i))=ibc(jnbr+ibc(jqn+i))+1
      end do                                                            4d22s21
      if(ldebug)write(6,*)(ibc(jnbr+i),i=1,nuq)
      ivr=ibcoff                                                        4d22s21
      jvr=ivr+nuq                                                       4d22s21
      ibcoff=jvr+nuq                                                    4d22s21
      call enough('sortoutsym.  3',bc,ibc)
      do i=0,nuq-1                                                      4d22s21
       ibc(ivr+i)=ibcoff                                                4d22s21
       ibc(jvr+i)=ibcoff                                                4d22s21
       ibcoff=ibcoff+nbasdws*ibc(inbr+i)                                4d22s21
      end do                                                            4d22s21
      call enough('sortoutsym.  4',bc,ibc)
      if(iflag.ne.0)then                                                5d5s21
       if(ldebug)then                                                   5d5s21
        write(6,*)('vecao')                                             5d5s21
        call prntm2(vecob,nbasisc,nbasdws,nbasisc)                      5d5s21
       end if                                                           5d5s21
       iptr(2)=ibcoff                                                   5d5s21
       ibcoff=iptr(2)+nbasdws                                           5d5s21
       call enough('sortoutsym.  5',bc,ibc)
       do i=1,nbasdws                                                    4d22s21
        do j=0,nuq-1                                                     4d22s21
         bc(isz+j)=0d0                                                   4d22s21
        end do                                                           4d22s21
        do j=1,nbasisc                                                   4d22s21
         iadr=j+nbasisc*(i-1)                                           5d5s21
         bc(isz+ibc(jpt+j))=bc(isz+ibc(jpt+j))+vecob(iadr)**2           1d19s23
        end do                                                           4d22s21
        xmx=0d0                                                          4d22s21
        im=i-1                                                          5d5s21
        do j=1,nuq                                                      5d5s21
         if(bc(jsz+j).gt.xmx)then                                       5d5s21
          xmx=bc(jsz+j)                                                  5d5s21
          ibc(iptr(2)+im)=j                                              5d5s21
         end if
        end do
        if(ldebug)write(6,2)i,ibc(iptr(2)+im),(bc(isz+j),j=0,nuq-1)      5d5s21
       end do                                                            4d22s21
      else
       if(ldebug)then
        write(6,*)('vecob ')
        call prntm2(vecob,nbasdws,nbasdws,nbasdws)
       end if
       itestr=ibcoff                                                     4d23s21
       ibcoff=itestr+nbasdws*nbasdws                                     4d23s21
       call enough('sortoutsym.  6',bc,ibc)
       xlarg=0d0                                                        6d6s21
       do i=1,nbasdws                                                    4d23s21
        iad=itestr-1+nbasdws*(i-1)                                       4d23s21
        do j=1,nbasdws                                                   4d23s21
         ji=j+nbasdws*(i-1)                                             1d19s23
         if(abs(vecob(ji)).gt.1d-12)then                                1d19s23
          bc(iad+j)=1d0                                                  4d23s21
         else                                                            4d23s21
          xlarg=max(xlarg,abs(vecob(ji)))                               1d19s23
          vecob(ji)=0d0                                                 1d19s23
          bc(iad+j)=0d0                                                  4d23s21
         end if                                                          4d23s21
        end do                                                           4d23s21
       end do                                                            4d23s21
       if(ldebug)write(6,*)('largest component neglected = '),xlarg     6d7s21
       iptr(1)=ibcoff                                                    4d23s21
       iptr(2)=iptr(1)+nbasdws                                           4d23s21
       itemp=iptr(2)+nbasdws                                             4d23s21
       itemp2=itemp+nbasdws*nbasdws                                      4d23s21
       ibcoff=itemp2+nbasdws*nbasdws                                     4d23s21
       call enough('sortoutsym.  7',bc,ibc)
       do i=0,nbasdws-1                                                  4d23s21
        do j=0,nbasdws-1                                                 4d23s21
         ji=itestr+j+nbasdws*i                                           4d23s21
         ij=itemp+i+nbasdws*j                                            4d23s21
         bc(ij)=bc(ji)                                                   4d23s21
        end do                                                           4d23s21
       end do                                                            4d23s21
       do ipass=1,2                                                      4d23s21
        if(ipass.eq.1)then                                               4d23s21
         call dgemm('n','n',nbasdws,nbasdws,nbasdws,1d0,bc(itemp),       4d23s21
     $      nbasdws,bc(itestr),nbasdws,0d0,bc(itemp2),nbasdws,          4d23s21
     d' sortoutsym.  1')
        else                                                             4d23s21
         call dgemm('n','n',nbasdws,nbasdws,nbasdws,1d0,bc(itestr),      4d23s21
     $      nbasdws,bc(itemp),nbasdws,0d0,bc(itemp2),nbasdws,           4d23s21
     d' sortoutsym.  2')
        end if                                                           4d23s21
        ihit=ibcoff
        ibcoff=ihit+nbasdws
        ihit2=ibcoff                                                      8d18s14
        ibcoff=ihit2+nbasdws                                                 8d18s14
        call enough('sortoutsym.  8',bc,ibc)
        jhit=ihit-1
        jhit2=ihit2-1
        do i=1,nbasdws                                                    4d23s21
         ibc(jhit+i)=0
        end do                                                            4d23s21
        thrs=1d-12                                                        4d23s21
        ngrp=1                                                            8d18s14
        is=1                                                              8d18s14
        ibc(ihit)=1                                                       8d18s14
        ibc(ihit2)=1                                                      8d18s14
        nsofar=1                                                          8d18s14
 1908   continue                                                          8d18s14
        nadd=0
        do j=1,nsofar                                                     8d18s14
         do i=is+1,nbasdws                                                    8d18s14
          if(ibc(jhit+i).eq.0)then                                         8d18s14
           iad=itemp2+i-1+nbasdws*(ibc(jhit2+j)-1)                        4d23s21
           if(abs(bc(iad)).gt.thrs.and.i.ne.ibc(jhit2+j))                 4d23s21
     $       then                                                       4d23s21
            ibc(jhit+i)=ngrp                                              8d18s14
            ibc(ihit2+nsofar)=i                                           8d18s14
            nsofar=nsofar+1                                               8d18s14
            nadd=nadd+1                                                   8d18s14
           end if                                                         8d18s14
          end if                                                           8d18s14
         end do                                                           8d18s14
        end do                                                            8d18s14
        if(nadd.ne.0)go to 1908                                           8d18s14
        do i=is,nbasdws                                                      8d18s14
         ibc(iptr(ipass)+i-1)=ibc(jhit+i)                                4d23s21
         if(ibc(jhit+i).eq.0)then
          ngrp=ngrp+1                                                     8d18s14
          is=i                                                            8d18s14
          ibc(jhit+is)=ngrp                                               8d18s14
          ibc(ihit2)=is                                                   8d18s14
          nsofar=1                                                        8d18s14
          go to 1908                                                      8d18s14
         end if                                                           8d18s14
        end do                                                            8d18s14
        ibcoff=ihit                                                      4d23s21
       end do                                                            4d23s21
       ibcoff=itestr                                                     4d23s21
      if(ngrp.gt.100)then
       write(6,*)('ngrp has to be bad!!!')
       call dws_synca
       call dws_finalize
       stop
      end if
      if(ngrp.ne.nuq)then                                               5d27s21
       igrpnn=ibcoff                                                    5d27s21
       ibcoff=igrpnn+ngrp                                               5d27s21
       call enough('sortoutsym.  9',bc,ibc)
       do i=igrpnn,ibcoff-1                                             5d27s21
        ibc(i)=0                                                        5d27s21
       end do                                                           5d27s21
       do i=0,nbasdws-1                                                 5d27s21
        iad=igrpnn+ibc(iptr(2)+i)-1                                     5d27s21
        ibc(iad)=ibc(iad)+1                                             5d27s21
       end do                                                           5d27s21
       if(nuq.eq.2.and.ngrp.eq.4)then                                   10d13s22
        itry1=ibc(igrpnn)+ibc(igrpnn+2)                                 5d27s21
        itry2=ibc(igrpnn+1)+ibc(igrpnn+3)                               5d27s21
        itry3=ibc(igrpnn)+ibc(igrpnn+3)                                 5d27s21
        itry4=ibc(igrpnn+1)+ibc(igrpnn+2)                               5d27s21
        itry5=ibc(igrpnn)+ibc(igrpnn+1)                                 9d16s22
        itry6=ibc(igrpnn+2)+ibc(igrpnn+3)                               9d16s22
        if(itry1.ne.itry3)then                                          5d27s21
         if(itry1.eq.ibc(jnbr+1).and.itry2.eq.ibc(jnbr+2))then          5d27s21
          do i=0,nbasdws-1                                              5d27s21
           if(ibc(iptr(2)+i).eq.3)ibc(iptr(2)+i)=1                      5d27s21
           if(ibc(iptr(2)+i).eq.4)ibc(iptr(2)+i)=2                      5d27s21
          end do                                                        5d27s21
         else if(itry3.eq.ibc(jnbr+1).and.itry4.eq.ibc(jnbr+2))then     5d27s21
          do i=0,nbasdws-1                                              5d27s21
           if(ibc(iptr(2)+i).eq.3)ibc(iptr(2)+i)=2                      5d27s21
           if(ibc(iptr(2)+i).eq.4)ibc(iptr(2)+i)=1                      5d27s21
          end do                                                        5d27s21
         else if(itry5.eq.ibc(jnbr+1).and.itry6.eq.ibc(jnbr+2))then     9d16s22
          do i=0,nbasdws-1                                              9d16s22
           if(ibc(iptr(2)+i).eq.2)then                                  9d16s22
            ibc(iptr(2)+i)=1                                            9d16s22
           else if(ibc(iptr(2)+i).eq.3)then                             9d16s22
            ibc(iptr(2)+i)=2                                            9d16s22
           else if(ibc(iptr(2)+i).eq.4)then                             9d16s22
            ibc(iptr(2)+i)=2                                            9d16s22
           end if                                                       9d16s22
          end do                                                        9d16s22
         else                                                           5d27s21
          write(6,*)('I don''t know how to map')                        5d27s21
         end if                                                         5d27s21
        end if                                                          5d27s21
       else if(nuq.eq.2.and.ngrp.eq.3)then                              5d27s21
        itry1=ibc(igrpnn)                                               5d27s21
        itry2=ibc(igrpnn+1)+ibc(igrpnn+2)                               5d27s21
        itry3=ibc(igrpnn)+ibc(igrpnn+1)                                 5d27s21
        itry4=ibc(igrpnn+2)                                             5d27s21
        itry5=ibc(igrpnn)+ibc(igrpnn+2)                                 2d24s22
        itry6=ibc(igrpnn+1)                                             2d24s22
        if(itry1.ne.itry3)then                                          5d27s21
         if(itry1.eq.ibc(jnbr+1).and.itry2.eq.ibc(jnbr+2))then          5d27s21
          do i=0,nbasdws-1                                              5d27s21
           if(ibc(iptr(2)+i).eq.3)ibc(iptr(2)+i)=2                      5d27s21
          end do                                                        5d27s21
         else if(itry3.eq.ibc(jnbr+1).and.itry4.eq.ibc(jnbr+2))then     5d27s21
          do i=0,nbasdws-1                                              5d27s21
           if(ibc(iptr(2)+i).eq.2)ibc(iptr(2)+i)=1                      5d27s21
           if(ibc(iptr(2)+i).eq.3)ibc(iptr(2)+i)=2                      5d27s21
          end do                                                        5d27s21
         else if(itry5.eq.ibc(jnbr+1).and.itry6.eq.ibc(jnbr+2))then     5d27s21
          do i=0,nbasdws-1                                              2d24s22
           if(ibc(iptr(2)+i).eq.3)ibc(iptr(2)+i)=1                      2d24s22
          end do                                                        5d27s21
         else                                                           5d27s21
          write(6,*)('I don''t know how to map')                        5d27s21
         end if                                                         5d27s21
        else                                                            5d27s21
         write(6,*)('I don''t know how to map')                         5d27s21
        end if                                                          5d27s21
       else if(nuq.eq.1)then                                            10d11s22
        do i=0,nbasdws-1                                                10d11s22
         if(ibc(iptr(2)+i).ne.1)ibc(iptr(2)+i)=1                        10d11s22
        end do                                                          10d11s22
       else
        write(6,*)('I don''t know how to map')                          5d27s21
       end if                                                           5d27s21
       ibcoff=igrpnn                                                    5d27s21
      end if                                                            5d27s21
      end if                                                            5d5s21
      if(ldebug)then                                                    4d23s21
       write(6,*)('qn for rows ')                                        4d23s21
       write(6,88)(ibc(iptr(2)+i),i=0,nbasdws-1)                         4d23s21
       write(6,*)('vs. ')
       write(6,88)(ibc(jqn+i),i=1,nbasdws)                               4d23s21
      end if                                                            4d23s21
      imap=ibcoff                                                       4d23s21
      ibcoff=imap+nbasdws                                               4d23s21
      call enough('sortoutsym. 10',bc,ibc)
      do j=0,nbasdws-1                                                  4d23s21
       if(ibc(iptr(2)+j).gt.nuq)ibc(iptr(2)+j)=nuq                      4d23s21
      end do                                                            4d23s21
      do i=1,nbasdws                                                    4d23s21
       im=i-1                                                           4d23s21
       do j=0,nbasdws-1                                                 4d23s21
        if(ibc(jqn+i).eq.ibc(iptr(2)+j))then                            4d23s21
         ibc(imap+im)=j+1                                               4d23s21
         ibc(iptr(2)+j)=0                                               4d23s21
         go to 89                                                       4d23s21
        end if                                                          4d23s21
       end do                                                           4d23s21
       write(6,*)('unable to map '),i,ibc(jqn+i)
       write(6,*)('qn for rows ')                                        4d23s21
       write(6,88)(ibc(iptr(2)+j),j=0,nbasdws-1)                         4d23s21
       write(6,*)('vs. ')
       write(6,88)(ibc(jqn+j),j=1,nbasdws)                               4d23s21
       stop 'sortoutsym'
   89  continue                                                         4d23s21
      end do                                                            4d23s21
      if(ldebug)then                                                    4d23s21
       write(6,*)('mapping ')
       do i=1,nbasdws                                                    4d23s21
        im=i-1                                                           4d23s21
        write(6,*)i,ibc(imap+im)                                         4d23s21
       end do                                                            4d23s21
      end if                                                            4d23s21
      itmp=ibcoff                                                       4d23s21
      ibcoff=itmp+nbasdws                                               4d23s21
      call enough('sortoutsym. 11',bc,ibc)
      jtmp=itmp-1                                                       4d23s21
      if(iflag.eq.0)then                                                5d5s21
       do i=1,nbasdws                                                    4d23s21
        do j=0,nbasdws-1                                                 4d23s21
         ji=ibc(imap+j)+nbasdws*(i-1)                                    1d19s23
         bc(itmp+j)=vecob(ji)                                           1d19s23
        end do                                                           4d23s21
        do j=1,nbasdws                                                   4d23s21
         ji=j+nbasdws*(i-1)                                             1d19s23
         vecob(ji)=bc(jtmp+j)                                           1d19s23
        end do                                                           4d23s21
       end do                                                            4d23s21
       iherei=ibcoff                                                    6d7s21
       ibcoff=iherei+2*nxsb                                             6d7s21
       call enough('sortoutsym. 12',bc,ibc)
       iherei0=ibcoff                                                   6d7s21
       do i=1,nxsb                                                      6d7s21
        im=i-1                                                          6d7s21
        ibc(iherei+2*im)=ibcoff                                         6d7s21
        ibc(iherei+2*im+1)=ibcoff                                       6d7s21
        ibcoff=ibcoff+ixsb(2,i)*ixsb(2,i)                               6d7s21
       end do                                                           6d7s21
       call enough('sortoutsym. 13',bc,ibc)
       do i=iherei0,ibcoff-1                                            6d7s21
        bc(i)=0d0                                                       6d7s21
       end do                                                           6d7s21
       iuhere=ibcoff                                                    6d7s21
       ibcoff=iuhere+nuq                                                6d7s21
       juhere=iuhere-1                                                  6d7s21
       iqnc=ibcoff                                                      6d7s21
       ibcoff=iqnc+nbasdws                                              6d7s21
       jqnc=iqnc-1                                                      6d7s21
       call enough('sortoutsym. 14',bc,ibc)
       do icol=1,nbasdws                                                6d7s21
        do j=1,nuq                                                      6d7s21
         bc(juhere+j)=0d0                                               6d7s21
        end do                                                          6d7s21
        do j=1,nbasdws                                                  6d7s21
         lsrow=ibc(jqn+j)                                               6d7s21
         ji=j+nbasdws*(icol-1)
         bc(juhere+lsrow)=bc(juhere+lsrow)+vecob(ji)**2                 1d19s23
        end do                                                          6d7s21
        do j=1,nuq                                                      6d7s21
         if(abs(bc(juhere+j)-1d0).lt.1d-8)then                          6d7s21
          ibc(jqnc+icol)=ibc(jnuq+j)                                    6d7s21
          lscol=ibc(jnuq+j)                                             6d7s21
          go to 3132                                                    6d7s21
         end if                                                         6d7s21
        end do                                                          6d7s21
        write(6,*)('unable to set qns!!!!')                             6d7s21
        stop 'sortoutsym'                                                            6d7s21
 3132   continue                                                        6d7s21
        if(nlzz.eq.6)then                                               10d7s22
         lscolu=lscol/100                                               10d7s22
        else                                                            10d7s22
         lscolu=lscol                                                   10d7s22
        end if                                                          10d7s22
        do j=1,nxsb                                                     6d7s21
         if(lscolu.eq.ixsb(1,j))then                                    10d7s22
          jhere=iherei+2*(j-1)                                          6d7s21
          go to 3131                                                    6d7s21
         end if                                                         6d7s21
        end do                                                          6d7s21
        write(6,*)('lscolu did not match any ixsb')                     10d7s22
        write(6,*)('lscolu = '),lscolu,lscol                            10d7s22
        write(6,*)('vs. available: ')
        do j=1,nxsb
         write(6,*)j,ixsb(1,j)                                          10d7s22
        end do                                                          10d7s22
        stop 'sortoutsym'                                               10d7s22
 3131   continue                                                        6d7s21
        do j=1,nbasdws                                                  6d7s21
         jm=j-1                                                         6d7s21
         ltrial=ibc(jnuq+ibc(jqn+j))                                    6d7s21
         if(ltrial.eq.lscol)then                                        6d7s21
          ji=j+nbasdws*(icol-1)                                         1d19s23
          bc(ibc(jhere))=vecob(ji)                                      1d19s23
          ibc(jhere)=ibc(jhere)+1                                       6d7s21
         else                                                           6d7s21
         end if                                                         6d7s21
        end do                                                          6d7s21
       end do                                                           6d7s21
       do i=1,nxsb                                                      6d7s21
        im=i-1                                                          6d7s21
        if(ibc(iherei+2*im).ne.ibc(iherei+2*im+1))then                  6d7s21
         iad=ibc(iherei+2*im+1)                                         6d7s21
         ibc(iherei+2*im)=iad                                           6d7s21
         if(ixsb(4,i).eq.0)then                                         6d7s21
          do j=0,ixsb(2,i)*ixsb(2,i)-1                                  6d7s21
           bc(ixsb(6,i)+j)=bc(iad+j)                                    6d7s21
          end do                                                        6d7s21
          ixsb(4,i)=1                                                   6d7s21
         else                                                           6d7s21
          do j=0,ixsb(2,i)*ixsb(2,i)-1                                  6d7s21
           bc(iad+j)=bc(ixsb(6,i)+j)                                    6d7s21
          end do                                                        6d7s21
         end if                                                         6d7s21
        end if                                                          6d7s21
       end do                                                           6d7s21
       do icol=1,nbasdws                                                6d7s21
        lscol=ibc(jqnc+icol)                                            6d7s21
        if(nlzz.eq.6)then                                               10d7s22
         lscolu=lscol/100                                               10d7s22
        else                                                            10d7s22
         lscolu=lscol                                                   10d7s22
        end if                                                          10d7s22
        do j=1,nxsb                                                     6d7s21
         if(lscolu.eq.ixsb(1,j))then                                    10d7s22
          jhere=iherei+2*(j-1)                                          6d7s21
          go to 3130                                                    6d7s21
         end if                                                         6d7s21
        end do                                                          6d7s21
        write(6,*)('lscolu did not match any ixsb'),lscolu
        stop 'sortoutsym'                                               10d7s22
 3130   continue                                                        6d7s21
        do j=1,nbasdws                                                  6d7s21
         jm=j-1                                                         6d7s21
         ltrial=ibc(jnuq+ibc(jqn+j))                                    6d7s21
         ji=j+nbasdws*(icol-1)                                          1d19s23
         if(ltrial.eq.lscol)then                                        6d7s21
          vecob(ji)=bc(ibc(jhere))                                      1d19s23
          ibc(jhere)=ibc(jhere)+1                                       6d7s21
         else                                                           6d7s21
          vecob(ji)=0d0                                                 1d19s23
         end if                                                         6d7s21
        end do                                                          6d7s21
       end do                                                           6d7s21
       ibcoff=iherei                                                    6d7s21
      else                                                              5d5s21
       do i=1,nbasisc                                                   5d5s21
        do j=0,nbasdws-1                                                5d5s21
         iad=i+nbasisc*(ibc(imap+j)-1)                                  5d5s21
         bc(itmp+j)=vecob(iad)                                          1d19s23
        end do                                                          5d5s21
        do j=1,nbasdws                                                  5d5s21
         iad=i+nbasisc*(j-1)                                            5d5s21
         vecob(iad)=bc(jtmp+j)                                          1d19s23
        end do                                                          5d5s21
       end do                                                           5d5s21
      end if                                                            5d5s21
   88 format(20(10i4,x))                                                4d23s21
      if(ldebug)then
       if(iflag.eq.0)then                                               5d5s21
        write(6,*)('new vecob ')
        call prntm2(vecob,nbasdws,nbasdws,nbasdws)
       else                                                             5d5s21
        write(6,*)('new vecao ')
        call prntm2(vecob,nbasisc,nbasdws,nbasisc)
       end if                                                           5d5s21
      end if
      ibcoff=ipt                                                        4d22s21
      return
      end
