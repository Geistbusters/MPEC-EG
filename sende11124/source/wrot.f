c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine wrot(iwave,mdon,mdoop,nvirt,multh,nsymb,vecon,vecno,bc,11d9s22
     $     ibc)                                                         11d9s22
      implicit real*8 (a-h,o-z)
      integer*4 ipack4(2)
      integer*8 ipack8                                                  8d18s21
      equivalence (ipack8,ipack4)                                       8d18s21
      dimension iwave(*),nvirt(*),multh(8,8),vecon(*),vecno(*)          2d13s22
      include "common.store"                                            8d18s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      isymmrci=iwave(2)
      nroot=iwave(3)
      nff0=iwave(4)+iwave(9)                                            8d18s21
      jff0=iwave(4)+iwave(10)                                           8d18s21
      nec=iwave(7)                                                      7d27s21
      iveci=iwave(4)+iwave(13)                                          7d27s21
      ipack8=ibc(iveci)                                                 5d12s21
      ncsft=ipack4(1)                                                   8d18s21
      iveci=iveci+1                                                     8d18s21
      imff1=iveci+nroot*(ncsft+1)                                       8d18s21
      itmp=ibcoff                                                       2d13s22
      ibcoff=itmp+nroot*ncsft                                           2d13s22
      call enough('wrot.  1',bc,ibc)
      call dgemm('n','n',ncsft,nroot,nroot,1d0,bc(iveci),ncsft,         2d13s22
     $     vecon,nroot,0d0,bc(itmp),ncsft,                              2d13s22
     d' wrot.  1')
      do i=0,ncsft*nroot-1                                              8d21s21
       bc(iveci+i)=bc(itmp+i)                                           2d13s22
      end do                                                            8d21s21
      ibcoff=itmp                                                       2d13s22
      iovr=ibcoff                                                       8d18s21
      ibcoff=iovr+nroot*nroot
      call enough('wrot.  2',bc,ibc)
      do i=iovr,ibcoff-1                                                8d21s21
       bc(i)=0d0                                                        8d21s21
      end do                                                            8d21s21
      call ilimts(1,ncsft,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)     8d21s21
      il=il-1                                                           8d21s21
      ih=ih-1                                                           8d21s21
      do ik=0,nroot-1                                                   8d21s21
       iadk=iveci+ncsft*ik                                              8d21s21
       do ib=0,nroot-1                                                  8d21s21
        iadb=iveci+ncsft*ib                                             8d21s21
        jovr=iovr+ib+nroot*ik                                           8d21s21
        do i=il,ih                                                      8d21s21
         bc(jovr)=bc(jovr)+bc(iadb+i)*bc(iadk+i)                        8d21s21
        end do                                                          8d21s21
       end do                                                           8d21s21
      end do                                                            8d21s21
      irsum=ibcoff                                                       8d21s21
      ibcoff=irsum+nroot*nroot                                          8d21s21
      call enough('wrot.  3',bc,ibc)
      do i=0,nroot*nroot-1                                              8d21s21
       bc(irsum+i)=bc(iovr+i)                                           8d21s21
      end do                                                            8d21s21
      call dws_gsumf(bc(irsum),nroot*nroot)                             8d21s21
      mff1=ibc(imff1)                                                   7d9s21
      inext=imff1+1                                                     7d21s21
      if(mff1.gt.0)then                                                 8d18s21
       ihsdiag=imff1+1                                                  7d9s21
       nff1=ihsdiag+nsymb*mdoop*2                                       7d9s21
       iff1=nff1+nsymb*mdoop                                            7d9s21
       icsf=iff1+ibc(imff1)                                             7d9s21
       inext=icsf+mdoop-mdon                                            7d21s21
       call wrot1(ibc(ihsdiag),ibc(nff1),ibc(icsf),mdon,mdoop,          2d13s22
     $     isymmrci,nvirt,nsymb,multh,bc(iovr),nroot,vecno,bc,ibc)      11d9s22
       do i=0,nroot*nroot-1                                              8d21s21
        bc(irsum+i)=bc(iovr+i)                                           8d21s21
       end do                                                            8d21s21
       call dws_gsumf(bc(irsum),nroot*nroot)                             8d21s21
      end if                                                            8d18s21
      mff2=ibc(inext)                                                   8d19s21
      if(mff2.gt.0)then                                                 8d19s21
       ihddiag=inext+1                                                  7d21s21
       nff2=ihddiag+nsymb*mdoop*2                                       7d21s21
       iff2=nff2+nsymb*mdoop                                            7d21s21
       icsf=iff2+mff2                                                   7d21s21
       icsf2=icsf+mdoop-mdon                                            7d21s21
       call wrot2(ibc(ihddiag),ibc(nff2),ibc(icsf),ibc(icsf2),mdon,     2d13s22
     $     mdoop,isymmrci,nvirt,nsymb,multh,bc(iovr),nroot,vecno,bc,ibc)11d9s22
      else if(mff2.lt.0)then                                            8d19s21
       mdoubstore=inext+1
       ipack8=ibc(mdoubstore)                                           8d3s21
       ndoub=ipack4(1)                                                  8d12s21
       mdoub=ipack4(2)                                                  8d3s21
       mff2a=iabs(mff2)                                                 8d3s21
       iff2=mdoubstore+1                                                8d3s21
       nff2=iff2+mff2a                                                  8d3s21
       nfdat=nff2+mdoop*nsymb                                           8d3s21
       ivdk=nfdat+10*nsymb                                               8d3s21
       ivdknon=ivdk+nroot*ndoub                                         8d19s21
       call wrot2c(bc(ivdk),bc(ivdknon),ibc(nfdat),isymmrci,nvirt,nsymb,
     $      multh,bc(iovr),nroot,vecon,bc,ibc)                          11d9s22
      end if                                                            8d19s21
      do i=0,nroot*nroot-1                                              8d21s21
       bc(irsum+i)=bc(iovr+i)                                           8d21s21
      end do                                                            8d21s21
      call dws_gsumf(bc(irsum),nroot*nroot)                             8d21s21
      ibcoff=iovr
      return                                                            8d18s21
      end                                                               8d18s21
