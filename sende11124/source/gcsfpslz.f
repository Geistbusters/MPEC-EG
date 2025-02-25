c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine gcsfpslz(nfcn,hdiag,hdps,ecut,ibasis,ipsbase,ncsf,       7d11s19
     $     ipointf,mdon,nvv,hivs,hiv,ncsft,nlzzu,zzdig,zzpsdig,ipointp, 5d14s21
     $     ndigs,nps,npsf,iptrbit,mysym,nec,mdoo,nsymb,ixlzz,islz,multh,7d19s22
     $     igcode,lprint,norbu,irefou,ismu,irelu,hdiagps,bc,ibc)        11d14s22
      implicit real*8 (a-h,o-z)                                         6d11s19
      logical lprint                                                    7d19s22
      integer*8 ipack8,gandcc,gandco,gandcb                             1d19s23
      integer*4 ipack4(2)                                               8d24s22
      equivalence (ipack8,ipack4)                                       8d24s22
      dimension hdiag(*),hdps(*),ibasis(3,*),ipsbase(3,*),ncsf(*),      7d11s19
     $     ipointf(*),hivs(*),hiv(*),zzdig(nfcn,*),zzpsdig(npsf,*),     5d14s21
     $     ipointp(*),iptrbit(2,mdoo+1,*),nab4(2,3),ixlzz(8,*),islz(*), 7d19s22
     $     multh(8,8),irefou(*),ismu(*),irelu(*),hdiagps(*),ioxx(2)     1d19s23
      data loopx/1000/
      include "common.store"                                            7d19s22
      include "common.mrci"                                             6d10s19
      data icall/0/
      save
      icall=icall+1
      if(nlzzu.eq.6)then                                                2d27s23
       npass=3                                                          2d27s23
      else                                                              2d27s23
       npass=1                                                          2d27s23
      end if                                                            2d27s23
      loop=0
      ipt1=ibcoff
      ipt2=ipt1+nfcn
      ibcoff=ipt2+nfcn
      call enough('gcsfpslz.  1',bc,ibc)
      ncoup=0                                                           8d24s22
      icoup=ibcoff                                                      8d24s22
      jcoup=icoup                                                       8d24s22
      do i=1,nfcn-1                                                     8d24s22
       ncloi=ibasis(1,i)                                                7d19s22
       ncloip=ncloi+1                                                     7d19s22
       iic=iptrbit(1,ncloip,mysym)+ibasis(2,i)-1                        7d19s22
       iio=iptrbit(2,ncloip,mysym)+ibasis(3,i)-1                        7d19s22
       do j=i+1,nfcn                                                    8d24s22
        ncloj=ibasis(1,j)                                                7d19s22
        if(iabs(ncloj-ncloi).le.2)then                                  8d24s22
         nclojp=ncloj+1                                                     7d19s22
         jjc=iptrbit(1,nclojp,mysym)+ibasis(2,j)-1                        7d19s22
         jjo=iptrbit(2,nclojp,mysym)+ibasis(3,j)-1                        7d19s22
         call gandc4(ibc(jjc),ibc(jjo),ibc(iic),ibc(iio),nopenj,nopen,  11d13s20
     $         norbu,nnot4,nab4,bc,ibc)                                 11d14s22
         if(nnot4.ge.3)then
          is12=multh(ismu(nab4(1,1)),ismu(nab4(2,1)))
          is34=multh(ismu(nab4(1,2)),ismu(nab4(2,2)))
          xcoup=0d0                                                     8d24s22
          do ixyz=1,npass                                               2d27s23
           if(is12.eq.islz(ixyz).and.is34.eq.islz(ixyz))then            2d27s23
            iad1=ixlzz(ismu(nab4(1,1)),ixyz)+irelu(nab4(1,1))-1         2d27s23
     $          +irefou(ismu(nab4(1,1)))*(irelu(nab4(2,1))-1)              7d19s22
            iad2=ixlzz(ismu(nab4(1,2)),ixyz)+irelu(nab4(1,2))-1         2d27s23
     $          +irefou(ismu(nab4(1,2)))*(irelu(nab4(2,2))-1)              7d19s22
            xcoup=xcoup+abs(bc(iad1)*bc(iad2))                          2d27s23
           end if                                                       2d27s23
          end do                                                        2d27s23
           if(abs(xcoup).gt.1d-10)then                                  7d19s22
            ipack4(1)=i                                                 8d24s22
            ipack4(2)=j                                                 8d24s22
            ibc(jcoup)=ipack8                                           8d24s22
            jcoup=jcoup+1                                               8d24s22
           end if                                                       8d24s22
          if(nnot4.eq.4.and.abs(xcoup).lt.1d-10)then                    8d24s22
           is14=multh(ismu(nab4(1,1)),ismu(nab4(2,2)))
           is32=multh(ismu(nab4(1,2)),ismu(nab4(2,1)))
           do ixyz=1,npass                                              2d27s23
            if(is14.eq.islz(ixyz).and.is32.eq.islz(ixyz))then           2d27s23
             iad1=ixlzz(ismu(nab4(1,1)),ixyz)+irelu(nab4(1,1))-1        2d27s23
     $          +irefou(ismu(nab4(1,1)))*(irelu(nab4(2,2))-1)              7d19s22
             iad2=ixlzz(ismu(nab4(1,2)),ixyz)+irelu(nab4(1,2))-1        2d27s23
     $          +irefou(ismu(nab4(1,2)))*(irelu(nab4(2,1))-1)              7d19s22
             xcoup=xcoup+abs(bc(iad1)*bc(iad2))                         2d27s23
            end if                                                      2d27s23
           end do                                                       2d27s23
            if(abs(xcoup).gt.1d-10)then                                 8d24s22
             ipack4(1)=i                                                8d24s22
             ipack4(2)=j                                                8d24s22
             ibc(jcoup)=ipack8                                          8d24s22
             jcoup=jcoup+1                                              8d24s22
            end if                                                      8d24s22
          end if                                                        8d24s22
         end if
        end if                                                          8d24s22
       end do                                                           8d24s22
      end do
      ncoup=jcoup-icoup                                                 8d24s22
      ibcoff=jcoup                                                      8d3s22
      ihit=ibcoff                                                       8d3s22
      igroup=ihit+nfcn                                                  8d24s22
      ibcoff=igroup+nfcn                                                8d24s22
      call enough('gcsfpslz.  2',bc,ibc)
      jgroup=igroup-1                                                   8d3s22
      jhit=ihit-1
      do i=0,nfcn-1                                                     8d24s22
       ibc(ihit+i)=0                                                    8d3s22
      end do                                                            8d3s22
   57 continue                                                          8d3s22
       ngroup=0                                                         8d3s22
       do i=0,nfcn-1                                                    8d25s22
        if(ibc(ihit+i).eq.0)then                                        8d3s22
         ibc(igroup)=i+1                                                8d3s22
         ibc(ihit+i)=1                                                  8d3s22
         ngroup=1                                                       8d3s22
         go to 58                                                       8d3s22
        end if                                                          8d3s22
       end do                                                           8d3s22
       go to 59                                                         8d3s22
   58  continue                                                         8d3s22
       ngroupo=ngroup                                                   8d3s22
       do i=0,ncoup-1                                                   8d3s22
        if(ibc(icoup+i).ne.0)then                                       8d3s22
         ipack8=ibc(icoup+i)                                             8d3s22
         j1=-1                                                          8d3s22
         j2=-1                                                          8d3s22
         do j=0,ngroup-1                                                 8d3s22
          if(ipack4(1).eq.ibc(igroup+j))j1=j                            8d3s22
          if(ipack4(2).eq.ibc(igroup+j))j2=j                            8d3s22
         end do                                                         8d3s22
         if(max(j1,j2).gt.-1.and.min(j1,j2).eq.-1)then                  8d3s22
          if(j1.lt.0)then                                               8d3s22
           ibc(igroup+ngroup)=ipack4(1)                                 8d3s22
           ibc(jhit+ipack4(1))=1                                        8d3s22
          else                                                          8d3s22
           ibc(igroup+ngroup)=ipack4(2)                                 8d3s22
           ibc(jhit+ipack4(2))=1                                        8d3s22
          end if                                                        8d3s22
          ngroup=ngroup+1                                               8d3s22
          ibc(icoup+i)=0                                                8d3s22
         end if                                                         8d3s22
        end if                                                          8d3s22
       end do                                                           8d3s22
       if(ngroup.ne.ngroupo)go to 58                                    8d3s22
        eavg=0d0                                                         8d3s22
        do i=0,ngroup-1                                                  8d3s22
         in=ibc(igroup+i)                                                8d3s22
         nclo=ibasis(1,in)
         nclop=nclo+1
         iic=iptrbit(1,nclop,mysym)+ibasis(2,in)-1                        7d19s22
         iio=iptrbit(2,nclop,mysym)+ibasis(3,in)-1                        7d19s22
         eavg=eavg+hdiagps(in)                                          8d24s22
        end do                                                           8d3s22
        eavg=eavg/dfloat(ngroup)                                         8d3s22
        do i=0,ngroup-1                                                  8d3s22
         in=ibc(igroup+i)                                                8d3s22
         hdiagps(in)=eavg                                               8d24s22
        end do                                                           8d3s22
        go to 57                                                        8d24s22
   59  continue                                                         8d24s22
       ibcoff=icoup                                                     8d24s22
      ips=1                                                             6d11s19
      ivv=1                                                             9d4s19
      ivvs=1                                                            9d4s19
      irun=0                                                            7d11s19
      iprun=0                                                           7d11s19
      jpt1=ipt1-1                                                       7d19s22
      jpt2=ipt2-1                                                       7d19s22
      do i=1,nfcn                                                       6d11s19
       edel=hdiagps(i)-ecut                                             8d24s22
       nclo=ibasis(1,i)                                                 7d11s19
       nclop=nclo+1
       iarg=nclop-mdon                                                  7d11s19
       ibc(jpt2+i)=0
       if(edel.le.0d0)then                                              6d11s19
        iprun=iprun+1                                                   7d19s22
        ibc(jpt1+iprun)=i                                               7d19s22
        ibc(jpt2+i)=1                                                   7d19s22
       end if                                                           6d11s19
      end do                                                            6d11s19
      iprun=iprun-1                                                     7d19s22
    1 continue
      iprun=iprun+1                                                     7d19s22
      do i=1,iprun
       ii=ibc(jpt1+i)                                                   7d19s22
       nclo=ibasis(1,ii)                                                7d19s22
       nopen=nec-2*nclo
       nclop=nclo+1                                                     7d19s22
       iic=iptrbit(1,nclop,mysym)+ibasis(2,ii)-1                        7d19s22
       iio=iptrbit(2,nclop,mysym)+ibasis(3,ii)-1                        7d19s22
       do j=1,nfcn
        if(ibc(jpt2+j).eq.0)then                                        7d19s22
         ncloj=ibasis(1,j)                                                7d19s22
         if(iabs(ncloj-nclo).le.2)then
          loop=0
          nclojp=ncloj+1                                                     7d19s22
          nopenj=nec-2*ncloj
          jjc=iptrbit(1,nclojp,mysym)+ibasis(2,j)-1                        7d19s22
          jjo=iptrbit(2,nclojp,mysym)+ibasis(3,j)-1                        7d19s22
          gandcc=ieor(ibc(jjc),ibc(iic))                                1d19s23
          gandco=ieor(ibc(jjo),ibc(iio))                                1d19s23
          gandcb=ior(gandcc,gandco)                                     1d19s23
          ndifb=popcnt(gandcb)                                          1d19s23
          if(ndifb.le.4)then                                            1d19s23
           ndifs=popcnt(gandco)                                         1d19s23
           ndifd=popcnt(gandcc)                                         1d19s23
           nnot4=0                                                      1d19s23
           if(ndifs.eq.4.and.ndifb.eq.4)then                            1d19s23
            nnot4=4                                                     1d19s23
            ioxx(1)=1                                                   1d19s23
            ioxx(2)=1                                                   1d19s23
            do l=1,norbu                                                1d19s23
             if(btest(gandcb,l))then                                    1d19s23
              if((btest(ibc(iic),l).and.btest(ibc(jjo),l)).or.          1d19s23
     $            (btest(ibc(iio),l).and..not.btest(ibc(jjc),l)))then   1d19s23
               nab4(2,ioxx(2))=l                                        1d19s23
               ioxx(2)=ioxx(2)+1                                        1d19s23
              else                                                      1d19s23
               nab4(1,ioxx(1))=l                                        1d19s23
               ioxx(1)=ioxx(1)+1                                        1d19s23
              end if                                                    1d19s23
             end if                                                     1d19s23
            end do                                                      1d19s23
           else if(ndifb.eq.3)then                                      1d19s23
            nnot4=3                                                     1d19s23
            ioxx(1)=1                                                   1d19s23
            ioxx(2)=1                                                   1d19s23
            iswap=0                                                     1d19s23
            do l=1,norbu                                                1d19s23
             if(btest(gandcb,l))then                                    1d19s23
              if(btest(gandcc,l).and.                                   1d19s23
     $        ((btest(ibc(jjc),l).and..not.btest(ibc(iio),l)).or.       1d19s23
     $         (btest(ibc(iic),l).and..not.btest(ibc(jjo),l))))then     1d19s23
               if(btest(ibc(iic),l))iswap=1                             1d19s23
               nab4(1,1)=l                                              1d19s23
               nab4(1,2)=l                                              1d19s23
              else                                                      1d19s23
               nab4(2,ioxx(2))=l                                        1d19s23
               ioxx(2)=ioxx(2)+1                                        1d19s23
              end if                                                    1d19s23
             end if                                                     1d19s23
            end do                                                      1d19s23
            if(iswap.ne.0)then                                          1d19s23
             icpy=nab4(1,1)                                             1d19s23
             nab4(1,1)=nab4(2,1)                                        1d19s23
             nab4(2,1)=icpy                                             1d19s23
             icpy=nab4(1,2)                                             1d19s23
             nab4(1,2)=nab4(2,2)                                        1d19s23
             nab4(2,2)=icpy                                             1d19s23
             nbt=0                                                      1d19s23
             if(btest(ibc(jjc),nab4(1,2)).and.                          1d19s23
     $            .not.btest(ibc(jjc),nab4(1,1)))nbt=1                  1d19s23
            else                                                        1d19s23
             nbt=0                                                      1d19s23
             if(btest(ibc(iic),nab4(2,2)).and.                          1d19s23
     $            .not.btest(ibc(iic),nab4(2,1)))nbt=1                  1d19s23
            end if                                                      1d19s23
            if(nbt.ne.0)then                                            1d19s23
             nab4(1,1)=nab4(1,2)                                        1d19s23
             nab4(2,1)=nab4(2,2)                                        1d19s23
            end if                                                      1d19s23
           else if(ndifs.eq.0.and.ndifd.eq.2)then                       1d19s23
            nnot4=3                                                     1d19s23
            do l=1,norbu                                                1d19s23
             if(btest(gandcb,l))then                                    1d19s23
              if(btest(ibc(jjc),l))then                                 1d19s23
               nab4(1,1)=l                                              1d19s23
               nab4(1,2)=l                                              1d19s23
              else                                                      1d19s23
               nab4(2,1)=l                                              1d19s23
               nab4(2,2)=l                                              1d19s23
              end if                                                    1d19s23
             end if                                                     1d19s23
            end do                                                      1d19s23
           end if                                                       1d19s23
           if(nnot4.ge.3)then
            is12=multh(ismu(nab4(1,1)),ismu(nab4(2,1)))
            is34=multh(ismu(nab4(1,2)),ismu(nab4(2,2)))
            xcoup=0d0                                                   3d21s23
            do ixyz=1,npass                                             2d27s23
             if(is12.eq.islz(ixyz).and.is34.eq.islz(ixyz))then          2d27s23
              iad1=ixlzz(ismu(nab4(1,1)),ixyz)+irelu(nab4(1,1))-1       2d27s23
     $          +irefou(ismu(nab4(1,1)))*(irelu(nab4(2,1))-1)              7d19s22
              iad2=ixlzz(ismu(nab4(1,2)),ixyz)+irelu(nab4(1,2))-1       2d27s23
     $          +irefou(ismu(nab4(1,2)))*(irelu(nab4(2,2))-1)              7d19s22
              xcoup=xcoup+abs(bc(iad1)*bc(iad2))                        2d27s23
             end if                                                     2d27s23
            end do                                                      2d27s23
             if(abs(xcoup).gt.1d-10)then                                 7d19s22
              ibc(ipt1+iprun)=j                                          7d19s22
              ibc(jpt2+j)=1                                              7d19s22
              go to 1                                                    7d19s22
             end if                                                      7d19s22
            if(nnot4.eq.4)then
             is14=multh(ismu(nab4(1,1)),ismu(nab4(2,2)))
             is32=multh(ismu(nab4(1,2)),ismu(nab4(2,1)))
             do ixyz=1,npass                                            2d27s23
              if(is14.eq.islz(ixyz).and.is32.eq.islz(ixyz))then         2d27s23
               iad1=ixlzz(ismu(nab4(1,1)),ixyz)+irelu(nab4(1,1))-1      2d27s23
     $          +irefou(ismu(nab4(1,1)))*(irelu(nab4(2,2))-1)              7d19s22
               iad2=ixlzz(ismu(nab4(1,2)),ixyz)+irelu(nab4(1,2))-1      2d27s23
     $             +irefou(ismu(nab4(1,2)))*(irelu(nab4(2,1))-1)              7d19s22
               xcoup=xcoup+abs(bc(iad1)*bc(iad2))                       2d27s23
              end if                                                    2d27s23
             end do                                                     2d27s23
              if(abs(xcoup).gt.1d-10)then                                 7d19s22
               ibc(ipt1+iprun)=j                                          7d19s22
               ibc(jpt2+j)=1                                              7d19s22
               go to 1                                                    7d19s22
              end if                                                      7d19s22
            end if
           end if                                                       1d19s23
          end if
         end if                                                         7d19s22
        end if
       end do
      end do
      npsf=iprun                                                        7d19s22
      nps=0                                                             7d19s22
      do i=1,iprun                                                      7d19s22
       ii=ibc(jpt1+i)                                                   7d19s22
       nclo=ibasis(1,ii)                                                7d19s22
       nclop=nclo+1
       iarg=nclop-mdon                                                  7d11s19
       nps=nps+ncsf(iarg)                                               7d19s22
       iic=iptrbit(1,nclop,mysym)+ibasis(2,ii)-1                        7d19s22
       iio=iptrbit(2,nclop,mysym)+ibasis(3,ii)-1                        7d19s22
      end do                                                            7d19s22
      if(igcode.eq.0)return                                             7d19s22
      ips=1                                                             6d11s19
      ivv=1                                                             9d4s19
      ivvs=1                                                            9d4s19
      irun=0                                                            7d11s19
      iprun=0                                                           7d11s19
      do i=1,nfcn                                                       6d11s19
       nclo=ibasis(1,i)                                                 7d11s19
       nclop=nclo+1
       iarg=nclop-mdon                                                  7d11s19
       if(ibc(jpt2+i).ne.0)then                                         7d19s22
        do j=1,ncsf(iarg)                                               7d11s19
         irun=irun+1                                                    7d11s19
         iprun=iprun+1                                                  7d11s19
         ipointf(iprun)=irun                                            7d11s19
        end do                                                          7d11s19
        ipointp(ips)=i                                                  3d24s21
        do j=1,3                                                        6d11s19
         ipsbase(j,ips)=ibasis(j,i)                                     6d11s19
        end do                                                          6d11s19
        hdps(ips)=hdiag(i)                                              6d11s19
        do izz=1,ndigs                                                  5d14s21
         zzpsdig(ips,izz)=zzdig(i,izz)                                  5d14s21
        end do                                                          5d14s21
        do j=0,nvv-1                                                    9d4s19
         do k=0,ncsf(iarg)-1                                            9d4s19
          iadf=ivv+k+ncsft*j                                             9d4s19
          iadt=ivvs+j+nvv*k                                             9d4s19
          hivs(iadt)=hiv(iadf)                                          9d4s19
         end do                                                         9d4s19
        end do                                                          9d4s19
        ivvs=ivvs+nvv*ncsf(iarg)                                        9d4s19
        ips=ips+1                                                       6d11s19
       else                                                             7d11s19
        irun=irun+ncsf(iarg)                                            7d11s19
       end if                                                           6d11s19
       ivv=ivv+ncsf(iarg)                                               9d4s19
      end do                                                            6d11s19
      ibcoff=ipt1                                                       7d19s22
      return                                                            6d11s19
      end                                                               6d11s19
