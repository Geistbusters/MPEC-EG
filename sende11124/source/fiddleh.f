c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine fiddleh(hdig,nconf,iaorb,nalpha,numa,iborb,nbeta,      8d3s22
     $     numb,ilc,ihc,nsymb,nsbeta,pshdig,iconfig,econfig,norb,mdoo,  8d3s22
     $     nadet,nbdet,nclodet,nlzz,ixlzze,islz,multh,iacto,bc,ibc)     11d14s22
      implicit real*8 (a-h,o-z)
      integer*1 iaorb(nalpha,numa),iborb(nbeta,numb),iboth(32,2)        8d3s22
      integer*8 ipack8,last8,jjc,jjo,iic,iio                            8d3s22
      integer*4 ipack4(2)                                               8d3s22
      integer*2 ipack2(4)                                               8d3s22
      logical lprt                                                      2d27s23
      equivalence (ipack8,ipack4),(ipack8,ipack2)                       8d3s22
      dimension hdig(*),ilc(*),ihc(*),nsbeta(*),iconfig(2,*),econfig(*),8d3s22
     $     nadet(*),nbdet(*),nclodet(*),ixlzze(8,*),islz(*),multh(8,8), 8d3s22
     $     iacto(*),nab4(2,3),iorbqq(32,2)                                           8d3s22
c
c     average det diagonals over all functions with same c:o list
c     for use in chosing ps functions and guaranteeing full spin
c     set by just looking at energies
c     if nlzz ne 0, look at couplings as well so lz2 fcns as well!
c
      include "common.store"
      data icall/0/                                                     8d9s22
      save                                                              8d9s22
      lprt=.true.
      icall=icall+1                                                     8d9s22
      ism=ibcoff                                                        8d3s22
      irel=ism+norb                                                     8d3s22
      ibcoff=irel+norb                                                  8d3s22
      ncount=ibcoff                                                     8d3s22
      ibcoff=ncount+mdoo+1                                              8d3s22
      call enough('fiddleh.  1',bc,ibc)
      jsm=ism-1                                                         8d3s22
      jrel=irel-1                                                       8d3s22
      ioff=0                                                            8d3s22
      do isb=1,nsymb                                                    8d3s22
       do i=0,iacto(isb)-1                                              8d3s22
        ip=i+ioff                                                       8d3s22
        ibc(ism+ip)=isb                                                 8d3s22
        ibc(irel+ip)=i                                                  8d3s22
       end do                                                           8d3s22
       ioff=ioff+iacto(isb)                                             8d3s22
      end do                                                            8d3s22
      do i=0,mdoo                                                       8d3s22
       ibc(ncount+i)=0                                                  8d3s22
      end do                                                            8d3s22
      ioffa=0                                                           8d3s22
      idet=1
      ih=0                                                              8d3s22
      do isa=1,nsymb                                                    8d3s22
       isb=nsbeta(isa)                                                  8d3s22
       ioffb=0                                                          8d3s22
       do i=1,isb-1                                                     8d3s22
        ioffb=ioffb+nbdet(i)                                            8d3s22
       end do                                                           8d3s22
       do iad=1,nadet(isa)                                              8d3s22
        iadp=iad+ioffa
        do i=1,norb                                                     8d3s22
         iboth(i,1)=0                                                   8d3s22
        end do                                                          8d3s22
        do i=1,nalpha                                                   8d3s22
         iboth(iaorb(i,iadp),1)=1                                       8d3s22
        end do                                                          8d3s22
        do ibd=1,nbdet(isb)                                             8d3s22
         if(iad.ge.ilc(isa).and.iad.le.ihc(isa))then                    8d3s22
          ih=ih+1                                                       8d3s22
          huse=hdig(ih)                                                 8d3s22
         else                                                           8d3s22
          huse=0d0                                                      8d3s22
         end if                                                         8d3s22
         ibdp=ibd+ioffb                                                 8d3s22
         do i=1,norb                                                    8d3s22
          iboth(i,2)=iboth(i,1)                                         8d3s22
         end do                                                         8d3s22
         do i=1,nbeta                                                   8d3s22
          iboth(iborb(i,ibdp),2)=iboth(iborb(i,ibdp),2)+1               8d3s22
         end do                                                         8d3s22
         iconfig(1,idet)=0                                              8d3s22
         iconfig(2,idet)=0                                              8d3s22
         econfig(idet)=huse                                             8d3s22
         n2=0                                                           8d3s22
         do i=1,norb                                                    8d3s22
          if(iboth(i,2).eq.2)then                                       8d3s22
           iconfig(1,idet)=ibset(iconfig(1,idet),i)                     8d3s22
           n2=n2+1                                                      8d3s22
          else if(iboth(i,2).eq.1)then                                  8d3s22
           iconfig(2,idet)=ibset(iconfig(2,idet),i)                     8d3s22
          end if                                                        8d3s22
         end do                                                         8d3s22
         nclodet(idet)=n2                                               8d3s22
         ibc(ncount+n2)=ibc(ncount+n2)+1                                8d3s22
         idet=idet+1                                                    8d3s22
        end do                                                          8d3s22
       end do                                                           8d3s22
       ioffa=ioffa+nadet(isa)                                           8d3s22
      end do                                                            8d3s22
      icount=ibcoff                                                     8d3s22
      ibcoff=icount+mdoo+1                                              8d3s22
      call enough('fiddleh.  2',bc,ibc)
      ibase=ibcoff                                                      8d3s22
      ibc(icount)=ibcoff                                                8d3s22
      do i=1,mdoo                                                       8d3s22
       im=i-1                                                            8d3s22
       ibc(icount+i)=ibc(icount+im)+ibc(ncount+im)                      8d3s22
      end do                                                            8d3s22
      ibcoff=ibc(icount+mdoo)+ibc(ncount+mdoo)                          8d3s22
      idpoint=ibcoff                                                    8d3s22
      ibcoff=idpoint+nconf                                              8d3s22
      call enough('fiddleh.  3',bc,ibc)
      idelta=idpoint-ibase                                              8d3s22
      do i=1,idet-1                                                     8d3s22
       ipack4(1)=iconfig(1,i)                                           8d3s22
       ipack4(2)=iconfig(2,i)                                           8d3s22
       iad=ibc(icount+nclodet(i))                                       8d3s22
       ibc(iad)=ipack8                                                  8d3s22
       ibc(iad+idelta)=i                                                8d3s22
       ibc(icount+nclodet(i))=ibc(icount+nclodet(i))+1                  8d3s22
      end do                                                            8d3s22
      ioff=ibase                                                        8d3s22
      iuniq=ibcoff                                                      8d3s22
      ipointq=iuniq+nconf
      ibcoff=ipointq+nconf                                              8d3s22
      call enough('fiddleh.  4',bc,ibc)
      nuniq=0                                                           8d3s22
      last8=0                                                           8d3s22
      do i=0,mdoo                                                       8d3s22
       if(ibc(ncount+i).gt.0)then                                       8d3s22
        isort=ibcoff                                                    8d3s22
        ibcoff=isort+ibc(ncount+i)                                      8d3s22
        call enough('fiddleh.  5',bc,ibc)
        do ii=0,ibc(ncount+i)-1                                          8d3s22
         ipack8=ibc(ibase+ii)                                            8d3s22
        end do
        ndsort=ibc(ncount+i)                                            1d18s23
        call idsortdws(ibc(ibase),ibc(isort),ndsort)
        do ii=0,ibc(ncount+i)-1                                         8d3s22
         jj=ibc(isort+ii)-1                                             8d3s22
         ipack8=ibc(ibase+ii)                                           8d3s22
         if(ibc(ibase+ii).ne.last8)then                                 8d3s22
          last8=ibc(ibase+ii)                                           8d3s22
          ibc(iuniq+nuniq)=last8                                        8d3s22
          nuniq=nuniq+1                                                 8d3s22
         end if                                                         8d3s22
         iad=ibase+jj+idelta                                            8d3s22
         iad=ipointq+ibc(iad)-1                                         8d3s22
         ibc(iad)=nuniq                                                 8d3s22
        end do                                                          8d3s22
        ibase=ibase+ibc(ncount+i)                                       8d3s22
       end if                                                           8d3s22
      end do
      muniq=ibcoff                                                      8d3s22
      ieuniq=muniq+nuniq                                                8d3s22
      ibcoff=ieuniq+nuniq                                               8d3s22
      call enough('fiddleh.  6',bc,ibc)
      do iz=muniq,ieuniq-1                                              8d3s22
       ibc(iz)=0                                                        8d3s22
      end do                                                            8d3s22
      do iz=ieuniq,ibcoff-1                                             8d3s22
       bc(iz)=0d0                                                       8d3s22
      end do                                                            8d3s22
      do i=1,nconf                                                      8d3s22
       im=i-1                                                           8d3s22
       iad=muniq+ibc(ipointq+im)-1                                      8d3s22
       ibc(iad)=ibc(iad)+1                                              8d3s22
       if(econfig(i).ne.0d0)then                                        8d3s22
        iad=ieuniq+ibc(ipointq+im)-1                                    8d3s22
        bc(iad)=bc(iad)+econfig(i)                                      8d3s22
       end if                                                           8d3s22
      end do                                                            8d3s22
      call dws_gsumf(bc(ieuniq),nuniq)                                  8d3s22
      do i=1,nuniq                                                      8d3s22
       im=i-1                                                           8d3s22
       ipack8=ibc(iuniq+im)
       eavg=bc(ieuniq+im)/dfloat(ibc(muniq+im))                         8d3s22
       bc(ieuniq+im)=eavg                                               8d3s22
       ipack8=ibc(iuniq+im)
      end do                                                            8d3s22
      if(nlzz.ne.0)then                                                 8d3s22
       if(nlzz.eq.6)then                                                2d27s23
        npass=3                                                         2d27s23
       else                                                             2d27s23
        npass=1                                                         2d27s23
       end if                                                           2d27s23
       nec=nalpha+nbeta                                                  8d3s22
       ncoup=0                                                          8d3s22
       icoup=ibcoff                                                     8d3s22
       jcoup=icoup                                                      8d3s22
       do if=1,nuniq-1                                                  8d3s22
        ifm=if-1                                                        8d3s22
        ipack8=ibc(iuniq+ifm)                                           8d3s22
        iic=ipack4(1)                                                   8d3s22
        ncloi=popcnt(iic)                                                8d3s22
        nopeni=nec-2*ncloi                                              8d3s22
        iio=ipack4(2)
        do jf=if+1,nuniq                                                8d3s22
         jfm=jf-1                                                       8d3s22
         ipack8=ibc(iuniq+jfm)                                           8d3s22
         jjc=ipack4(1)                                                   8d3s22
         jjo=ipack4(2)
         ncloj=popcnt(jjc)                                              8d3s22
         if(iabs(ncloj-ncloi).le.2)then                                 8d3s22
          nopenj=nec-2*ncloj                                            8d3s22
          call gandc4(jjc,jjo,iic,iio,nopenj,nopeni,norb,nnot,nab4,bc,  11d14s22
     $         ibc)                                                     11d14s22
          if(nnot.ge.3)then                                             8d3s22
           is12=multh(ibc(jsm+nab4(1,1)),ibc(jsm+nab4(2,1)))            8d3s22
           is34=multh(ibc(jsm+nab4(1,2)),ibc(jsm+nab4(2,2)))
           xcoup=0d0                                                    8d9s22
           do ixyz=1,npass                                              2d27s23
            if(is12.eq.islz(ixyz).and.is34.eq.islz(ixyz))then           2d27s23
             iad1=ixlzze(ibc(jsm+nab4(1,1)),ixyz)+ibc(jrel+nab4(1,1))   2d27s23
     $          +iacto(ibc(jsm+nab4(1,1)))*ibc(jrel+nab4(2,1))          8d3s22
            iad2=ixlzze(ibc(jsm+nab4(1,2)),ixyz)+ibc(jrel+nab4(1,2))    2d27s23
     $          +iacto(ibc(jsm+nab4(1,2)))*ibc(jrel+nab4(2,2))          8d3s22
             xcoup=xcoup+abs(bc(iad1)*bc(iad2))                         2d27s23
            end if                                                      2d27s23
           end do                                                       2d27s23
            if(abs(xcoup).gt.1d-10)then                                 8d3s22
             ipack4(1)=if                                               8d3s22
             ipack4(2)=jf                                               8d3s22
             ibc(jcoup)=ipack8                                          8d3s22
             jcoup=jcoup+1                                              8d3s22
            end if                                                      8d3s22
           if(nnot.eq.4.and.abs(xcoup).lt.1d-10)then                    8d9s22
            is14=multh(ibc(jsm+nab4(1,1)),ibc(jsm+nab4(2,2)))            8d3s22
            is32=multh(ibc(jsm+nab4(1,2)),ibc(jsm+nab4(2,1)))
            do ixyz=1,npass                                             2d27s23
             if(is14.eq.islz(ixyz).and.is32.eq.islz(ixyz))then          2d27s23
              iad1=ixlzze(ibc(jsm+nab4(1,1)),ixyz)+ibc(jrel+nab4(1,1))  2d27s23
     $          +iacto(ibc(jsm+nab4(1,1)))*ibc(jrel+nab4(2,2))          8d3s22
              iad2=ixlzze(ibc(jsm+nab4(1,2)),ixyz)+ibc(jrel+nab4(1,2))  2d27s23
     $          +iacto(ibc(jsm+nab4(1,2)))*ibc(jrel+nab4(2,1))          8d9s22
              xcoup=xcoup+abs(bc(iad1)*bc(iad2))                        2d27s23
             end if                                                     2d27s23
            end do                                                      2d27s23
             if(abs(xcoup).gt.1d-10)then                                 8d3s22
              ipack4(1)=if                                               8d3s22
              ipack4(2)=jf                                               8d3s22
              ibc(jcoup)=ipack8                                          8d3s22
              jcoup=jcoup+1                                              8d3s22
             end if                                                      8d3s22
           end if                                                       8d9s22
          end if                                                        8d3s22
         end if                                                         8d3s22
        end do                                                          8d3s22
       end do                                                           8d3s22
       ncoup=jcoup-icoup                                                 8d3s22
       ibcoff=jcoup                                                      8d3s22
       ihit=ibcoff                                                       8d3s22
       igroup=ihit+nuniq                                                 8d3s22
       ibcoff=igroup+nuniq                                               8d3s22
       call enough('fiddleh.  7',bc,ibc)
       jgroup=igroup-1                                                   8d3s22
       jhit=ihit-1
       do i=0,nuniq-1                                                    8d3s22
        ibc(ihit+i)=0                                                    8d3s22
       end do                                                            8d3s22
   57  continue                                                          8d3s22
        ngroup=0                                                         8d3s22
        do i=0,nuniq-1                                                   8d3s22
         if(ibc(ihit+i).eq.0)then                                        8d3s22
          ibc(igroup)=i+1                                                8d3s22
          ibc(ihit+i)=1                                                  8d3s22
          ngroup=1                                                       8d3s22
          go to 58                                                       8d3s22
         end if                                                          8d3s22
        end do                                                           8d3s22
        go to 59                                                         8d3s22
   58   continue                                                         8d3s22
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
         im=in-1                                                         8d3s22
         eavg=eavg+bc(ieuniq+im)                                         8d3s22
        end do                                                           8d3s22
        eavg=eavg/dfloat(ngroup)                                         8d3s22
        do i=0,ngroup-1                                                  8d3s22
         in=ibc(igroup+i)                                                8d3s22
         im=in-1                                                         8d3s22
         bc(ieuniq+im)=eavg                                              8d3s22
        end do                                                           8d3s22
        go to 57                                                          8d3s22
   59  continue                                                          8d3s22
      end if                                                            8d3s22
      ioffa=0                                                           8d3s22
      idet=1
      ih=0                                                              8d3s22
      do isa=1,nsymb                                                    8d3s22
       isb=nsbeta(isa)                                                  8d3s22
       ioffb=0                                                          8d3s22
       do i=1,isb-1                                                     8d3s22
        ioffb=ioffb+nbdet(i)                                            8d3s22
       end do                                                           8d3s22
       do iad=1,nadet(isa)                                              8d3s22
        iadp=iad+ioffa
        do ibd=1,nbdet(isb)                                             8d3s22
         ibdp=ibd+ioffb                                                 8d3s22
         ipack2(1)=iad                                                  8d3s22
         ipack2(2)=ibd                                                  8d3s22
         ipack2(3)=isa                                                  8d3s22
         im=idet-1                                                      8d3s22
         iadd=ieuniq+ibc(ipointq+im)-1                                     8d3s22
         econfig(idet)=bc(iadd)                                         8d3s22
         iconfig(1,idet)=ipack4(1)                                      8d3s22
         iconfig(2,idet)=ipack4(2)                                      8d3s22
         idet=idet+1                                                    8d3s22
        end do                                                          8d3s22
       end do                                                           8d3s22
      end do                                                            8d3s22
      ibcoff=ism                                                        8d3s22
      return                                                            8d3s22
      end                                                               8d3s22
