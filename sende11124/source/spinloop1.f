c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine spinloop1(i2sb,i2smb,i2sk,i2smk,nopen,ncsfb,imap,      10d20s21
     $     nrootz,itype,imatx,ntype,irw0,veck,ndim,ieoro,bc,ibc)        11d14s22
      implicit real*8 (a-h,o-z)                                         10d20s21
      include "common.store"                                            10d20s21
      integer*1 imap(*),ipackc1(8)                                      10d20s21
      integer*2 ipack2(4),ipack2m(4)                                    10d20s21
      integer*8 ipackc,ipack,ipackm                                     10d20s21
      equivalence (ipackc,ipackc1),(ipack,ipack2),(ipackm,ipack2m)      10d20s21
      dimension mcsf(2),veck(ndim,*)                                    10d20s21
      mslow=max(i2smb-2,i2smk-2)                                        10d20s21
      mshi=min(i2smb+2,i2smk+2)                                         10d20s21
      ishi=min(nopen,i2sb+2,i2sk+2)                                     10d20s21
      mxmat=nopen*nopen*4                                               10d20s21
      ntype=0                                                           10d20s21
      itype=ibcoff                                                      10d19s21
      imatx=itype+mxmat                                                  10d19s21
      ibcoff=imatx+mxmat*ncsfb*nrootz                                    10d19s21
      if(ibcoff.lt.0)write(6,*)('enoughe '),ibcoff
      call enough('spinloop1.  1',bc,ibc)
      do ms=mslow,mshi,2                                                10d20s21
       if(iabs(ms).le.nopen)then                                        10d20s21
        do is=iabs(ms),ishi,2                                           10d20s21
         call getcup(irw0,nopen,i2sb,i2smb,is,ms,ntypeg,                10d20s21
     $                 ioutg,imatg,mcsf,imap,bc,ibc)                    11d14s22
         ntype1a=ntypeg                                                 10d20s21
         iouta=ioutg                                                    10d20s21
         imata=imatg                                                    10d20s21
         ncsfkk=mcsf(2)                                                 10d20s21
         call getcup(irw0,nopen,is,ms,i2sk,i2smk,ntypeg,                10d20s21
     $                 ioutg,imatg,mcsf,imap,bc,ibc)                    11d14s22
         ncsfk=mcsf(2)                                                  10d19s21
         ntype1b=ntypeg                                                 10d20s21
         ioutb=ioutg                                                    10d20s21
         imatb=imatg                                                    10d20s21
         if(ncsfkk.gt.0)then                                            10d20s21
          itmpb=ibcoff                                                  10d20s21
          ibcoff=itmpb+ncsfkk*nrootz                                    10d20s21
          call enough('spinloop1.  2',bc,ibc)
          do ib=0,ntype1b-1                                             10d20s21
           ipackc=ibc(ioutb+ib)                                         10d20s21
           ipack2(3)=ipackc1(1)                                         10d20s21
           ipack2(4)=ipackc1(2)                                         10d20s21
           call dgemm('n','n',ncsfkk,nrootz,ncsfk,1d0,                  10d19s21
     $          bc(imatb),ncsfkk,veck,ndim,0d0,bc(itmpb),ncsfkk,        10d20s21
     d' spinloop1.  1')
           jmata=imata                                                  10d20s21
           do ia=0,ntype1a-1                                            10d20s21
            ipackc=ibc(iouta+ia)                                        10d20s21
            ipack2(1)=ipackc1(1)                                        10d20s21
            ipack2(2)=ipackc1(2)                                        10d20s21
            do i=1,4                                                    10d19s21
             ipack2m(i)=-ipack2(i)                                      10d19s21
            end do                                                      10d19s21
            mtype=-1                                                    10d19s21
            ff=1d0                                                      10d19s21
            do i=0,ntype-1                                              10d19s21
             if(ipack.eq.ibc(itype+i))then                              10d19s21
              mtype=i                                                   10d19s21
              fact=1d0                                                  10d19s21
             else if(ipackm.eq.ibc(itype+i))then                        10d19s21
              mtype=i                                                   10d19s21
              fact=1d0                                                  10d19s21
              if(ieoro.eq.1)ff=-1d0                                     10d19s21
             end if                                                     10d19s21
            end do                                                      10d19s21
            if(mtype.lt.0)then                                          10d19s21
             ibc(itype+ntype)=ipack                                     10d19s21
             fact=0d0                                                   10d19s21
             mtype=ntype                                                10d19s21
             ntype=ntype+1                                              10d19s21
             if(ntype.gt.mxmat)then                                     10d19s21
              write(6,*)('tooooo many integral types!! ')               10d19s21
              do i=0,ntype-1                                            10d19s21
               ipack=ibc(itype+i)                                       10d19s21
               write(6,*)i,('type '),ipack2                             10d19s21
              end do                                                    10d19s21
              call dws_synca                                            10d19s21
              call dws_finalize                                         10d19s21
              stop 'spinloop1'                                           10d19s21
             end if                                                     10d19s21
            end if                                                      10d19s21
            itmp=imatx+ncsfb*nrootz*mtype                               10d19s21
            call dgemm('n','n',ncsfb,nrootz,ncsfkk,ff,bc(jmata),ncsfb,  10d20s21
     $           bc(itmpb),ncsfkk,fact,bc(itmp),ncsfb,                  10d20s21
     d' spinloop1.  2')
            jmata=jmata+ncsfb*ncsfkk                                    10d20s21
           end do                                                       10d20s21
           imatb=imatb+ncsfk*ncsfkk                                     10d20s21
          end do
          ibcoff=itmpb                                                  10d20s21
         end if
         ibcoff=iouta                                                   10d20s21
        end do                                                          10d20s21
       end if                                                           10d20s21
      end do                                                            10d20s21
      return                                                            10d20s21
      end                                                               10d20s21
