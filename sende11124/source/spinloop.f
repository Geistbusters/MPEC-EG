c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine spinloop(i2sb,i2smb,i2sk,i2smk,nopenb,nopenk,nopenkk,  10d19s21
     $     ncsfb,nrootz,itype,imatx,ntype,nab1,iwpb1,iwpk1,nab2,iwpb2,  10d19s21
     $     iwpk2,veck,ndim,ieoro,bc,ibc)                                11d14s22
      implicit real*8 (a-h,o-z)
      integer*1 nab1(2),nab2(2),ipackc1(8)                              10d19s21
      integer*2 ipack2(4),ipack2m(4)                                    10d19s21
      integer*8 ipackc,ipack,ipackm                                     10d19s21
      logical ldebug                                                    10d28s21
      equivalence (ipackc,ipackc1),(ipack,ipack2),(ipackm,ipack2m)      10d19s21
      dimension mcsf(2),veck(ndim,*)                                    10d19s21
      include "common.store"
      data icall/0/
      save
      icall=icall+1
      ldebug=.false.                                                    3d19s24
      mslow=max(i2smb-2,i2smk-2)                                        10d19s21
      mshi=min(i2smb+2,i2smk+2)                                         10d19s21
      ishi=min(nopenkk,i2sb+2,i2sk+2)                                   10d19s21
      mxmat=4                                                           10d19s21
      itype=ibcoff
      imatx=itype+mxmat
      ibcoff=imatx+mxmat*ncsfb*nrootz                                   10d19s21
      if(ibcoff.lt.0)write(6,*)('enoughf '),ibcoff
      call enough('spinloop.  1',bc,ibc)
      do iz=imatx,ibcoff-1                                              10d19s21
       bc(iz)=0d0                                                       10d19s21
      end do                                                            10d19s21
      ntype=0                                                           10d19s21
      ivec=0                                                            10d28s21
      if(ldebug)then
       write(6,*)('spin for bra: '),i2sb,i2smb,iwpb1,iwpb2
       write(6,*)('spin for ket: '),i2sk,i2smk,iwpk1,iwpk2
       write(6,*)('ms range: '),mslow,mshi,nopenkk
       write(6,*)('ishi: '),ishi
      end if
      do ms=mslow,mshi,2                                                10d19s21
       if(iabs(ms).le.nopenkk)then                                      10d19s21
        jbcoff=ibcoff                                                   10d19s21
        do is=iabs(ms),ishi,2                                           10d19s21
         if(iabs(is-i2sb).le.2.and.iabs(is-i2sk).le.2.and.              10d19s21
     $                iabs(ms-i2smb).le.2.and.iabs(ms-i2smk).le.2)then  10d19s21
          if(ldebug)write(6,*)('for is,ms = '),is,ms
          call gencup(is,ms,i2sk,i2smk,nopenkk,nopenk,nab2,             10d19s21
     $                  iwpb2,iwpk2,ioutg,imatg,ntypeg,mcsf,veck,       10d19s21
     $                  ndim,nrootz,bc,ibc)                             11d14s22
          if(ldebug)write(6,*)('back from gencup ')
          if(ivec.eq.0.and.ldebug)then                                  10d28s21
           write(6,*)('input vectors: ')
           call prntm2(veck,mcsf(2),nrootz,ndim)
           ivec=1                                                       10d28s21
          end if                                                        10d28s21
          nrow=mcsf(1)                                                  10d19s21
          ncol=ntypeg*nrootz                                            10d19s21
          ntype1b=ntypeg                                                10d19s21
          ioutb=ioutg                                                   10d19s21
          imatb=imatg                                                   10d19s21
          call gencup(i2sb,i2smb,is,ms,nopenb,nopenkk,nab1,             10d19s21
     $                  iwpb1,iwpk1,ioutg,imatg,ntypeg,mcsf,            10d19s21
     $                   bc(imatb),nrow,ncol,bc,ibc)                    11d14s22
          if(ldebug)write(6,*)('back from gencup ')
          ntype1a=ntypeg                                                10d19s21
          iouta=ioutg                                                   10d19s21
          do ib=0,ntype1b-1                                             10d19s21
           ipackc=ibc(ioutb+ib)                                         10d19s21
           ipack2(3)=ipackc1(1)                                         10d19s21
           ipack2(4)=ipackc1(2)                                         10d19s21
           do ia=0,ntype1a-1                                            10d19s21
            ipackc=ibc(iouta+ia)                                        10d19s21
            ipack2(1)=ipackc1(1)                                        10d19s21
            ipack2(2)=ipackc1(2)                                        10d19s21
            if(ldebug)write(6,*)('integral type '),ipack2               10d28s21
            do i=1,4                                                    10d19s21
             ipack2m(i)=-ipack2(i)                                      10d19s21
            end do                                                      10d19s21
            mtype=-1                                                    10d19s21
            ff=1d0                                                      10d19s21
            if(ldebug)write(6,*)('ipack,ipackm,ieoro '),ipack,
     $           ipackm,ieoro
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
              stop 'spinloop'                                           10d19s21
             end if                                                     10d19s21
            end if                                                      10d19s21
            itmp=imatx+ncsfb*nrootz*mtype                               10d19s21
            do ir=0,nrootz-1                                            10d19s21
             ifrm=imatg+ncsfb*(ir+nrootz*(ib+ntype1b*ia))               10d20s21
             jtmp=itmp+ncsfb*ir                                         10d19s21
             do i=0,ncsfb-1                                             10d19s21
              orig=bc(jtmp+i)
              bc(jtmp+i)=bc(jtmp+i)+ff*bc(ifrm+i)                       10d19s21
             end do                                                     10d19s21
            end do                                                      10d19s21
            if(ldebug)then
             write(6,*)('we are type '),mtype,ff                            10d19s21
             write(6,*)('running sum')
             call prntm2(bc(itmp),ncsfb,nrootz,ncsfb)
            end if
           end do                                                       10d19s21
          end do
          ibcoff=ioutb                                                  10d19s21
         end if                                                         10d19s21
        end do                                                          10d19s21
        ibcoff=jbcoff                                                   10d19s21
       end if                                                           10d19s21
      end do
      return
      end
