c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hsscsf12(nhand,ih0av,nh0av,iff1,nff1,ncsf,mdon,mdoo,   11d19s20
     $     ioooo,jmats,kmats,nec,nvirt,shift,isbv,iffoff,bc,ibc)        11d10s22
      implicit real*8 (a-h,o-z)                                         7d15s19
      external second                                                   11d1s19
      integer*8 nhand(*),i1,i12,in                                      12d7s20
      dimension ih0av(*),nh0av(*),ioooo(*),jmats(*),kmats(*),ncsf(*),   11d19s20
     $     iff1(*),nff1(*),idorb(32),isorb(32)                          11d19s20
      common/fnd2cm/inv(2,8,8,8)                                        4d9s18
      common/kmfind/invk1(2,8,8,8,2)                                    6d30s18
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      include "common.mrci"                                             7d15s19
      include "common.store"                                            7d15s19
      call ilimts(nvirt,nvirt,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e) 7d15s19
      nhere=ih+1-il                                                     7d15s19
      loop=0
      i1=1                                                              11d19s20
      do nclop=1,mdoo+1                                                 11d19s20
       nclo=nclop-1                                                     11d19s20
       if(nff1(nclop).gt.0)then                                         11d19s20
        nopen=nec-nclo*2-1                                              11d19s20
        iarg=nclop-mdon                                                 11d19s20
        igg=ibcoff                                                       11d19s20
        ibcoff=igg+nff1(nclop)*nvirt                                    11d19s20
        call enough('hsscsf12.  1',bc,ibc)
        do i=igg,ibcoff-1                                               11d19s20
         bc(i)=0d0                                                      11d19s20
        end do                                                          11d19s20
        in=nff1(nclop)                                                  12d7s20
        jgg=igg                                                         11d19s20
        do if=1,nff1(nclop)                                             11d19s20
         ncheck=popcnt(iff1(iffoff))                                    11d19s20
         if(ncheck.ne.nclo)stop 'ncheckc'
         ii=1                                                           11d19s20
         do i=1,norb                                                    11d19s20
          if(btest(iff1(iffoff),i))then                                 11d19s20
           idorb(ii)=i                                                  11d19s20
           ii=ii+1                                                      11d19s20
          end if                                                        11d19s20
         end do                                                         11d19s20
         iffoff=iffoff+1                                                11d19s20
         ii=1                                                           11d19s20
         ncheck=popcnt(iff1(iffoff))
         if(ncheck.ne.nopen)stop 'nhecko'
         do i=1,norb                                                    11d19s20
          if(btest(iff1(iffoff),i))then                                 11d19s20
           isorb(ii)=i                                                  11d19s20
           ii=ii+1                                                      11d19s20
          end if                                                        11d19s20
         end do                                                         11d19s20
         iffoff=iffoff+1                                                11d19s20
         sumc=0d0                                                         7d15s19
         do i=1,nclo                                                    11d19s20
          ig=irel(idorb(i))-1                                           11d19s20
          is=ism(idorb(i))                                              11d19s20
          iad=ih0av(is)+ig*(nh0av(is)+1)                                  7d15s19
          h0=bc(iad)                                                      7d15s19
          do j=1,nclo                                                   11d19s20
           jg=irel(idorb(j))-1                                          11d19s20
           js=ism(idorb(j))                                             11d19s20
           xj=getint(ioooo,is,is,js,js,ig,ig,jg,jg,bc,ibc)              11d15s22
           xk=getint(ioooo,is,js,js,is,ig,jg,jg,ig,bc,ibc)              11d15s22
           sumc=sumc+2d0*xj-xk                                            7d15s19
          end do                                                          7d15s19
          sumc=sumc+h0*2d0                                                7d15s19
         end do                                                           7d15s19
         sumo=0d0                                                         7d15s19
         do i=1,nopen                                                   11d19s20
          ig=irel(isorb(i))-1                                           11d19s20
          is=ism(isorb(i))                                              11d19s20
          iad=ih0av(is)+ig*(nh0av(is)+1)                                  7d15s19
          h0=bc(iad)                                                      7d15s19
          sumo=sumo+h0                                                    7d15s19
          do j=1,nopen                                                  11d19s20
           if(i.ne.j)then                                                 7d15s19
            jg=irel(isorb(j))-1                                         11d19s20
            js=ism(isorb(j))                                            11d19s20
            xj=getint(ioooo,is,is,js,js,ig,ig,jg,jg,bc,ibc)             11d15s22
            sumo=sumo+xj*0.5d0                                            7d15s19
           end if                                                         7d15s19
          end do                                                          7d15s19
         end do                                                           7d15s19
         sum=sumc+sumo+shift                                              7d15s19
         do i=1,nclo                                                    11d19s20
          is=ism(idorb(i))                                              11d19s20
          ig=irel(idorb(i))-1                                           11d19s20
          do j=1,nopen                                                  11d19s20
           js=ism(isorb(j))                                             11d19s20
           jg=irel(isorb(j))-1                                          11d19s20
           xj=getint(ioooo,is,is,js,js,ig,ig,jg,jg,bc,ibc)              11d15s22
           xk=getint(ioooo,is,js,js,is,ig,jg,jg,ig,bc,ibc)              11d15s22
           sum=sum+xj*2d0-xk                                              7d15s19
          end do                                                          7d15s19
         end do                                                           7d15s19
         do iv=0,nvirt-1                                                  7d15s19
          icol=iv*(nvirt+1)+1                                             7d15s19
          if(icol.ge.il.and.icol.le.ih)then                               7d15s19
           ivp=iv+irefo(isbv)                                              7d15s19
           icolp=icol-il                                                  7d15s19
           iad=ih0av(isbv)+ivp*(nh0av(isbv)+1)                            7d15s19
           h0=bc(iad)                                                     7d15s19
           here=h0+sum                                                  11d19s20
           do j=1,nopen                                                 11d19s20
            jg=irel(isorb(j))-1                                         11d19s20
            js=ism(isorb(j))                                            11d19s20
            i2eu=inv(1,js,js,isbv)                                         7d15s19
            irow=((jg*(jg+1))/2)+jg                                       7d15s19
            iad=jmats(i2eu)+icolp+nhere*irow                              7d15s19
            here=here+bc(iad)                                             7d15s19
           end do                                                         7d15s19
           do j=1,nclo                                                  11d19s20
            jg=irel(idorb(j))-1                                         11d19s20
            js=ism(idorb(j))                                            11d19s20
            i2eu=inv(1,js,js,isbv)                                        7d15s19
            irow=((jg*(jg+1))/2)+jg                                       7d15s19
            iadj=jmats(i2eu)+icolp+nhere*irow                             7d15s19
            i2eu=invk1(1,js,js,isbv,1)                                    7d15s19
            irow=jg+irefo(js)*jg                                          7d15s19
            iadk=kmats(i2eu)+icolp+nhere*irow                             7d15s19
            here=here+2d0*bc(iadj)-bc(iadk)                               7d15s19
           end do                                                         7d15s19
           bc(jgg+iv)=here                                              11d19s20
          end if                                                        11d19s20
         end do                                                          11d19s20
         jgg=jgg+nvirt                                                  11d19s20
        end do                                                          11d19s20
        i12=nvirt                                                       12d7s20
        call ddi_acc(bc,ibc,nhand(nclop),i1,i12,i1,in,bc(igg))          11d15s22
        ibcoff=igg                                                       11d19s20
       end if                                                           11d19s20
      end do                                                            11d19s20
      return                                                            7d15s19
      end                                                               7d15s19
