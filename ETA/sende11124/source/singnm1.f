c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine singnm1(ham,nvrt,vec,nadet,nbdet,ma,idata,mb,idatb,ism,5d7s18
     $     irelo,jmats,jtype,kmats,ktype,il,ih,nvirtc,nrj,nrk,noc,iroot,5d10s18
     $     nroot,wwww,bc,ibc)                                           11d10s22
      implicit real*8 (a-h,o-z)
      integer*8 ism(*),irelo(*)
      integer*2 idata(4,ma),idatb(4,mb),idwstmp2                        5d7s18
      integer*1 idwstmp1(2)
      equivalence (idwstmp2,idwstmp1)                                   5d7s18
      dimension ham(*),vec(nroot,nadet,*),phs(2),jmats(*),kmats(*),     5d10s18
     $     jtype(*),wwww(*),                                            1d31s22
     $     ktype(*),nrj(*),nrk(*),noc(*)                                       5d7s18
      include "common.store"
      data phs/1d0,-1d0/                                                3d10s17
      if(ma*mb.eq.0)return                                              5d7s18
      idwstmp2=idatb(3,1)
      idwstmp2=idata(3,1)
      do ib=1,mb                                                        5d7s18
       idwstmp2=idatb(3,ib)                                             5d7s18
       ito=idatb(1,ib)                                                  5d7s18
       ifr=idatb(2,ib)                                                  5d7s18
       nto=irelo(idwstmp1(1))                                           5d7s18
       isb=ism(idwstmp1(1))                                             5d7s18
       nfr=irelo(idwstmp1(2))                                           5d7s18
       ix=max(nto,nfr)                                                  5d7s18
       in=min(nto,nfr)                                                  5d7s18
       jj=jmats(jtype(isb))+in-1+((ix*(ix-1))/2)                        5d8s18
       p=phs(idatb(4,ib))                                               5d7s18
       do iva=0,nvrt-1                                                  5d7s18
        do ivb=0,iva                                                    5d7s18
         icol=ivb+1+nvirtc*iva                                          5d7s18
         if(icol.ge.il.and.icol.le.ih)then                              5d7s18
          jju=jj+nrj(isb)*(icol-il)                                     5d7s18
          trm=p*bc(jju)*2d0                                             5d8s18
          iadh=((iva*(iva+1))/2)+ivb+1                                  5d7s18
          do ia=1,nadet                                                 5d7s18
           ham(iadh)=ham(iadh)+vec(iroot,ia,ito)*vec(iroot,ia,ifr)*trm  5d10s18
     $          *wwww(iroot)                                            1d31s22
          end do                                                        5d7s18
         end if                                                         5d7s18
        end do                                                          5d7s18
       end do                                                           5d7s18
      end do                                                            5d7s18
      do ia=1,ma                                                        5d7s18
       idwstmp2=idata(3,ia)                                             5d7s18
       ito=idata(1,ia)                                                  5d7s18
       ifr=idata(2,ia)                                                  5d7s18
       nto=irelo(idwstmp1(1))                                           5d7s18
       isb=ism(idwstmp1(1))                                             5d7s18
       nfr=irelo(idwstmp1(2))                                           5d7s18
       ix=max(nto,nfr)                                                  5d7s18
       in=min(nto,nfr)                                                  5d7s18
       jj=jmats(jtype(isb))+in-1+((ix*(ix-1))/2)                        5d8s18
c
c     or is it nfr,nto?
       kk=kmats(ktype(isb))+nto-1+noc(isb)*(nfr-1)                      5d8s18
       kkp=kmats(ktype(isb))+nfr-1+noc(isb)*(nto-1)                     5d8s18
       p=phs(idata(4,ia))                                               5d7s18
       do iva=0,nvrt-1                                                  5d7s18
        do ivb=0,iva                                                    5d7s18
         icol=ivb+1+nvirtc*iva                                          5d7s18
         if(icol.ge.il.and.icol.le.ih)then                              5d7s18
          jju=jj+nrj(isb)*(icol-il)                                     5d7s18
          kku=kk+nrk(isb)*(icol-il)                                     5d7s18
          kkpu=kkp+nrk(isb)*(icol-il)                                   5d8s18
          trm=p*(2d0*bc(jju)-bc(kku)-bc(kkpu))                          5d8s18
          iadh=((iva*(iva+1))/2)+ivb+1                                  5d7s18
          do ib=1,nbdet                                                 5d7s18
           ham(iadh)=ham(iadh)+vec(iroot,ito,ib)*vec(iroot,ifr,ib)*trm  5d10s18
     $          *wwww(iroot)                                            1d31s22
          end do                                                        5d7s18
         end if                                                         5d7s18
        end do                                                          5d7s18
       end do                                                           5d7s18
      end do                                                            5d7s18
      return                                                            5d7s18
      end                                                               5d7s18
