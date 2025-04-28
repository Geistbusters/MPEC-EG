c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine setbqn(ngaus,ibdat,ibstor,isstor,iorbsym,iorbsymz,     4d19s21
     $     idorel,nlzz,nsymb,nbasisp,nbasisc,iapair,bc,ibc)             11d9s22
      implicit real*8 (a-h,o-z)                                         4d19s21
      integer*8 ibstor(*),isstor(*)                                     4d19s21
      include "common.store"                                            4d19s21
      include "common.spher"                                            4d19s21
      dimension iorbsym(*),iorbsymz(*),nbasisp(*),nbasisc(*),iapair(3,*)9d30s21
c
c     set lz and perhaps l^2 qs for primative and then contracted basis 4d19s21
c     functions.                                                        4d19s21
c
c
c     for primative fcns
c
      jbdat=ibdat-1                                                     4d19s21
      kbdat=jbdat+3*ngaus                                               4d19s21
      n7=ngaus*7                                                        9d30s21
      do ig=1,ngaus                                                     4d19s21
       lhere=ibc(jbdat+ig)                                              4d19s21
       ig8=jbdat+ig+n7                                                  9d30s21
       lp=lhere+1                                                       4d19s21
       ioff=ibc(kbdat+ig)                                               4d19s21
       nl=2*lhere+1                                                     4d19s21
       iioff=0                                                          9d30s21
       do ipass=1,max(1,iapair(1,ibc(ig8)))                             10d1s21
        do i=1,nl                                                        4d19s21
         ip=i+iioff                                                     9d30s21
         iad=iorbsym(isstor(ioff+ip))+ibstor(ioff+ip)-1                 9d30s21
         if(nlzz.eq.2.and.nsymb.eq.8)then                               10d1s21
          if(isstor(ioff+ip).eq.2.or.isstor(ioff+ip).eq.3.or.
     $       isstor(ioff+ip).eq.5.or.isstor(ioff+ip).eq.8)then          10d1s21
           ibc(iad)=mlzb(i,lp)+50                                       10d1s21
          else                                                          10d1s21
           ibc(iad)=mlzb(i,lp)                                          10d1s21
          end if                                                        10d1s21
         else                                                           10d1s21
          ibc(iad)=mlzb(i,lp)                                             4d19s21
         end if                                                         10d1s21
        end do                                                           4d19s21
        iioff=iioff+nl                                                  9d30s21
       end do                                                           9d30s21
       if(nlzz.ne.2)then                                                4d19s21
        do i=1,nl                                                        4d19s21
         iad=iorbsymz(isstor(ioff+i))+ibstor(ioff+i)-1                  4d19s21
         ibc(iad)=lhere                                                 4d19s21
        end do                                                          4d19s21
       end if                                                           4d19s21
      end do                                                            4d19s21
      if(nlzz.eq.2)then                                                 4d19s21
       do isb=1,nsymb
        do i=1,nbasisp(isb)
         iad=iorbsym(isb)+i-1
        end do
        if(idorel.ne.0)then                                             4d19s21
         do i=1,nbasisp(isb)
          iad=iorbsym(isb)+i-1
          iadp=iad+nbasisp(isb)                                         4d19s21
          ibc(iadp)=ibc(iad)                                            4d19s21
         end do
        end if                                                          4d19s21
       end do                                                            4d19s21
      else                                                              4d19s21
       do isb=1,nsymb
        do i=1,nbasisp(isb)
         iad=iorbsym(isb)+i-1
         iadz=iorbsymz(isb)+i-1
        end do
        if(idorel.ne.0)then                                             4d19s21
         do i=1,nbasisp(isb)
          iad=iorbsym(isb)+i-1
          iadp=iad+nbasisp(isb)                                         4d19s21
          ibc(iadp)=ibc(iad)                                            4d19s21
          iad=iorbsymz(isb)+i-1
          iadp=iad+nbasisp(isb)                                         4d19s21
          ibc(iadp)=ibc(iad)                                            4d19s21
         end do
        end if                                                          4d19s21
       end do                                                            4d19s21
      end if                                                            4d19s21
      return                                                            4d19s21
      end                                                               4d19s21
