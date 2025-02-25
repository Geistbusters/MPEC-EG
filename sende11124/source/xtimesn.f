c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine xtimesn(ncsfb,ncsfk,ncsfmid,mcol,iwpb,iwpk,xin,ndim,   12d4s20
     $     xout,ndim2,f1,f2,bc,ibc)                                     11d10s22
      implicit real*8 (a-h,o-z)                                         7d17s20
      logical lnew
      dimension xin(ndim,*),xout(ndim2,*)                               12d4s20
      include "common.store"                                            7d17s20
      if(mcol.le.0)return                                               11d26s20
      ibcoffo=ibcoff                                                    2d10s23
c
c     wpb*(wpk*xin)
c
      if(iwpk.gt.0)then                                                 2d10s23
       n1=ncsfmid*ncsfk*mcol                                             2d10s23
      else                                                              2d10s23
       nusedi=ibc(-iwpk)                                                2d10s23
       n1=nusedi*mcol                                                   2d10s23
      end if                                                            2d10s23
      if(iwpb.gt.0)then                                                 2d10s23
       n2=ncsfb*ncsfmid*mcol                                             2d10s23
      else                                                              2d10s23
       nusedi=ibc(-iwpb)                                                2d10s23
       n2=nusedi*mcol                                                   2d10s23
      end if                                                            2d10s23
      n23=n1+n2                                                         2d10s23
c
c     (wpb*wpk)*xin
c
      if(iwpk.lt.0)then                                                 2d10s23
       nusedi=ibc(-iwpk)                                                2d10s23
       n1=nusedi*ncsfb                                                  2d10s23
      else                                                              2d10s23
       if(iwpb.lt.0)then                                                2d10s23
        nusedi=ibc(-iwpb)                                               2d10s23
        n1=nusedi*ncsfk                                                 2d10s23
       else                                                             2d10s23
        n1=ncsfb*ncsfmid*ncsfk                                          2d10s23
       end if                                                           2d10s23
      end if                                                            2d10s23
      n2=ncsfb*ncsfk*mcol                                               2d10s23
      n12=n1+n2                                                         2d10s23
      if(n12.lt.n23)then                                                2d10s23
       iblock=1
       ifirst=ibcoff                                                    2d10s23
       ibcoff=ifirst+ncsfb*ncsfk                                        2d10s23
       call enough('xtimesn.12',bc,ibc)                                 2d10s23
       call prodn(iwpb,iwpk,ncsfb,ncsfk,ncsfmid,bc(ifirst),bc,ibc,      2d13s23
     $      1d0,0d0)                                                    2d13s23
       call dgemm('n','n',ncsfb,mcol,ncsfk,f1,bc(ifirst),ncsfb,         2d10s23
     $      xin,ndim,f2,xout,ndim2)                                     4d19s23
      else                                                              2d10s23
       iblock=2
       ixin=ibcoff                                                      2d13s23
       ifirst=ixin+ncsfk*mcol                                           2d13s23
       ibcoff=ifirst+ncsfmid*mcol                                       2d10s23
       call enough('xtimesn.23',bc,ibc)                                 2d10s23
       jxin=ixin-1                                                      2d13s23
       do i=1,mcol                                                      2d13s23
        do j=1,ncsfk                                                    2d13s23
         bc(jxin+j)=xin(j,i)                                            2d13s23
        end do                                                          2d13s23
        jxin=jxin+ncsfk                                                 2d13s23
       end do                                                           2d13s23
       call prodn(iwpk,ixin,ncsfmid,mcol,ncsfk,bc(ifirst),bc,ibc,f1,    2d13s23
     $      0d0)                                                        2d13s23
       if(ndim2.eq.ncsfb)then                                           4d19s23
        call prodn(iwpb,ifirst,ncsfb,mcol,ncsfmid,xout,bc,ibc,           2d13s23
     $       1d0,f2)                                                     2d13s23
       else                                                             4d19s23
        ixout=ibcoff                                                    4d19s23
        ibcoff=ixout+ncsfb*mcol                                         4d19s23
        call enough('xtimesn.xout',bc,ibc)                              4d19s23
        call prodn(iwpb,ifirst,ncsfb,mcol,ncsfmid,bc(ixout),bc,ibc,           2d13s23
     $       1d0,0d0)                                                   4d19s23
        jxout=ixout-1                                                   4d19s23
        if(f2.eq.0d0)then                                               4d19s23
         do i=1,mcol                                                    4d19s23
          do j=1,ncsfb                                                  4d19s23
           xout(j,i)=bc(jxout+j)                                        4d19s23
          end do                                                        4d19s23
          jxout=jxout+ncsfb                                             4d19s23
         end do                                                         4d19s23
        else                                                            4d19s23
         do i=1,mcol                                                    4d19s23
          do j=1,ncsfb                                                  4d19s23
           xout(j,i)=f2*xout(j,i)+bc(jxout+j)                           4d19s23
          end do                                                        4d19s23
          jxout=jxout+ncsfb                                             4d19s23
         end do                                                         4d19s23
        end if                                                          4d19s23
        ibcoff=ixout                                                    4d19s23
       end if
      end if                                                            2d10s23
      ibcoff=ibcoffo                                                    2d10s23
      return                                                            7d17s20
      end                                                               7d17s20
