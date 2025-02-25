c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine checkps(vec,ncsft,nrootz,nps,ipointf,iok,bc,ibc)       11d14s22
      implicit real*8 (a-h,o-z)                                         2d25s21
      dimension vec(ncsft,*),ipointf(*)                                 3d25s21
      include "common.store"
      ivcpy=ibcoff                                                      2d25s21
      ibcoff=ivcpy+ncsft                                                2d25s21
      jvcpy=ivcpy-1                                                     2d25s21
      call enough('checkps.  1',bc,ibc)
      iok=0                                                             2d25s21
      do ir=1,nrootz
       dot=0d0                                                          2d25s21
       do j=1,ncsft                                                     2d25s21
        dot=dot+vec(j,ir)**2                                            2d25s21
        bc(jvcpy+j)=vec(j,ir)                                           2d25s21
       end do                                                           2d25s21
       do j=1,nps                                                       2d25s21
        jj=ipointf(j)                                                   2d25s21
        bc(jvcpy+jj)=0d0                                                2d25s21
       end do                                                           2d25s21
       dota=0d0                                                         2d25s21
       do j=1,ncsft                                                     2d25s21
        dota=dota+bc(jvcpy+j)**2                                        2d25s21
       end do                                                           2d25s21
       if(dota.gt.1d-2)iok=iok+1                                        2d25s21
      end do                                                            2d25s21
      ibcoff=ivcpy                                                      2d25s21
      return
      end
