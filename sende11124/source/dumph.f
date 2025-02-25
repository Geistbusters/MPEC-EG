c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dumph(htmp,nr,nc,irow,icol)                            7d17s19
      implicit real*8 (a-h,o-z)                                         7d17s19
      dimension htmp(nr,nc)                                             7d17s19
      do i=1,nc                                                         7d17s19
       jcol=icol+i-1                                                    7d17s19
       do j=1,nr                                                        7d17s19
        jrow=irow+j-1                                                   7d17s19
        if(abs(htmp(j,i)).gt.1d-14)write(12)jrow,jcol,htmp(j,i)         7d17s19
       end do                                                           7d17s19
      end do                                                            7d17s19
      return                                                            7d17s19
      end                                                               7d17s19
