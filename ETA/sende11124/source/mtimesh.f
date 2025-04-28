c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine mtimesh(xin,xout,istart,nhere,nroot,nps,ham,lsym,v,    4d10s23
     $     e,nwds)                                                      4d10s23
      implicit real*8 (a-h,o-z)                                         4d10s23
      logical lsym                                                      4d10s23
      dimension xin(*),xout(*),ham(*),v(nps,*),e(*)                     4d10s23
      do i=1,nwds                                                       4d10s23
       xout(i)=0d0                                                      4d10s23
      end do                                                            4d10s23
      call dgemm('n','n',nhere,nroot,nps,1d0,ham,nhere,xin,nps,0d0,     4d10s23
     $     xout(istart),nps,'mtimesh.1')                                4d10s23
      call dws_gsumf(xout,nwds)                                         4d10s23
      if(lsym)then                                                      4d10s23
       do ir=1,nroot                                                    4d10s23
        dot=0d0                                                         4d10s23
        iad=nps*(ir-1)                                                  4d10s23
        do i=1,nps                                                      4d10s23
         dot=dot+xin(iad+i)*v(i,ir)                                     4d10s23
        end do                                                          4d10s23
        do i=1,nps                                                      4d10s23
         xout(iad+i)=xout(iad+i)+e(ir)*(dot*v(i,ir)-xin(iad+i))         4d10s23
        end do                                                          4d10s23
       end do                                                           4d10s23
      else                                                              4d10s23
       do ir=1,nroot                                                    4d10s23
        iad=nps*(ir-1)                                                  4d10s23
        do i=1,nps                                                      4d10s23
         xout(iad+i)=xout(iad+i)-e(ir)*xin(iad+i)                       4d10s23
        end do                                                          4d10s23
       end do                                                           4d10s23
      end if                                                            4d10s23
      return                                                            4d10s23
      end                                                               4d10s23
