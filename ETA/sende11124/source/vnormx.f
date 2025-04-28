c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine vnormx(vin,vout,nrow,ncol,ldebug)                      1d25s21
      implicit real*8 (a-h,o-z)
      logical ldebug                                                    1d25s21
      dimension vin(nrow,ncol),vout(nrow,ncol)
      do i=1,ncol
       do j=1,nrow
        vout(j,i)=vin(j,i)
       end do
      end do
      if(ldebug)then                                                    1d25s21
       write(6,*)('starting vectors')
       call prntm2(vin,nrow,ncol,nrow)
      end if                                                            1d25s21
      do i=1,ncol
       do j=1,i-1
        dot=0d0
        do k=1,nrow
         dot=dot+vout(k,j)*vout(k,i)
        end do
        do k=1,nrow
         vout(k,i)=vout(k,i)-dot*vout(k,j)
        end do
       end do
       dot=0d0
       do k=1,nrow
        dot=dot+vout(k,i)**2
       end do
       if(ldebug)write(6,*)('sz of vector '),i,('is '),dot              1d25s21
       xnorm=1d0/sqrt(dot)
       do k=1,nrow
        vout(k,i)=vout(k,i)*xnorm
       end do
      end do
      if(ldebug)then                                                    1d25s21
       write(6,*)('final vectors ')
       call prntm2(vout,nrow,ncol,nrow)
      end if                                                            1d25s21
      return
      end
