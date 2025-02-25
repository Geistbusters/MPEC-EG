c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine cart2spher(cart,spher,nrow,lorb,bc,ibc)                11d14s22
      implicit real*8 (a-h,o-z)
      dimension cart(nrow,1),spher(nrow,1)
      include "common.spher"                                            2d19s10
      include "common.store"                                            2d19s10
      idat1=ippn(lorb+1)
      jj=idat1+1
      do i=1,ibc(idat1)
       nt=ibc(jj)
       jj0=jj
       jj=jj+1
       nz=ibc(jj)
       jj=jj+1
       kk=ibc(jj)
       jj=jj+1
       do j=1,nz
        ito=ibc(jj+nt)
        ifm=ibc(jj)
        if(nt.eq.1)then
         do l=1,nrow
          spher(l,ito)=bc(kk)*cart(l,ifm)
         end do
         jj=jj+1
         kk=kk+1
        else if(nt.eq.2)then
         jj=jj+1
         ifm2=ibc(jj)
         jj=jj+1
         kk2=kk+1
         do l=1,nrow
          spher(l,ito)=bc(kk)*cart(l,ifm)+bc(kk2)*cart(l,ifm2)
         end do
         kk=kk2+1
        else if(nt.eq.3)then
         jj=jj+1
         ifm2=ibc(jj)
         jj=jj+1
         ifm3=ibc(jj)
         jj=jj+1
         kk2=kk+1
         kk3=kk2+1
         do l=1,nrow
          spher(l,ito)=bc(kk)*cart(l,ifm)+bc(kk2)*cart(l,ifm2)
     $         +bc(kk3)*cart(l,ifm3)
         end do
         kk=kk3+1
        else if(nt.eq.4)then
         jj=jj+1
         ifm2=ibc(jj)
         jj=jj+1
         ifm3=ibc(jj)
         jj=jj+1
         ifm4=ibc(jj)
         jj=jj+1
         kk2=kk+1
         kk3=kk2+1
         kk4=kk3+1
         do l=1,nrow
          spher(l,ito)=bc(kk)*cart(l,ifm)+bc(kk2)*cart(l,ifm2)
     $         +bc(kk3)*cart(l,ifm3)+bc(kk4)*cart(l,ifm4)
         end do
         kk=kk4+1
        else if(nt.eq.5)then
         jj=jj+1
         ifm2=ibc(jj)
         jj=jj+1
         ifm3=ibc(jj)
         jj=jj+1
         ifm4=ibc(jj)
         jj=jj+1
         ifm5=ibc(jj)
         jj=jj+1
         kk2=kk+1
         kk3=kk2+1
         kk4=kk3+1
         kk5=kk4+1
         do l=1,nrow
          spher(l,ito)=bc(kk)*cart(l,ifm)+bc(kk2)*cart(l,ifm2)
     $         +bc(kk3)*cart(l,ifm3)+bc(kk4)*cart(l,ifm4)
     $         +bc(kk5)*cart(l,ifm5)
         end do
         kk=kk5+1
        else if(nt.eq.6)then
         jj=jj+1
         ifm2=ibc(jj)
         jj=jj+1
         ifm3=ibc(jj)
         jj=jj+1
         ifm4=ibc(jj)
         jj=jj+1
         ifm5=ibc(jj)
         jj=jj+1
         ifm6=ibc(jj)
         jj=jj+1
         kk2=kk+1
         kk3=kk2+1
         kk4=kk3+1
         kk5=kk4+1
         kk6=kk5+1
         do l=1,nrow
          spher(l,ito)=bc(kk)*cart(l,ifm)+bc(kk2)*cart(l,ifm2)
     $         +bc(kk3)*cart(l,ifm3)+bc(kk4)*cart(l,ifm4)
     $         +bc(kk5)*cart(l,ifm5)+bc(kk6)*cart(l,ifm6)
         end do
         kk=kk6+1
        else
         do l=1,nrow
          spher(l,ito)=bc(kk)*cart(l,ifm)
         end do
         jj=jj+1
         kk=kk+1
         mlow=2
         if(mod(nt,2).eq.0.or.nt.eq.2)then
          mlow=3
          ifm=ibc(jj)
          do l=1,nrow
           spher(l,ito)=spher(l,ito)+bc(kk)*cart(l,ifm)
          end do
          jj=jj+1
          kk=kk+1
         end if
         do m=mlow,nt-1,2
          ifm=ibc(jj)
          ifmn=ibc(jj+1)
          kkp=kk+1
          do l=1,nrow
           spher(l,ito)=spher(l,ito)
     $          +bc(kk)*cart(l,ifm)
     $          +bc(kkp)*cart(l,ifmn)
          end do
          jj=jj+2
          kk=kk+2
         end do
        end if
        jj=jj+1
       end do
      end do
      return
      end
