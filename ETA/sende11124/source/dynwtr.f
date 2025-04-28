c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dynwtr(nroot,eig,wgt,eref,deref,dynw,ewgt,wsum,dwsum,  7d20s22
     $     shift,dshift,icode)                                          7d20s22
      implicit real*8 (a-h,o-z)                                         7d20s22
c
c     icode=0: compute dynamic weights.
c     icode ne 0: compute dynamics weights and energy derivatives
c
      dimension eig(nroot,2),wgt(nroot,2),dynw(3),ewgt(*)               7d20s22
      if(icode.eq.0)then                                                7d20s22
       do ir=1,nroot                                                     7d20s22
        if(dynw(2).eq.0d0)then                                           7d20s22
         ex=exp(-dynw(1)*(eig(ir,1)-eref))                               7d20s22
         exi=1d0/ex                                                      7d20s22
         wwww=(2d0/(ex+exi))**2                                          7d20s22
        else                                                             7d20s22
         wde=dynw(1)*(eig(ir,1)+shift-dynw(2))                           7d20s22
         wwww=4d0/(dynw(3)+exp(2d0*wde))                                 7d20s22
        end if                                                           7d20s22
        wgt(ir,1)=ewgt(ir)*max(1d-7,wwww)                                7d20s22
        wsum=wsum+wgt(ir,1)                                              7d20s22
       end do                                                            7d20s22
      else                                                              7d20s22
       ffact=1d0
       do ir=1,nroot                                                     7d20s22
        if(dynw(2).eq.0d0)then                                           7d20s22
         ex=exp(-dynw(1)*(eig(ir,1)-eref))                               7d20s22
         dex=-ex                                                        7d20s22
         exi=1d0/ex                                                      7d20s22
         dexi=-dex*exi*exi                                              7d20s22
         esum=ex+exi                                                    7d20s22
         esumi=1d0/esum                                                 7d20s22
         esumi2=esumi*esumi                                             7d20s22
         desum=dex+dexi                                                 7d20s22
         wwww=4d0*esumi2                                                7d20s22
         dwww0=-2d0*wwww*esumi*desum                                    7d20s22
         dwww=dynw(1)*dwww0*(eig(ir,2)-deref)                           7d27s22
        else                                                             7d20s22
         wde=dynw(1)*(eig(ir,1)+shift-dynw(2))                           7d20s22
         dwde=dynw(1)                                                   7d20s22
         ex=exp(2d0*wde)                                                7d20s22
         dex=2d0*dwde*ex                                                7d20s22
         bot=1d0/(dynw(3)+ex)                                           7d20s22
         wwww=4d0*bot                                                   7d20s22
         dwww=-wwww*bot*dex*eig(ir,2)                                   7d20s22
        end if                                                           7d20s22
        wgt(ir,1)=ewgt(ir)*max(1d-7,wwww)                                7d20s22
        wgt(ir,2)=ewgt(ir)*dwww*ffact                                         7d20s22
        wsum=wsum+wgt(ir,1)                                              7d20s22
        dwsum=dwsum+wgt(ir,2)                                           7d20s22
       end do                                                            7d20s22
      end if                                                            7d20s22
      return                                                            7d20s22
      end                                                               7d20s22
