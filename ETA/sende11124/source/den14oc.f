c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine den14oc(jsa,jsb,jsc,jsd,jga,jgb,jgc,jgd,irefo,itmp3,   9d29s23
     $     nli,iad1,nfdat,iden14o,l,bc,ibc,fmul,nlj,jad1)               9d29s23
      implicit real*8 (a-h,o-z)                                         9d29s23
      include "common.store"                                            9d29s23
      dimension irefo(*),nli(*),nfdat(5,*),iden14o(4,*),nlj(*)          9d29s23
      i2eu=ifind2(jsa,jsb,jsc,jsd,icase)                                7d7s19
      nrow=irefo(jsa)*irefo(jsb)                                        7d28s22
      ncol=irefo(jsc)*irefo(jsd)                                        7d28s22
      if(jsa.eq.jsb)then                                                7d7s19
       nrow=(irefo(jsa)*(irefo(jsa)+1))/2                               7d7s19
       ix=max(jga,jgb)                                                  7d7s19
       in=min(jga,jgb)                                                  7d7s19
       irow=((ix*(ix+1))/2)+in                                          7d7s19
       ix=max(jgc,jgd)                                                  7d7s19
       in=min(jgc,jgd)                                                  7d7s19
       icol=((ix*(ix+1))/2)+in                                          7d7s19
       igot=irow+nrow*icol                                              9d29s23
      else if(icase.eq.1)then
       igot=jga+irefo(jsa)*(jgb+irefo(jsb)*(jgc+irefo(jsc)*jgd))        9d29s23
      else if(icase.eq.2)then
       igot=jgb+irefo(jsb)*(jga+irefo(jsa)*(jgd+irefo(jsd)*jgc))        9d29s23
      else if(icase.eq.3)then
       igot=jga+irefo(jsa)*(jgb+irefo(jsb)*(jgd+irefo(jsd)*jgc))        9d29s23
      else if(icase.eq.4)then
       igot=jgb+irefo(jsb)*(jga+irefo(jsa)*(jgc+irefo(jsc)*jgd))        9d29s23
      end if                                                            7d7s19
      jtmp3=itmp3                                                        2d24s21
      do j=0,nlj(l)-1                                                     2d23s21
       jj=ibc(jad1+j)-1                                                    1d8s21
       jden=iden14o(l,i2eu)+nfdat(2,l)*(jj                              9d29s23
     $      +nfdat(2,l)*igot)-1                                         9d29s23
       do k=0,nli(l)-1                                                   2d24s21
        kk=ibc(iad1+k)                                                   2d24s21
        bc(jden+kk)=bc(jden+kk)+bc(jtmp3+k)*fmul                        9d29s23
       end do                                                            2d24s21
       jtmp3=jtmp3+nli(l)                                                2d24s21
      end do                                                             2d24s21
      return                                                            9d29s23
      end                                                               9d29s23
