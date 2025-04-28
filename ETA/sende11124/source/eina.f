c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      real*8 function eina(tfi,gi,we,ci,icode,aaa)
      implicit real*8 (a-h,o-z)
c
c     compute einstein a
c     icode=1 e1
c     icode=2 m1
c     icode=3 e2
c     ci is 1/c
c     aaa is conversion from 1/au time to 1/sec.
c     gi is 1/g initial
c
      eina=tfi**2
      iq=max(1,icode-1)
      ff=2d0/dfloat(iq*2+1)
      eina=eina*ff*gi
      wec=we*ci
      if(icode.eq.1)then
       eina=eina*2d0*(wec**3)
      else if(icode.eq.2)then
       eina=eina*(wec**3)*ci*ci/8d0
      else
       eina=eina*(wec**5)/4d0
      end if
      eina=eina*aaa
      if(eina.ne.eina)then                                              11d21s22
      write(6,*)('einstein A is nan! '),eina
      write(6,*)('tfi,icode,we,ci,aaa '),tfi,icode,we,ci,aaa
      stop 'eina'
      end if                                                            11d21s22
      return
      end
