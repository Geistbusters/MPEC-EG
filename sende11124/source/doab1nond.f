c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine doab1nond(idat1a,n1a,idat1b,n1b,g,nhere,ilc,ihc,cvec,  7d25s22
     $     nadet,i2e,ism,irelo,iacto,iu1,iu2,ju1,ju2,ipxder,bc,ibc,wgt, 3d15s23
     $     nroot)                                                       3d15s23
      implicit real*8 (a-h,o-z)
      integer*8 iarg1,iarg2,ism,irelo                                   12d13s06
      integer*2 idat1b(4,*),idat1a(4,*),idwstmp2                        6d15s22
      integer*1 idwstmp1(2)                                             6d15s22
      equivalence (idwstmp2,idwstmp1)                                   5d30s06
      include "common.store"                                            10d6s06
      dimension cvec(nroot,nadet,*),g(nroot,nhere,*),wgt(*),            3d15s23
     $     i2e(*),ism(*),irelo(*),iacto(*),ipxder(4,8,8,8)              6d28s22
      dimension phs(2)
      data phs/1d0,-1d0/
      do ja=1,n1a
       jjb=idat1a(iu1,ja)                                               7d25s22
       jjk=idat1a(iu2,ja)                                               7d25s22
       idwstmp2=idat1a(3,ja)                                              5d30s06
       jsbb=ism(idwstmp1(iu1))                                            8d29s06
       jsbk=ism(idwstmp1(iu2))                                            8d29s06
       jnb=irelo(idwstmp1(iu1))-1                                         6d13s22
       jnk=irelo(idwstmp1(iu2))-1                                         6d13s22
       pa=phs(idat1a(4,ja))                                             6d15s22
       if(jjb.ge.ilc.and.jjb.le.ihc)then                                6d15s22
        jjbr=jjb+1-ilc                                                  6d15s22
        do ib=1,n1b
         iib=idat1b(ju1,ib)                                             7d25s22
         iik=idat1b(ju2,ib)                                             7d25s22
         idwstmp2=idat1b(3,ib)                                              5d30s06
         isbb=ism(idwstmp1(ju1))                                            8d29s06
         isbk=ism(idwstmp1(ju2))                                            8d29s06
         inb=irelo(idwstmp1(ju1))-1                                         6d13s22
         ink=irelo(idwstmp1(ju2))-1                                         6d13s22
         i2eu=ipxder(1,jsbb,jsbk,isbb)                                  6d15s22
         ixa=jnb+iacto(jsbb)*jnk                                        6d15s22
         ixb=ixa+ipxder(2,jsbb,jsbk,isbb)*(jnk+iacto(jsbk)*jnb-ixa)     6d15s22
         ixc=inb+iacto(isbb)*ink                                         6d15s22
         ixd=ixc+ipxder(3,jsbb,jsbk,isbb)*(ink+iacto(isbk)*inb-ixc)     6d15s22
         ix12=ixb+iacto(jsbb)*iacto(jsbk)*ixd                           6d28s22
         ix21=ixd+iacto(isbb)*iacto(isbk)*ixb-ix12                      6d28s22
         iadc=i2e(i2eu)+ix12+ipxder(4,jsbb,jsbk,isbb)*ix21              6d28s22
         vv=0d0                                                         3d15s23
         do ir=1,nroot                                                  3d15s23
          vv=vv+g(ir,jjbr,iib)*cvec(ir,jjk,iik)*wgt(ir)                 3d15s23
         end do                                                         3d15s23
         bc(iadc)=bc(iadc)
     $        +vv*pa*phs(idat1b(4,ib))                                  3d15s23
        end do                                                          6d15s22
       end if                                                           6d15s22
      end do
      return
      end
