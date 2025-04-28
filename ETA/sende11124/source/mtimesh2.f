c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine mtimesh2(xin,xout,nsymb,nroot,nadet,nbdet,nsbeta,      4d10s23
     $     nherec,nherect,icode,ilc,ihc,ilct,ihct,hdig,ih0e,aorb,borb,  4d7s23
     $     nalpha,nbeta,i4oa,m12sym,norb,idata1,idata2,idatb1,idatb2,   4d7s23
     $     m1c,idatac,idatbc,mst,mnz,nz,lbail,bc,ibc,ivin,eig)                   4d7s23
      implicit real*8 (a-h,o-z)
      dimension nsbeta(*),ig(8),igt(8),nadet(*),nbdet(*),nherec(*),     4d7s23
     $     nherect(*),iva(8),iv(8),ivt(8),xin(*),xout(*),               4d7s23
     $     ivin(*),eig(*),iga(8),ilc(*),ihc(*),ilct(*),ihct(*)          4d7s23
      include "common.store"
      ibcoffo=ibcoff                                                    4d7s23
      igoal=8
      ncona=0d0                                                         4d7s23
      do isb=1,nsymb                                                    5d27s22
       jsb=nsbeta(isb)                                                  5d27s22
       ig(isb)=ibcoff                                                   5d27s22
       igt(isb)=ig(isb)+nroot*nbdet(jsb)*nherec(isb)                    5d27s22
       iva(isb)=igt(isb)+nroot*nadet(jsb)*nherect(isb)                    5d27s22
       iv(isb)=iva(isb)+nroot*nadet(isb)*nbdet(jsb)                     4d7s23
       ivt(isb)=iv(isb)+nroot*nbdet(jsb)*nherec(isb)                    4d7s23
       ibcoff=ivt(isb)+nroot*nadet(jsb)*nherect(isb)                    4d7s23
       ncona=ncona+nadet(isb)*nbdet(jsb)                                4d7s23
      end do                                                            5d27s22
      idot=ibcoff                                                       4d7s23
      ibcoff=idot+nroot                                                 4d7s23
      call enough('mtimesh.  1',bc,ibc)
      nrootm=nroot-1                                                    5d27s22
      do irt=0,nrootm                                                   4d7s23
       bc(idot+irt)=0d0                                                 4d7s23
      end do                                                            4d7s23
      ii=1                                                              4d7s23
      do isb=1,nsymb                                                    4d7s23
       jsb=nsbeta(isb)                                                  4d7s23
       ncol=nadet(isb)*nbdet(jsb)                                       4d7s23
       ncolr=ncol*nroot                                                 4d7s23
       iad1=iva(isb)                                                    4d7s23
       do i=0,ncolr-1                                                   4d7s23
        bc(iad1+i)=xin(ii+i)                                            4d7s23
       end do                                                           4d7s23
       jvin=ivin(isb)
       do i=0,ncol-1                                                    4d7s23
        do irt=0,nrootm                                                 4d7s23
         bc(idot+irt)=bc(idot+irt)+xin(ii+irt)*bc(jvin+irt)             4d7s23
        end do                                                          4d7s23
        ii=ii+nroot                                                     4d7s23
        jvin=jvin+nroot                                                 4d7s23
       end do                                                           4d7s23
      end do                                                            4d7s23
      do isb=1,nsymb                                                    5d27s22
       jsb=nsbeta(isb)                                                  5d27s22
       do icol=0,nadet(jsb)-1                                           5d27s22
        do irow=0,nherect(isb)-1                                        5d27s22
         jrow=irow+ilct(isb)-1                                          5d27s22
         iad1=ivt(isb)+nroot*(irow+nherect(isb)*icol)                   4d7s23
         iad2=iva(jsb)+nroot*(icol+nadet(jsb)*jrow)                     4d7s23
         do irt=0,nrootm                                                 5d27s22
          bc(iad1+irt)=bc(iad2+irt)                                       5d27s22
         end do                                                         5d27s22
        end do                                                          5d27s22
       end do                                                           5d27s22
       do icol=0,nbdet(jsb)-1                                           5d27s22
        do irow=0,nherec(isb)-1                                         5d27s22
         jrow=irow+ilc(isb)-1                                           5d27s22
         iad1=iv(isb)+nroot*(irow+nherec(isb)*icol)                     4d7s23
         iad2=iva(isb)+nroot*(jrow+nadet(isb)*icol)                     4d7s23
         do irt=0,nrootm                                                 5d27s22
          bc(iad1+irt)=bc(iad2+irt)                                       5d27s22
         end do                                                         5d27s22
        end do                                                          5d27s22
       end do                                                           5d27s22
      end do                                                            5d27s22
      call hc1c(iv,ivt,iva,ig,igt,ilc,ihc,ilct,ihct,nroot,              4d7s23
     $       hdig,ih0e,nherec,nherect,aorb,borb,                        4d7s23
     $       nalpha,nbeta,i4oa,m12sym,norb,nsymb,idata1,idata2,idatb1,  5d27s22
     $       idatb2,nsbeta,m1c,idatac,idatbc,mst,mnz,nz,lbail,bc,ibc)   12d19s22
      do i=1,ncona*nroot                                                4d7s23
       xout(i)=0d0                                                      4d7s23
      end do                                                            4d6s23
      ii=1                                                              4d7s23
      do isb=1,nsymb                                                    4d7s23
       jsb=nsbeta(isb)                                                  4d7s23
       iga(isb)=ii                                                      4d7s23
       do icol=0,nbdet(jsb)-1                                           3d14s17
        do irow=0,nherec(isb)-1                                         3d14s17
         iad1=ig(isb)+nroot*(irow+nherec(isb)*icol)                     5d27s22
         iad2=ii+nroot*(irow+ilc(isb)-1+nadet(isb)*icol)                4d4s23
         do i3=0,nrootm                                                 5d27s22
          xout(iad2+i3)=bc(iad1+i3)                                       3d14s17
         end do                                                         3d14s17
        end do                                                          3d14s17
       end do                                                           3d14s17
       ii=ii+nadet(isb)*nroot*nbdet(jsb)                                4d4s23
      end do                                                            3d12s17
      do jsb=1,nsymb                                                    3d12s17
       isb=nsbeta(jsb)                                                  3d12s17
       do icol=0,nadet(isb)-1                                           3d14s17
        do irow=0,nherect(jsb)-1                                        3d14s17
         iad1=igt(jsb)+nroot*(irow+nherect(jsb)*icol)                   5d27s22
         iad2=iga(isb)+nroot*(icol+nadet(isb)*(irow+ilct(jsb)-1))       5d27s22
         do i3=0,nrootm                                                 5d27s22
          xout(iad2+i3)=xout(iad2+i3)+bc(iad1+i3)                           3d14s17
         end do                                                         3d14s17
        end do                                                          3d14s17
       end do                                                           3d14s17
      end do
      call dws_gsumf(xout,ncona*nroot)                                      4d4s23
      do isb=1,nsymb                                                    5d27s22
       jsb=nsbeta(isb)                                                  5d27s22
       ncol=nadet(isb)*nbdet(jsb)                                       5d27s22
       ncolm=ncol-1                                                     5d27s22
       iad1=iga(isb)                                                    5d27s22
       iad3=iva(isb)                                                    5d27s22
       if(icode.eq.1)then                                               6d15s22
        iad2=ivin(isb)                                                    5d27s22
        do i=0,ncolm                                                    5d27s22
         do irt=0,nrootm                                                 5d27s22
          irtp=irt+1                                                    4d7s23
          xout(iad1+irt)=xout(iad1+irt)                                 4d7s23
     $          +eig(irtp)*(bc(idot+irt)*bc(iad2+irt)-bc(iad3+irt))     4d7s23
         end do                                                         5d27s22
         iad1=iad1+nroot                                                5d27s22
         iad2=iad2+nroot                                                5d27s22
         iad3=iad3+nroot                                                5d27s22
        end do                                                          5d27s22
       else                                                             6d15s22
        do i=0,ncolm                                                    5d27s22
         do irt=0,nrootm                                                 5d27s22
          irtp=irt+1                                                    4d7s23
          xout(iad1+irt)=xout(iad1+irt)-eig(irtp)*bc(iad3+irt)             4d4s23
         end do                                                         5d27s22
         iad1=iad1+nroot                                                5d27s22
         iad3=iad3+nroot                                                5d27s22
        end do                                                          5d27s22
       end if                                                           6d15s22
      end do                                                            5d27s22
      ibcoff=ibcoffo                                                    4d7s23
      return                                                            4d7s23
      end                                                               4d7s23
