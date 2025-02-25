c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine sandtd(xin,nprim,ncont,ncol,trans,bc,ibc)              11d10s22
      implicit real*8 (a-h,o-z)
c
c     square and transform to diagonals
c     overwrite xin with result
c
      dimension xin(*),trans(nprim,ncont)                               11d23s20
      include "common.store"                                            11d23s20
      itmp1=ibcoff                                                      11d23s20
      itmp2=itmp1+nprim*nprim                                           11d23s20
      ibcoff=itmp2+nprim*ncont                                          11d23s20
      call enough('sandtd.  1',bc,ibc)
      ntri=(nprim*(nprim+1))/2                                          11d23s20
      ntrim=ntri-1                                                      11d23s20
      ioff=1                                                            11d23s20
      koff=1                                                            11d23s20
      do icol=1,ncol                                                    11d23s20
       do i=0,ntrim                                                     11d23s20
        bc(itmp1+i)=xin(ioff+i)                                         11d23s20
       end do                                                           11d23s20
       call square(bc(itmp1),nprim)                                     11d23s20
       call dgemm('n','n',nprim,ncont,nprim,1d0,bc(itmp1),nprim,        11d23s20
     $      trans,nprim,0d0,bc(itmp2),nprim,                            11d23s20
     d' sandtd.  1')
       jtmp2=itmp2-1                                                    11d23s20
       ioff=ioff+ntri                                                   11d23s20
       do i=1,ncont                                                     11d23s20
        sum=0d0                                                         11d23s20
        do j=1,nprim                                                    11d23s20
         sum=sum+trans(j,i)*bc(jtmp2+j)                                 11d23s20
        end do                                                          11d23s20
        xin(koff)=sum                                                   11d23s20
        koff=koff+1                                                     11d23s20
        jtmp2=jtmp2+nprim                                               11d23s20
       end do                                                           11d23s20
      end do                                                            11d23s20
      ibcoff=itmp1                                                      11d23s20
c
c     transpose if ncol gt 1
c
      if(ncol.gt.1)then                                                 11d23s20
       itmp=ibcoff                                                      11d23s20
       ibcoff=itmp+ncont*ncol                                           11d23s20
       call enough('sandtd.  2',bc,ibc)
       do i=0,ncol-1                                                    11d23s20
        do j=0,ncont-1                                                  11d23s20
         ji=1+j+ncont*i                                                 11d23s20
         ij=itmp+i+ncol*j                                               11d23s20
         bc(ij)=xin(ji)                                                 11d23s20
        end do                                                          11d23s20
       end do                                                           11d23s20
       jtmp=itmp-1                                                      11d23s20
       do i=1,ncol*ncont                                                11d23s20
        xin(i)=bc(jtmp+i)                                               11d23s20
       end do                                                           11d23s20
       ibcoff=itmp                                                      11d23s20
      end if                                                            11d23s20
      return                                                            11d23s20
      end                                                               11d23s20
