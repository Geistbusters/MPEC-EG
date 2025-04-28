c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dumpdoubc(nff22,nfdat,vd,nsymb,ndoub,mdoub,            8d3s21
     $     nroot,mdon,mdoo,ncsf,ncsf2,ioffset,bc,ibc,nder)              9d18s24
      implicit real*8 (a-h,o-z)                                         8d3s21
      integer*8 ipack8
      dimension vd(*),nfdat(5,4,*),nff22(mdoo+1,2,*),ncsf(*),ncsf2(4,*),8d3s21
     $     ipack4(2)                                                    8d3s21
      equivalence (ipack8,ipack4)
      include "common.store"
c
c     nfdat(1,1,isb) = pointer to first word?
c     nfdat(2,l,isb) = no. non n-2 functions
c     nfdat(3,l,isb) = no. contracted n-2 functions
c     nfdat(4,l,isb) = address of transformation
c     nfdat(5,1,isb) = pointer to bc
c     nff22(ii,1,isb) = no. prim n-2 functions
c     nff22(ii,2,isb) = offset in bc
      nderp=nder+1                                                      9d18s24
      mff2=0                                                            8d3s21
      ifirst=0                                                          8d3s21
      do isb=1,nsymb                                                    8d3s21
       nftot=nfdat(2,1,isb)+nfdat(2,2,isb)+nfdat(2,3,isb)+nfdat(2,4,isb)8d3s21
       if(nftot.gt.0)then                                               8d3s21
        nhere=0                                                         8d3s21
        do ii=mdon+1,mdoo+1                                             8d3s21
         if(nff22(ii,1,isb).gt.0)then                                   8d3s21
          iarg=ii-mdon
          ivcv=nfdat(5,1,isb)+ioffset                                   2d14s22
          jvcv=ivcv+nff22(ii,2,isb)                                     8d3s21
          if(ifirst.eq.0)ifirst=ivcv                                    8d3s21
          do if=1,nff22(ii,1,isb)                                       8d3s21
           ipack8=ibc(jvcv)
           nspace=ibc(jvcv+1)                                           8d3s21
           do l=1,4
            nl=ibc(jvcv+1+l)
            if(nl.gt.0)then
             iad1=jvcv+ibc(jvcv+5+l)
             iad2=iad1+nl
            end if
           end do
           jvcv=jvcv+nspace                                             8d3s21
           mff2=mff2+nspace                                             8d3s21
          end do                                                        8d3s21
         end if                                                         8d3s21
        end do                                                          8d3s21
        do l=1,4                                                        8d3s21
         if(nfdat(2,l,isb).gt.0)then
          mff2=mff2+nderp*nfdat(2,l,isb)*nfdat(3,l,isb)                 9d18s24
         end if
        end do
       end if                                                           8d3s21
      end do                                                            8d3s21
      mff2m=-mff2                                                        8d3s21
      write(2)mff2m
      write(2)ndoub,mdoub                                               8d3s21
      ifirstm=ifirst-1                                                  8d3s21
      write(2)(bc(ifirstm+i),i=1,mff2)                                  8d3s21
      write(2)(((nff22(ii,j,isb),ii=mdon+1,mdoo+1),isb=1,nsymb),j=1,2)   7d8s21
      write(2)(((nfdat(j,l,isb),j=1,5),l=1,4),isb=1,nsymb)              8d3s21
      write(2)(ncsf(iarg),iarg=1,mdoo+1-mdon)                           7d8s21
      write(2)((ncsf2(l,iarg),l=1,4),iarg=1,mdoo+1-mdon)                8d3s21
      write(2)(vd(i),i=1,ndoub*nroot)                                   8d3s21
      return
      end
