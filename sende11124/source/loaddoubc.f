c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine loaddoubc(iff22,nff22,nfdat,ncsf,ncsf2,mdon,mdoo,      8d3s21
     $     nsymb,mff2a,vd,nall,isymmrci,nroot,multh,nvirt,inorm)        8d10s22
      implicit real*8 (a-h,o-z)                                         8d3s21
c
c     inorm=
c     1 read from unit 1 and store vectors. This is for prop.           8d18s22
c     2 read from unit 2 and store vectors. This is for mrci restart.   8d18s22
c                                                                       8d18s22
      integer*8 iff22(mff2a),ipack8                                     8d18s22
      dimension nff22(mdoo+1,2,*),nfdat(5,4,*),ncsf(*),ncsf2(4,*),      8d3s21
     $     vd(nall),multh(8,*),nvirt(*),ipack4(2)                       8d18s22
      equivalence (ipack8,ipack4)                                       8d18s22
      if(inorm.eq.2)then                                                8d18s22
       iunit=2                                                          8d18s22
      else                                                              8d18s22
       iunit=1                                                          8d18s22
      end if                                                            8d18s22
      read(iunit)iff22                                                  8d18s22
      ipack8=iff22(1)
      read(iunit)                                                       8d18s22
     $     (((nff22(ii,j,isb),ii=mdon+1,mdoo+1),isb=1,nsymb),j=1,2)     8d18s22
      read(iunit)(((nfdat(j,l,isb),j=1,5),l=1,4),isb=1,nsymb)           8d18s22
      read(iunit)(ncsf(iarg),iarg=1,mdoo+1-mdon)                        8d18s22
      read(iunit)((ncsf2(l,iarg),l=1,4),iarg=1,mdoo+1-mdon)             8d18s22
      read(iunit)vd                                                     8d18s22
      ifirst=0                                                          8d3s21
      do isb=1,nsymb                                                    8d3s21
       nftot=nfdat(2,1,isb)+nfdat(2,2,isb)+nfdat(2,3,isb)               8d18s22
     $       +nfdat(2,4,isb)                                            8d18s22
       if(nftot.gt.0.and.ifirst.eq.0)ifirst=nfdat(5,1,isb)              8d3s21
      end do                                                            8d3s21
      nshift=ifirst-1                                                   8d3s21
      do isb=1,nsymb                                                    8d3s21
       nfdat(5,1,isb)=nfdat(5,1,isb)-nshift                             8d3s21
       do l=1,4                                                         8d3s21
        if(nfdat(2,l,isb).gt.0)then                                     8d3s21
         nfdat(4,l,isb)=nfdat(4,l,isb)-nshift                           8d3s21
        end if                                                          8d3s21
       end do                                                           8d3s21
      end do                                                            8d3s21
      return                                                            8d3s21
      end                                                               8d3s21
