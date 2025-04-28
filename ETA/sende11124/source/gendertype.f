c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine gendertype(ipsym,multh,noc,isblkder,isblkxder,nsblkder,8d3s16
     $     nsblkxder,isblkkder,nsblkkder,ipxder)                        11d28s22
      implicit real*8 (a-h,o-z)
c
c     generate symmetry type indicies for (oooo)' and (ooox)'
c     when ' has symmetry ipsym
c     and K.
c
      include "common.hf"
      dimension multh(8,8),noc(8),isblkder(4,idbk),isblkxder(4,idbk),   8d3s16
     $     isblkkder(4,idbk),ipxder(4,8,8,8)                            6d27s22
      nsblkder=0
      do isd=1,nsymb
       isdd=multh(isd,ipsym)
       do isc=1,nsymb
        iscd=multh(isdd,isc)
        do isb=1,nsymb
         isad=multh(iscd,isb)                                           4d18s16
         isa=multh(isad,ipsym)                                          4d18s16
         do i=1,nsblkder
          if(noc(isad)*noc(isb)*noc(isc)*noc(isd).eq.0)go to 1          7d7s22
          if(isad.eq.isblkder(1,i).and.isb.eq.isblkder(2,i))then        4d18s16
           if((isc.eq.isblkder(3,i).and.isd.eq.isblkder(4,i)).or.
     $         (isd.eq.isblkder(3,i).and.isc.eq.isblkder(4,i)))then
            go to 1
           end if
          else if(isb.eq.isblkder(1,i).and.isad.eq.isblkder(2,i))then   4d18s16
           if((isc.eq.isblkder(3,i).and.isd.eq.isblkder(4,i)).or.
     $          (isd.eq.isblkder(3,i).and.isc.eq.isblkder(4,i)))then
            go to 1
           end if
          else if(isc.eq.isblkder(1,i).and.isd.eq.isblkder(2,i))then
           if((isad.eq.isblkder(3,i).and.isb.eq.isblkder(4,i)).or.      4d18s16
     $          (isb.eq.isblkder(3,i).and.isad.eq.isblkder(4,i)))then
           end if
          else if(isd.eq.isblkder(1,i).and.isc.eq.isblkder(2,i))then
           if((isad.eq.isblkder(3,i).and.isb.eq.isblkder(4,i)).or.      4d18s16
     $          (isb.eq.isblkder(3,i).and.isad.eq.isblkder(4,i)))then
            go to 1
           end if
          end if
c     final ders will show full 8 way symmetry, but when computing
c     ders, it is convenient to make the 4th index the sole der,
c     so we only get 2 way symmetry. this is the same as for onex.
         end do
         nsblkder=nsblkder+1
         if(nsblkder.gt.idbk)then
          write(6,*)('nsblkder exceeds idbk '),idbk
          call dws_sync
          call dws_finalize
          stop
         end if
         isblkder(1,nsblkder)=isad                                      4d18s16
         isblkder(2,nsblkder)=isb
         isblkder(3,nsblkder)=isc
         isblkder(4,nsblkder)=isd
         ipxder(1,isad,isb,isc)=nsblkder                                6d13s22
         ipxder(2,isad,isb,isc)=0                                       6d13s22
         ipxder(3,isad,isb,isc)=0                                       6d13s22
         ipxder(4,isad,isb,isc)=0                                       6d13s22
         ipxder(1,isc,isd,isad)=nsblkder                                6d13s22
         ipxder(2,isc,isd,isad)=0                                       6d13s22
         ipxder(3,isc,isd,isad)=0                                       6d13s22
         ipxder(4,isc,isd,isad)=1                                       6d13s22
         if(isad.ne.isb)then                                            6d13s22
          ipxder(1,isb,isad,isc)=nsblkder                               6d13s22
          ipxder(2,isb,isad,isc)=1                                      6d13s22
          ipxder(3,isb,isad,isc)=0                                      6d13s22
          ipxder(4,isb,isad,isc)=0                                      6d13s22
          ipxder(1,isc,isd,isb)=nsblkder                                6d27s22
          ipxder(2,isc,isd,isb)=0                                       6d29s22
          ipxder(3,isc,isd,isb)=1                                       6d29s22
          ipxder(4,isc,isd,isb)=1                                       6d27s22
         end if                                                         6d13s22
         if(isc.ne.isd)then                                             6d13s22
          ipxder(1,isad,isb,isd)=nsblkder                                6d13s22
          ipxder(2,isad,isb,isd)=0                                       6d13s22
          ipxder(3,isad,isb,isd)=1                                       6d13s22
          ipxder(4,isad,isb,isd)=0                                       6d13s22
          ipxder(1,isd,isc,isad)=nsblkder                               6d27s22
          ipxder(2,isd,isc,isad)=1                                      6d29s22
          ipxder(3,isd,isc,isad)=0                                      6d29s22
          ipxder(4,isd,isc,isad)=1                                      6d27s22
          if(isad.ne.isb)then                                           6d13s22
           ipxder(1,isb,isad,isd)=nsblkder                              6d13s22
           ipxder(2,isb,isad,isd)=1                                     6d13s22
           ipxder(3,isb,isad,isd)=1                                     6d13s22
           ipxder(4,isb,isad,isd)=0                                     6d13s22
           ipxder(1,isd,isc,isb)=nsblkder                               6d27s22
           ipxder(2,isd,isc,isb)=1                                      6d27s22
           ipxder(3,isd,isc,isb)=1                                      6d27s22
           ipxder(4,isd,isc,isb)=1                                      6d27s22
          end if                                                        6d13s22
         end if                                                         6d13s22
    2    format(i3,2x,4i1,i3)
    1    continue
        end do
       end do
      end do
      nsblkxder=0
      do isd=1,nsymb
       isdd=multh(isd,ipsym)
       do isc=1,nsymb
        iscd=multh(isdd,isc)
        do isb=1,nsymb
         isad=multh(iscd,isb)                                           4d18s16
         if(min(noc(isad),noc(isb),noc(isc),nvirt(isd)).eq.0)go to 3    7d1s22
c     we need to distinguish between 3rd and 4th index for we
c     are only differentiating the 4th index.
         do i=1,nsblkxder
          if(isc.eq.isblkxder(3,i).and.isd.eq.isblkxder(4,i))then
           if((isad.eq.isblkxder(1,i).and.isb.eq.isblkxder(2,i)).or.    4d18s16
     $        (isb.eq.isblkxder(1,i).and.isad.eq.isblkxder(2,i)))go to 34d18s16
          end if
         end do
         nsblkxder=nsblkxder+1
         if(nsblkxder.gt.idbk)then
          write(6,*)('nsblkxder exceeds idbk '),idbk
          call dws_sync
          call dws_finalize
          stop
         end if
         isblkxder(1,nsblkxder)=isad                                    4d18s16
         isblkxder(2,nsblkxder)=isb
         isblkxder(3,nsblkxder)=isc
         isblkxder(4,nsblkxder)=isd
    3    continue
        end do
       end do
      end do
      if(nsblkkder.lt.0)return                                          8d4s16
      do isd=1,nsymb
       isdd=multh(isd,ipsym)
       do isc=1,nsymb
        iscd=multh(isdd,isc)
        do isb=1,nsymb
         isad=multh(iscd,isb)                                           4d18s16
         isbd=multh(isb,ipsym)                                          8d24s16
         isa=multh(isad,ipsym)                                          4d18s16
         if(noc(isa)*noc(isad)*noc(isb)*nbasdws(isc)*nbasdws(isd).eq.0  8d24s16
     $      .and.noc(isad)*noc(isb)*noc(isbd)*nbasdws(isc)*nbasdws(isd) 8d24s16
     $        .eq.0.and.                                                8d24s16
     $      noc(isad)*noc(isb)*nbasdws(iscd)*nbasdws(isd).eq.0.and.      8d24s16
     $      noc(isad)*noc(isb)*nbasdws(isc)*nbasdws(isdd).eq.0)go to 4  8d24s16
         do i=1,nsblkkder
          if(isc.eq.isblkkder(3,i).and.isd.eq.isblkkder(4,i))then
           if((isad.eq.isblkkder(1,i).and.isb.eq.isblkkder(2,i)).or.    4d18s16
     $        (isb.eq.isblkkder(1,i).and.isad.eq.isblkkder(2,i)))go to 44d18s16
          end if
         end do
         nsblkkder=nsblkkder+1
         if(nsblkkder.gt.idbk)then
          write(6,*)('nsblkkder exceeds idbk '),idbk
          call dws_sync
          call dws_finalize
          stop
         end if
         isblkkder(1,nsblkkder)=isad                                    4d18s16
         isblkkder(2,nsblkkder)=isb
         isblkkder(3,nsblkkder)=isc
         isblkkder(4,nsblkkder)=isd
    4    continue
        end do
       end do
      end do
      return
      end
