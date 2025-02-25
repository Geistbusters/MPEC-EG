c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine derh01b(h0,h0d,h0da,idarot,nsymb,nbasdws,idwsdeb,multh,6d16s22
     $     ipuse,bc,ibc)                                                11d14s22
      implicit real*8 (a-h,o-z)                                         5d16s22
      dimension h0(*),h0d(*),h0da(*),idarot(*),nbasdws(*),multh(8,8)    6d16s22
      include "common.store"                                            5d16s22
      ioffd=0                                                           6d16s22
      do isb=1,nsymb                                                    6d16s22
       jsb=multh(isb,ipuse)                                             6d16s22
       ioffd=ioffd+nbasdws(isb)*nbasdws(jsb)                            6d16s22
      end do                                                            6d16s22
      do i=1,ioffd                                                      6d16s22
       h0da(i)=0d0                                                      6d16s22
      end do                                                            6d16s22
      ioff=1                                                            5d16s22
      ioffd=1                                                           6d16s22
      do isb=1,nsymb                                                    5d16s22
       jsb=multh(isb,ipuse)                                             6d16s22
       nn=nbasdws(isb)**2                                               6d16s22
       nnd=nbasdws(isb)*nbasdws(jsb)                                     6d16s22
       if(min(nbasdws(isb),nbasdws(jsb)).gt.0)then                      7d11s22
        if(idwsdeb.gt.10)then                                           5d19s22
         write(6,*)('h0 for sym '),isb,ioff
         call prntm2(h0(ioff),nbasdws(isb),nbasdws(isb),nbasdws(isb))    5d16s22
         write(6,*)('dh from kinematics '),jsb,ioffd
         call prntm2(h0d(ioffd),nbasdws(isb),nbasdws(jsb),nbasdws(isb)) 6d16s22
         write(6,*)('darot ')
         call prntm2(bc(idarot(isb)),nbasdws(jsb),nbasdws(isb),         6d16s22
     $       nbasdws(isb))                                              5d16s22
        end if                                                          5d19s22
        itmp=ibcoff                                                     5d16s22
        ibcoff=itmp+nnd                                                  5d16s22
        call enough('derh01b.  1',bc,ibc)
        call dgemm('n','n',nbasdws(isb),nbasdws(jsb),nbasdws(isb),1d0,  6d16s22
     $       h0(ioff),nbasdws(isb),bc(idarot(isb)),nbasdws(isb),0d0,    5d16s22
     $       bc(itmp),nbasdws(isb),'derh01.6',                          7d11s22
     d'derh01b.  1')
        if(idwsdeb.gt.10)then                                           5d19s22
         sum=0d0
         do i=0,nbasdws(isb)-1
          iad1=ioff+nbasdws(isb)*i
          iad2=idarot(isb)+i
          term=h0(iad1)*bc(iad2)
          sum=sum+term
          if(abs(term).gt.1d-10)write(6,*)i,bc(iad1),bc(iad2),sum
         end do
         write(6,*)('h0*darot ')
         call prntm2(bc(itmp),nbasdws(isb),nbasdws(jsb),nbasdws(isb))    5d16s22
        end if                                                          5d19s22
        if(isb.eq.jsb)then                                              6d16s22
         do i=0,nbasdws(isb)-1                                           5d16s22
          do j=0,nbasdws(isb)-1                                          5d16s22
           ji=j+nbasdws(isb)*i                                           5d16s22
           ij=i+nbasdws(isb)*j                                           5d16s22
           h0da(ioffd+ij)=h0d(ioffd+ij)+bc(itmp+ij)+bc(itmp+ji)         6d16s22
          end do                                                         5d16s22
         end do                                                          5d16s22
        else                                                            6d16s22
         do j=0,nbasdws(jsb)-1                                          6d16s22
          do i=0,nbasdws(isb)-1                                         6d16s22
           ij=i+nbasdws(isb)*j                                          6d16s22
           h0da(ioffd+ij)=h0d(ioffd+ij)+bc(itmp+ij)                      6d16s22
          end do                                                        6d16s22
         end do                                                         6d16s22
         joff=1                                                         6d16s22
         do ii=1,jsb-1                                                  6d16s22
          joff=joff+nbasdws(ii)**2                                      6d16s22
         end do                                                         6d16s22
         call dgemm('n','n',nbasdws(jsb),nbasdws(isb),nbasdws(jsb),1d0, 6d16s22
     $        h0(joff),nbasdws(jsb),bc(idarot(jsb)),nbasdws(jsb),0d0,   6d16s22
     $        bc(itmp),nbasdws(jsb),'derh01.7',                                    6d16s22
     d'derh01b.  2')
         if(idwsdeb.gt.10)then
          write(6,*)('other h0*darot '),joff
          call prntm2(bc(itmp),nbasdws(jsb),nbasdws(isb),nbasdws(jsb))
         end if                                                         6d16s22
         do j=0,nbasdws(jsb)-1                                          6d16s22
          do i=0,nbasdws(isb)-1                                         6d16s22
           ij=i+nbasdws(isb)*j                                          6d16s22
           ji=j+nbasdws(jsb)*i                                          6d16s22
           h0da(ioffd+ij)=h0da(ioffd+ij)+bc(itmp+ji)                    6d16s22
          end do                                                        6d16s22
         end do                                                         6d16s22
        end if                                                          6d16s22
        if(idwsdeb.gt.10)then                                           5d19s22
         write(6,*)('full dh ')
         call prntm2(h0da(ioffd),nbasdws(isb),nbasdws(jsb),nbasdws(isb))  5d16s22
        end if                                                          5d19s22
        ibcoff=itmp                                                     5d16s22
        ioffd=ioffd+nnd                                                 6d16s22
       end if                                                           5d16s22
       ioff=ioff+nn                                                     7d11s22
      end do                                                            5d16s22
      return                                                            5d16s22
      end                                                               5d16s22
