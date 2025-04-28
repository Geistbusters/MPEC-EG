c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine corelz2(ixlzz,shift,noc4,ixlzze,islz,multh,nlzz,bc,ibc)11d14s22
      implicit real*8 (a-h,o-z)                                         12d22s19
c
c     take care of double occupied orbitals for lz2 operator.
c
      dimension ixlzz(8,nlzz),noc4(8),ixlzze(8,*),multh(8,8),islz(*)    12d31s19
      include "common.store"
      include "common.hf"
      include "common.cas"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      shift=0d0                                                         12d22s19
      if(nlzz.eq.2)then                                                 12d31s19
       ione=2                                                           12d31s19
      else                                                              12d31s19
       ione=4                                                           12d31s19
      end if                                                            12d31s19
      do isb=1,nsymb                                                    12d22s19
       if(noc4(isb).gt.0)then                                           12d31s19
        if(nlzz.eq.2)then                                                12d31s19
         ione=2                                                          12d31s19
        else                                                             12d31s19
         ione=4                                                          12d31s19
        end if                                                           12d31s19
        ixlzze(isb,ione)=ibcoff                                          12d31s19
        ibcoff=ixlzze(isb,ione)+iacto(isb)*iacto(isb)                    12d31s19
        call enough('corelz2.  1',bc,ibc)
        do i=ixlzze(isb,ione),ibcoff-1                                   12d31s19
         bc(i)=0d0                                                       12d31s19
        end do                                                           12d31s19
        do ipass=1,nlzz/2                                                12d31s19
         if(nlzz.eq.2)then                                               12d31s19
          iuse=2                                                         12d31s19
         else                                                            12d31s19
          iuse=3+ipass                                                   12d31s19
         end if                                                          12d31s19
         do i=0,idoub(isb)-1                                              12d22s19
          iad=ixlzz(isb,iuse)+i*(noc4(isb)+1)                            12d31s19
          shift=shift+2d0*bc(iad)                                         12d22s19
         end do                                                           12d22s19
         do i=0,iacto(isb)-1                                                    12d22s19
          ip=i+idoub(isb)                                                 12d22s19
          iad1=ixlzze(isb,ione)+iacto(isb)*i                             12d31s19
          iad2=ixlzz(isb,iuse)+idoub(isb)+noc4(isb)*ip                   12d31s19
          do j=0,iacto(isb)-1                                                   12d22s19
           bc(iad1+j)=bc(iad1+j)+bc(iad2+j)                              12d31s19
          end do                                                          12d22s19
         end do                                                           12d22s19
        end do                                                           12d31s19
       end if
      end do                                                            12d22s19
c
c     2-electron part.
c
      do ipass=1,nlzz/2                                                 12d31s19
       do idws=1,nsdlk                                                   12d22s19
        n1=isblk(1,idws)                                                 12d22s19
        n2=isblk(2,idws)                                                 12d22s19
        n3=isblk(3,idws)                                                 12d22s19
        n4=isblk(4,idws)                                                 12d22s19
        if(multh(n1,n2).eq.islz(ipass))then                             12d31s19
         if(n1.eq.n2.and.n3.eq.n4)then                                    12d22s19
          do i34=0,idoub(n4)-1                                            12d22s19
           iad34=ixlzz(n3,ipass)+i34*(noc4(n3)+1)                       12d31s19
           do i12=0,idoub(n1)-1
            iad12=ixlzz(n1,ipass)+i12*(noc4(n1)+1)                      12d31s19
            term=-2d0*bc(iad34)*bc(iad12)                                12d22s19
            shift=shift+2d0*term                                         12d22s19
           end do                                                        12d22s19
          end do                                                         12d22s19
         end if                                                          12d22s19
        end if                                                           12d22s19
        if(multh(n1,n3).eq.islz(ipass))then                             12d31s19
         if(n1.eq.n3.and.n2.eq.n4)then                                   12d22s19
          fact=1d0                                                       12d22s19
          if(n1.ne.n2)fact=2d0                                           12d22s19
          do i24=0,idoub(n4)-1                                           12d22s19
           iad23=ixlzz(n2,ipass)+i24*(noc4(n2)+1)                       12d31s19
           do i13=0,idoub(n3)-1                                          12d22s19
            iad13=ixlzz(n1,ipass)+i13*(noc4(n1)+1)                      12d31s19
            term=-2d0*bc(iad23)*bc(iad13)                                12d22s19
            shift=shift-term*fact                                        12d22s19
           end do                                                        12d22s19
          end do                                                         12d22s19
         end if                                                          12d22s19
        end if                                                           12d22s19
c
c     fold into effective lzlz.
c     J part is zero becase lz is skew symmetric.
c
        if(multh(n1,n2).eq.islz(ipass).and.n1.eq.n3.and.n2.eq.n4)then   12d31s19
         do i4=0,iacto(n4)-1                                             12d22s19
          i4p=i4+idoub(n4)                                               12d22s19
          do i2=0,iacto(n2)-1                                            12d22s19
           i3p=i2+idoub(n2)                                              12d22s19
           iad2=ixlzze(n2,ione)+i2+iacto(n2)*i4                         12d31s19
           do i13=0,idoub(n1)-1                                          12d22s19
            iad13=ixlzz(n1,ipass)+i13+noc4(n1)*i2p                      12d31s19
            iad24=ixlzz(n3,ipass)+i13+noc4(n3)*i4p                      12d31s19
            term=-2d0*bc(iad13)*bc(iad24)                                12d22s19
            bc(iad2)=bc(iad2)-term                                       12d22s19
           end do                                                        12d22s19
          end do                                                         12d22s19
         end do                                                          12d22s19
        end if                                                           12d22s19
       end do                                                            12d22s19
      end do                                                            12d31s19
      do isb=1,nsymb                                                    12d22s19
       if(iacto(isb).gt.0)then
        do ipass=1,nlzz/2                                                12d31s19
         isk=multh(isb,islz(ipass))                                      12d31s19
         ixlzze(isb,ipass)=ibcoff                                        12d31s19
         ibcoff=ixlzze(isb,ipass)+iacto(isk)*iacto(isb)                  12d31s19
         call enough('corelz2.  2',bc,ibc)
         if(iacto(isk).gt.0)then                                        12d31s19
          do i2=0,iacto(isk)-1                                             12d22s19
           i2p=i2+idoub(isk)                                               12d22s19
           iad1=ixlzze(isb,ipass)+iacto(isb)*i2                           12d31s19
           iad2=ixlzz(isb,ipass)+idoub(isb)+noc4(isb)*i2p                 12d31s19
           do i1=0,iacto(isb)-1                                            12d22s19
            bc(iad1+i1)=bc(iad2+i1)                                        12d22s19
           end do                                                          12d22s19
          end do                                                           12d22s19
         end if                                                         12d31s19
        end do                                                           12d31s19
       end if                                                           12d31s19
      end do                                                            12d22s19
      return                                                            12d22s19
      end                                                               12d22s19
