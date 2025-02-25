c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine buildcasgrad(i4o,ionex,iamat,idoub,iacto,noc,nvirtc,
     $     ipuse,ih0,multh,idwsdeb,iden1,j2den,nsdlkp,isblkp,nsdlk1p,   12d15s17
     $     isblk1p,iamatb,natot,rmssz,idden,iduse,nsblkder,isblkder,    11d9s22
     $     bc,ibc)                                                      11d9s22
      implicit real*8 (a-h,o-z)
c
c     this routine computes the gradient of the energy wrt orbital
c     rotations for a cas wavefunction with fixed densities
c
      include "common.store"
      include "common.hf"
      dimension i4o(*),ionex(*),iamat(8),idoub(8),iacto(8),noc(8),
     $     nvirtc(8),multh(8,8),idoit(5),iden1(8),j2den(*),isblkp(4,*), 12d15s17
     $     isblk1p(4,*),iamatb(8),rmssz(8),isblkder(4,*)                6d16s22
      data idoit/5*1/
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      data icall/0/                                                     5d9s22
      save icall                                                        5d9s22
      loopx=100000
      loop=0
      icall=icall+1                                                     5d9s22
      idoit(1)=idden                                                    5d17s22
      idoit(2)=idden                                                    5d17s22
      irms=0                                                            2d16s18
      natot=0
      do isb=1,nsymb
       iamat(isb)=ibcoff
       jsb=multh(isb,multh(iduse,ipuse))                                6d16s22
       nda=idoub(jsb)*iacto(isb)                                        6d7s22
       nov=noc(jsb)*nvirtc(isb)                                         6d7s22
       ibcoff=iamat(isb)+nda+nov
       natot=natot+nda+nov
      end do
      if(idwsdeb.ne.0)then
       jh0=ih0
       do isb=1,nsymb
        nb=nvirt(isb)+noc(isb)
        if(noc(isb).gt.0)then
         write(6,*)('for symmetry '),isb
         write(6,*)('density: ')
         call prntm2(bc(iden1(isb)),iacto(isb),iacto(isb),iacto(isb))
         write(6,*)('h0 ')
         call prntm2(bc(jh0),nb,nb,nb)
        end if
        jh0=jh0+nb*nb
       end do
       if(nsymb.eq.1)then
        igoal=iamat(1)+1+4*17+6
       else
        igoal=iamat(2)
       end if
       write(6,*)('2e density ')
       do isd=1,nsdlk
        write(6,*)('symmetry type: '),(isblk(j,isd),j=1,4)
        if(isblk(1,isd).eq.isblk(2,isd))then
         nrow=(iacto(isblk(1,isd))*(iacto(isblk(1,isd))+1))/2
        else
         nrow=iacto(isblk(1,isd))*iacto(isblk(2,isd))
        end if
        if(isblk(3,isd).eq.isblk(4,isd))then
         ncol=(iacto(isblk(3,isd))*(iacto(isblk(3,isd))+1))/2
        else
         ncol=iacto(isblk(3,isd))*iacto(isblk(4,isd))
        end if
        call prntm2(bc(j2den(isd)),nrow,ncol,nrow)
        if(nsymb.eq.1)then
         call printa(bc(j2den(isd)),iacto(isblk(1,isd)),
     $        idoub(isblk(1,isd)),iacto(isblk(2,isd)),
     $        idoub(isblk(2,isd)),iacto(isblk(3,isd)),
     $        idoub(isblk(3,isd)),iacto(isblk(4,isd)),
     $        idoub(isblk(4,isd)),bc(ibcoff))
        end if
       end do
      else
       igoal=1
      end if
      call enough('buildcasgrad.  1',bc,ibc)
      do i=0,natot-1
       bc(iamat(1)+i)=0d0
      end do
c
c     dd part
c
      if(idoit(1).ne.0.and.mynowprog.eq.0)then                          5d17s22
       jh0=ih0
       do isb=1,nsymb
        jsb=multh(isb,ipuse)                                            6d7s22
        nb=nvirtc(isb)+noc(isb)
        nj=nvirtc(jsb)+noc(jsb)                                         6d8s22
        jamat=iamat(isb)
        kh0=jh0+idoub(isb)
        do io=0,idoub(jsb)-1                                            6d7s22
         do ia=0,iacto(isb)-1
          bc(jamat+ia)=bc(jamat+ia)+4d0*bc(kh0+ia)
         end do
         jamat=jamat+iacto(isb)
         kh0=kh0+nb
        end do
        kh0=jh0+noc(isb)
        do io=0,idoub(jsb)-1                                            6d7s22
         do iv=0,nvirtc(isb)-1
          bc(jamat+iv)=bc(jamat+iv)+4d0*bc(kh0+iv)
         end do
         jamat=jamat+iv
         kh0=kh0+nb
        end do
        jh0=jh0+nb*nj                                                   6d8s22
       end do
      end if
      if(idoit(2).ne.0)then                                             5d17s22
c
c     dddd j part
c
       do isb=1,nsymb                                                   12d15s17
        jsb=multh(isb,ipuse)                                            6d7s22
        do is=1,nsdlkp                                                    12d15s17
c     do we need jsb isb as well?
         if(isblkp(3,is).eq.isb.and.isblkp(4,is).eq.jsb)then            6d7s22
          is12=isblkp(1,is)                                              12d15s17
          call ilimts(noc(isb),noc(jsb),mynprocg,mynowprog,il,ih,i1s,   6d7s22
     $        i1e,i2s,i2e)                                              12d15s17
          i10=i1s                                                        12d15s17
          i1n=noc(isb)                                                   12d15s17
          if(ipuse.eq.1)then                                            6d7s22
           nrow=(noc(is12)*(noc(is12)+1))/2                               12d15s17
           isw=0                                                        6d7s22
          else                                                          6d7s22
           nrow=noc(is12)*noc(is12)                                     6d7s22
           isw=1                                                        6d7s22
          end if                                                        6d7s22
          ioooo=i4o(is)-1                                               12d15s17
          do i2=i2s,i2e                                                  12d15s17
           if(i2.eq.i2e)i1n=i1e                                          12d15s17
           i2m=i2-1                                                      12d15s17
           do i1=i10,i1n                                                 12d15s17
            if(i2.le.idoub(jsb).and.i1.gt.idoub(isb))then               6d7s22
             i1m=i1-1-idoub(isb)                                         12d15s17
             iad1=iamat(isb)+i1m+iacto(isb)*i2m                          12d15s17
             do i34=1,idoub(is12)                                        12d15s17
              irec=i34+noc(is12)*(i34-1)                                6d7s22
              itri=((i34*(i34+1))/2)                                    6d7s22
              itri=itri+isw*(irec-itri)                                 6d7s22
              iad2=ioooo+itri                                           6d7s22
              bc(iad1)=bc(iad1)+8d0*bc(iad2)                             12d15s17
             end do                                                      12d15s17
            end if                                                       12d15s17
            ioooo=ioooo+nrow                                             12d15s17
           end do                                                        12d15s17
           i10=1                                                         12d15s17
          end do                                                         12d15s17
         end if                                                          12d15s17
        end do                                                           12d15s17
        jamat=iamat(isb)+iacto(isb)*idoub(jsb)                          6d7s22
        do is=1,nsdlk1p                                                  12d15s17
         if(isblk1p(3,is).eq.jsb.and.isblk1p(4,is).eq.isb)then          6d7s22
          is12=isblk1p(1,is)                                             12d15s17
          call ilimts(noc(jsb),nvirtc(isb),mynprocg,mynowprog,il,ih,i1s,6d7s22
     $         i1e,i2s,i2e)                                              12d15s17
          i10=i1s                                                        12d15s17
          i1n=noc(jsb)                                                  6d9s22
          if(ipuse.eq.1)then                                            6d7s22
           nrow=(noc(is12)*(noc(is12)+1))/2                               12d15s17
           isw=0                                                        6d7s22
          else                                                          6d7s22
           nrow=noc(is12)*noc(is12)                                     6d7s22
           isw=1                                                        6d7s22
          end if                                                        6d7s22
          i1x=ionex(is)-1                                               12d15s17
          do i2=i2s,i2e                                                  12d15s17
           if(i2.eq.i2e)i1n=i1e                                          12d15s17
           i2m=i2-1                                                      12d15s17
           do i1=i10,i1n                                                 12d15s17
            if(i1.le.idoub(jsb))then                                     12d15s17
             i1m=i1-1                                                    12d15s17
             iad1=jamat+i2m+nvirtc(isb)*i1m                             12d15s17
             do i34=1,idoub(is12)                                        12d15s17
              irec=i34+noc(is12)*(i34-1)                                6d9s22
              itri=((i34*(i34+1))/2)                                    6d7s22
              itri=itri+isw*(irec-itri)                                 6d7s22
              iad2=i1x+itri                                             6d7s22
              bc(iad1)=bc(iad1)+8d0*bc(iad2)                             12d15s17
             end do                                                      12d15s17
            end if                                                       12d15s17
            i1x=i1x+nrow                                                 12d15s17
           end do                                                        12d15s17
           i10=1                                                         12d15s17
          end do                                                         12d15s17
         end if                                                          12d15s17
        end do                                                           12d15s17
c
c     k part
c
        do is=1,nsdlkp                                                  12d15s17
         if(isblkp(1,is).eq.jsb.and.isblkp(3,is).eq.isb)then            12d15s17
          is24=isblkp(2,is)                                             12d15s17
          call ilimts(noc(isb),noc(is24),mynprocg,mynowprog,il,ih,i1s,  6d7s22
     $         i1e,i2s,i2e)                                             12d15s17
          i10=i1s                                                       12d15s17
          i1n=noc(isb)                                                  12d15s17
          if(isb.eq.is24.and.ipuse.eq.1)then                            6d7s22
           nrow=(noc(isb)*(noc(isb)+1))/2                               12d15s17
          else                                                          12d15s17
           nrow=noc(jsb)*noc(is24)                                      6d8s22
          end if                                                        12d15s17
          ioooo=i4o(is)                                                 12d15s17
          do i2=i2s,i2e                                                 12d15s17
           i2m=i2-1                                                     12d15s17
           if(i2.eq.i2e)i1n=i1e                                         12d15s17
           do i1=i10,i1n                                                12d15s17
            if(i2.le.idoub(is24).and.i1.gt.idoub(isb))then              12d15s17
             i1m=i1-1-idoub(isb)
             if(isb.eq.is24.and.ipuse.eq.1)then                         6d7s22
              do id=0,idoub(isb)-1                                       12d15s17
               ix=max(id,i2m)                                           12d15s17
               in=min(id,i2m)                                           12d15s17
               iad1=ioooo+((ix*(ix+1))/2)+in                            12d15s17
               iad2=iamat(isb)+i1m+iacto(isb)*id                        12d15s17
               bc(iad2)=bc(iad2)-4d0*bc(iad1)                           12d15s17
              end do                                                     12d15s17
             else
              do id=0,idoub(jsb)-1                                      6d7s22
               iad1=ioooo+id+noc(jsb)*i2m                               6d7s22
               iad2=iamat(isb)+i1m+iacto(isb)*id                        12d15s17
               bc(iad2)=bc(iad2)-4d0*bc(iad1)                           12d15s17
              end do                                                     12d15s17
             end if
            end if                                                      12d15s17
            ioooo=ioooo+nrow                                            12d15s17
           end do                                                       12d15s17
           i10=1                                                        12d15s17
          end do                                                        12d15s17
         else if(isblkp(1,is).eq.jsb.and.isblkp(4,is).eq.isb)then       6d7s22
          is23=isblkp(2,is)                                             12d15s17
          call ilimts(noc(is23),noc(isb),mynprocg,mynowprog,il,ih,i1s,  12d15s17
     $         i1e,i2s,i2e)                                             12d15s17
          i10=i1s                                                       12d15s17
          i1n=noc(is23)                                                 12d15s17
          nrow=noc(jsb)*noc(is23)                                       6d7s22
          ioooo=i4o(is)                                                 12d15s17
          do i2=i2s,i2e                                                 12d15s17
           i2m=i2-1-idoub(isb)                                          12d15s17
           if(i2.eq.i2e)i1n=i1e                                         12d15s17
           do i1=i10,i1n                                                12d15s17
            if(i1.le.idoub(is23).and.i2.gt.idoub(isb))then              12d15s17
             i1m=i1-1                                                   12d15s17
             do id=0,idoub(jsb)-1                                       6d7s22
              iad1=ioooo+id+noc(jsb)*i1m                                6d7s22
              iad2=iamat(isb)+i2m+iacto(isb)*id                         12d15s17
              bc(iad2)=bc(iad2)-4d0*bc(iad1)                            12d15s17
             end do                                                     12d15s17
            end if                                                      12d15s17
            ioooo=ioooo+nrow                                            12d15s17
           end do                                                       12d15s17
           i10=1                                                        12d15s17
          end do                                                        12d15s17
         else if(isblkp(2,is).eq.jsb.and.isblkp(4,is).eq.isb)then       6d7s22
          is13=isblkp(1,is)                                             12d15s17
          call ilimts(noc(is13),noc(isb),mynprocg,mynowprog,il,ih,i1s,  12d15s17
     $         i1e,i2s,i2e)                                             12d15s17
          i10=i1s                                                       12d15s17
          i1n=noc(is13)                                                 12d15s17
          nrow=noc(jsb)*noc(is13)                                       6d7s22
          ioooo=i4o(is)                                                 12d15s17
          do i2=i2s,i2e                                                 12d15s17
           i2m=i2-1-idoub(isb)                                          12d15s17
           if(i2.eq.i2e)i1n=i1e                                         12d15s17
           do i1=i10,i1n                                                12d15s17
            if(i1.le.idoub(is13).and.i2.gt.idoub(isb))then              12d15s17
             i1m=i1-1                                                   12d15s17
             do id=0,idoub(jsb)-1                                       6d7s22
              iad1=ioooo+i1m+noc(is13)*id                               12d15s17
              iad2=iamat(isb)+i2m+iacto(isb)*id                         12d15s17
              bc(iad2)=bc(iad2)-4d0*bc(iad1)                            12d15s17
             end do                                                     12d15s17
            end if                                                      12d15s17
            ioooo=ioooo+nrow                                            12d15s17
           end do                                                       12d15s17
           i10=1                                                        12d15s17
          end do                                                        12d15s17
         else if(isblkp(2,is).eq.jsb.and.isblkp(3,is).eq.isb)then       6d7s22
          is14=isblkp(1,is)                                             12d15s17
          call ilimts(noc(isb),noc(is14),mynprocg,mynowprog,il,ih,i1s,  12d15s17
     $         i1e,i2s,i2e)                                             12d15s17
          i10=i1s                                                       12d15s17
          i1n=noc(isb)                                                  12d15s17
          nrow=noc(jsb)*noc(is14)                                       6d7s22
          ioooo=i4o(is)                                                 12d15s17
          do i2=i2s,i2e                                                 12d15s17
           i2m=i2-1                                                     12d15s17
           if(i2.eq.i2e)i1n=i1e                                         12d15s17
           do i1=i10,i1n                                                12d15s17
            if(i2.le.idoub(is14).and.i1.gt.idoub(isb))then              12d15s17
             i1m=i1-1-idoub(isb)                                        12d15s17
             do id=0,idoub(jsb)-1                                       6d7s22
              iad1=ioooo+i2m+noc(is14)*id                               12d15s17
              iad2=iamat(isb)+i1m+iacto(isb)*id                         12d15s17
              bc(iad2)=bc(iad2)-4d0*bc(iad1)                            12d15s17
             end do                                                     12d15s17
            end if                                                      12d15s17
            ioooo=ioooo+nrow                                            12d15s17
           end do                                                       12d15s17
           i10=1                                                        12d15s17
          end do                                                        12d15s17
         end if                                                         12d15s17
        end do                                                          12d15s17
        jamat=iamat(isb)+iacto(isb)*idoub(jsb)                          6d7s22
        do is=1,nsdlk1p                                                 12d15s17
         if(isblk1p(1,is).eq.jsb.and.isblk1p(4,is).eq.isb)then          6d7s22
          is23=isblk1p(2,is)                                            12d15s17
          call ilimts(noc(is23),nvirtc(isb),mynprocg,mynowprog,il,ih,   12d15s17
     $         i1s,i1e,i2s,i2e)                                         12d15s17
          i10=i1s                                                       12d15s17
          i1n=noc(is23)                                                 12d15s17
          if(isb.eq.is23.and.ipuse.eq.1)then                            6d7s22
           nrow=(noc(isb)*(noc(isb)+1))/2                               12d15s17
           iswitch=0                                                    12d16s17
          else                                                          12d15s17
           nrow=noc(jsb)*noc(is23)                                      6d9s22
           iswitch=1                                                    12d16s17
          end if                                                        12d15s17
          i1x=ionex(is)                                                 12d16s17
          do i2=i2s,i2e                                                 12d16s17
           if(i2.eq.i2e)i1n=i1e                                         12d16s17
           i2m=i2-1                                                     12d16s17
           do i1=i10,i1n                                                12d16s17
            if(i1.le.idoub(is23))then                                   12d16s17
             i1m=i1-1                                                   12d16s17
             do i3=0,idoub(jsb)-1                                       12d16s17
              ix=max(i3,i1m)                                            12d16s17
              in=min(i3,i1m)                                            12d16s17
              ieq=((ix*(ix+1))/2)+in                                    12d16s17
              inot=i3+noc(jsb)*i1m                                      6d7s22
              iad1=i1x+(inot-ieq)*iswitch+ieq                           12d16s17
              iad2=jamat+i2m+nvirtc(isb)*i3                             12d16s17
              bc(iad2)=bc(iad2)-4d0*bc(iad1)                            12d16s17
             end do                                                     12d16s17
            end if                                                      12d16s17
            i1x=i1x+nrow                                                12d16s17
           end do                                                       12d16s17
           i10=1                                                        12d16s17
          end do                                                        12d16s17
         else if(isblk1p(2,is).eq.jsb.and.isblk1p(4,is).eq.isb)then     6d7s22
          is13=isblk1p(1,is)                                            12d16s17
          call ilimts(noc(is13),nvirtc(isb),mynprocg,mynowprog,il,ih,   12d15s17
     $         i1s,i1e,i2s,i2e)                                         12d15s17
          i10=i1s                                                       12d15s17
          i1n=noc(is13)                                                 12d15s17
          nrow=noc(jsb)*noc(is13)                                       6d7s22
          i1x=ionex(is)                                                 12d16s17
          do i2=i2s,i2e                                                 12d16s17
           if(i2.eq.i2e)i1n=i1e                                         12d16s17
           i2m=i2-1                                                     12d16s17
           do i1=i10,i1n                                                12d16s17
            if(i1.le.idoub(is13))then                                   12d16s17
             i1m=i1-1                                                   12d16s17
             do i4=0,idoub(jsb)-1                                       6d7s22
              inot=i3+noc(isb)*i1m                                      12d16s17
              iad1=i1x+i1m+noc(is13)*i4                                 12d16s17
              iad2=jamat+i2m+nvirtc(isb)*i4                             12d16s17
              bc(iad2)=bc(iad2)-4d0*bc(iad1)                            12d16s17
             end do                                                     12d16s17
            end if                                                      12d16s17
            i1x=i1x+nrow                                                12d16s17
           end do                                                       12d16s17
           i10=1                                                        12d16s17
          end do                                                        12d16s17
         end if                                                         12d15s17
        end do                                                          12d15s17
       end do                                                           12d15s17
      end if                                                            12d15s17
      if(idoit(4).ne.0)then                                             5d17s22
c
c     ddaa part
c
       do isb=1,nsymb                                                   12d16s17
        jsb=multh(isb,ipuse)                                            6d16s22
        ksb=multh(isb,iduse)                                            6d16s22
        do is=1,nsdlkp                                                  12d16s17
         if(iduse.eq.1)then                                             6d16s22
          if(isblkp(1,is).eq.isb.and.isblkp(2,is).eq.jsb)then            6d8s22
           is34=isblkp(3,is)                                             12d16s17
           call ilimts(noc(is34),noc(is34),mynprocg,mynowprog,il,ih,     12d16s17
     $         i1s,i1e,i2s,i2e)                                         12d16s17
           i10=i1s                                                       12d16s17
           i1n=noc(is34)                                                 12d16s17
           if(jsb.eq.isb)then                                            6d7s22
            nrow=(noc(isb)*(noc(isb)+1))/2                                12d16s17
            isw=0                                                        6d7s22
           else                                                          6d7s22
            nrow=noc(isb)*noc(jsb)                                       6d7s22
            isw=1                                                        6d7s22
           end if                                                        6d7s22
           ioooo=i4o(is)                                                 12d16s17
           do i2=i2s,i2e                                                 12d16s17
            if(i2.eq.i2e)i1n=i1e                                         12d16s17
            i2m=i2-1-idoub(is34)                                         12d16s17
            do i1=i10,i1n                                                12d16s17
             if(i1.gt.idoub(is34).and.i2.gt.idoub(is34))then             6d8s22
              i1m=i1-1-idoub(is34)                                       12d16s17
              jden=iden1(is34)+i1m+iacto(is34)*i2m                       12d16s17
              fact=4d0*bc(jden)                                          12d16s17
              do i4=0,idoub(jsb)-1                                       6d7s22
               do i3=0,iacto(isb)-1                                      12d16s17
                i3p=i3+idoub(isb)                                        12d16s17
                irec=i3p+noc(isb)*i4                                     6d7s22
                ix=max(i4,i3p)                                           12d16s17
                in=min(i4,i3p)                                           12d16s17
                itri=((ix*(ix+1))/2)+in                                  6d7s22
                itri=itri+isw*(irec-itri)                                6d7s22
                iad1=ioooo+itri                                          6d7s22
                iad2=iamat(isb)+i3+iacto(isb)*i4                         12d16s17
                bc(iad2)=bc(iad2)+fact*bc(iad1)                          12d16s17
               end do                                                    12d16s17
              end do                                                     12d16s17
             else if(i1.le.idoub(is34).and.i1.eq.i2)then                 6d8s22
              do i4=0,iacto(isb)-1                                       12d16s17
               i4p=i4+idoub(isb)                                         12d17s17
               jden=iden1(isb)+iacto(isb)*i4                             12d17s17
               do i3=0,idoub(jsb)-1                                      6d7s22
                irec=i4p+noc(isb)*i3                                     6d8s22
                ix=max(i3,i4p)                                           12d16s17
                in=min(i3,i4p)                                           12d16s17
                itri=((ix*(ix+1))/2)+in                                  6d7s22
                iad1=ioooo+itri                                          6d7s22
                fact=4d0*bc(iad1)                                        12d16s17
                jamat=iamat(isb)+iacto(isb)*i3                           12d16s17
                do ia=0,iacto(isb)-1                                     12d16s17
                 bc(jamat+ia)=bc(jamat+ia)-fact*bc(jden+ia)              12d16s17
                end do                                                   12d16s17
               end do                                                    12d16s17
              end do                                                     12d16s17
             end if                                                      12d16s17
             ioooo=ioooo+nrow                                            12d16s17
            end do                                                       12d16s17
            i10=1                                                        12d17s17
           end do                                                        12d16s17
          end if                                                         12d16s17
          if((isblkp(1,is).eq.jsb.and.isblkp(3,is).eq.isb).or.           6d7s22
     $      (isblkp(1,is).eq.isb.and.isblkp(3,is).eq.jsb))then          6d7s22
           is24=isblkp(2,is)                                             12d16s17
           if(isblkp(3,is).eq.isb)then                                   6d7s22
            call ilimts(noc(isb),noc(is24),mynprocg,mynowprog,il,ih,      12d16s17
     $         i1s,i1e,i2s,i2e)                                         12d16s17
            i10=i1s                                                       12d16s17
            i1n=noc(isb)                                                  12d16s17
            if(isb.eq.is24.and.ipuse.eq.1)then                            6d7s22
             nrow=(noc(isb)*(noc(isb)+1))/2                               12d16s17
             iswitch=0                                                    12d16s17
            else                                                          12d16s17
             nrow=noc(jsb)*noc(is24)                                      6d7s22
             iswitch=1                                                    12d16s17
            end if                                                        12d16s17
            ioooo=i4o(is)                                                 12d16s17
            do i2=i2s,i2e                                                 12d16s17
             if(i2.eq.i2e)i1n=i1e                                         12d16s17
             i2m=i2-1                                                     12d16s17
             do i1=i10,i1n                                                12d16s17
              if(i1.gt.idoub(isb).and.i2.le.idoub(is24))then              12d16s17
               i1m=i1-1-idoub(isb)                                        12d16s17
               jden=iden1(isb)+iacto(isb)*i1m                             12d16s17
               do i4=0,idoub(jsb)-1                                       6d7s22
                ix=max(i4,i2m)                                            12d16s17
                in=min(i4,i2m)                                            12d16s17
                ieq=((ix*(ix+1))/2)+in                                    12d16s17
                inot=i4+noc(jsb)*i2m                                      6d7s22
                iad1=ioooo+(inot-ieq)*iswitch+ieq                         12d16s17
                iad2=iamat(isb)+iacto(isb)*i4                             12d17s17
                fact=-2d0*bc(iad1)                                        12d16s17
                do ia=0,iacto(isb)-1                                      12d16s17
                 bc(iad2+ia)=bc(iad2+ia)-fact*bc(jden+ia)                 12d16s17
                end do                                                    12d16s17
               end do                                                     12d16s17
              else if(i1.le.idoub(isb).and.i2.gt.idoub(is24).and.        6d7s22
     $             ipuse.eq.1)then                                      6d7s22
               i1m=i1-1                                                   12d16s17
               iad3=iamat(isb)+iacto(isb)*i1m                             12d16s17
               do i4=0,iacto(is24)-1                                      12d16s17
                iad2=iden1(is24)+i4+iacto(is24)*(i2m-idoub(is24))         12d16s17
                fact=-2d0*bc(iad2)                                        12d16s17
                i4p=i4+idoub(is24)                                        12d16s17
                do i3=0,iacto(isb)-1                                      12d16s17
                 i3p=i3+idoub(isb)                                        12d16s17
                 ix=max(i3p,i4p)                                          12d16s17
                 in=min(i3p,i4p)                                          12d16s17
                 ieq=((ix*(ix+1))/2)+in                                   12d16s17
                 inot=i3p+noc(isb)*i4p                                    12d16s17
                 iad1=ioooo+(inot-ieq)*iswitch+ieq                        12d16s17
                 bc(iad3+i3)=bc(iad3+i3)+fact*bc(iad1)                    12d16s17
                end do                                                    12d16s17
               end do                                                     12d16s17
              end if                                                      12d16s17
              ioooo=ioooo+nrow                                            12d16s17
             end do                                                       12d16s17
             i10=1                                                        12d16s17
            end do                                                        12d16s17
           else                                                          6d7s22
            call ilimts(noc(jsb),noc(is24),mynprocg,mynowprog,il,ih,     6d7s22
     $         i1s,i1e,i2s,i2e)                                         12d16s17
            i10=i1s                                                       12d16s17
            i1n=noc(jsb)                                                 6d7s22
            nrow=noc(isb)*noc(is24)                                      6d7s22
            ioooo=i4o(is)                                                 12d16s17
            do i2=i2s,i2e                                                 12d16s17
             if(i2.eq.i2e)i1n=i1e                                         12d16s17
             i2m=i2-1                                                     12d16s17
             do i1=i10,i1n                                                12d16s17
              if(i1.le.idoub(jsb).and.i2.gt.idoub(is24))then             6d7s22
               i1m=i1-1                                                   12d16s17
               iad3=iamat(isb)+iacto(isb)*i1m                             12d16s17
               do i4=0,iacto(is24)-1                                      12d16s17
                iad2=iden1(is24)+i4+iacto(is24)*(i2m-idoub(is24))         12d16s17
                fact=-2d0*bc(iad2)                                        12d16s17
                i4p=i4+idoub(is24)                                        12d16s17
                do i3=0,iacto(isb)-1                                      12d16s17
                 i3p=i3+idoub(isb)                                        12d16s17
                 inot=i3p+noc(isb)*i4p                                    12d16s17
                 iad1=ioooo+inot                                         6d7s22
                 bc(iad3+i3)=bc(iad3+i3)+fact*bc(iad1)                    12d16s17
                end do                                                    12d16s17
               end do                                                     12d16s17
              end if                                                      12d16s17
              ioooo=ioooo+nrow                                            12d16s17
             end do                                                       12d16s17
             i10=1                                                        12d16s17
            end do                                                        12d16s17
           end if                                                        6d7s22
          else if((isblkp(1,is).eq.jsb.and.isblkp(4,is).eq.isb).or.      6d7s22
     $           (isblkp(1,is).eq.isb.and.isblkp(4,is).eq.jsb))then     6d7s22
           is23=isblkp(2,is)                                             12d16s17
           if(isblkp(1,is).eq.jsb.and.isblkp(4,is).eq.isb)then           6d7s22
            call ilimts(noc(is23),noc(isb),mynprocg,mynowprog,il,ih,      12d19s17
     $         i1s,i1e,i2s,i2e)                                         12d16s17
            i10=i1s                                                       12d16s17
            i1n=noc(is23)                                                 12d19s17
            nrow=noc(jsb)*noc(is23)                                       6d7s22
            ioooo=i4o(is)                                                 12d16s17
            do i2=i2s,i2e                                                 12d16s17
             if(i2.eq.i2e)i1n=i1e                                         12d16s17
             i2m=i2-1-idoub(isb)                                          12d19s17
             i2n=i2-1                                                     12d19s17
             do i1=i10,i1n                                                12d16s17
              if(i2.gt.idoub(isb).and.i1.le.idoub(is23))then              12d19s17
               i1m=i1-1                                                   12d19s17
               jden=iden1(isb)+iacto(isb)*i2m                             12d19s17
               do i4=0,idoub(jsb)-1                                       12d16s17
                iad1=ioooo+i4+noc(jsb)*i1m                                12d19s17
                iad2=iamat(isb)+iacto(isb)*i4                             12d17s17
                fact=-2d0*bc(iad1)                                        12d16s17
                do ia=0,iacto(isb)-1                                      12d16s17
                 bc(iad2+ia)=bc(iad2+ia)-fact*bc(jden+ia)                 12d16s17
                end do                                                    12d16s17
               end do                                                     12d16s17
              else if(i2.le.idoub(isb).and.i1.gt.idoub(is23).and.        6d7s22
     $             ipuse.eq.1)then                                      6d7s22
               i1m=i1-1-idoub(is23)                                       12d19s17
               iad3=iamat(isb)+iacto(isb)*i2n                             12d19s17
               do i4=0,iacto(is23)-1                                      12d16s17
                iad2=iden1(is23)+i4+iacto(is23)*i1m                       12d19s17
                fact=-2d0*bc(iad2)                                        12d16s17
                i4p=i4+idoub(is23)                                        12d16s17
                do i3=0,iacto(isb)-1                                      12d16s17
                 i3p=i3+idoub(isb)                                        12d16s17
                 iad1=ioooo+i3p+noc(isb)*i4p                              12d19s17
                 bc(iad3+i3)=bc(iad3+i3)+fact*bc(iad1)                    12d16s17
                end do                                                    12d16s17
               end do                                                     12d16s17
              end if                                                      12d16s17
              ioooo=ioooo+nrow                                            12d16s17
             end do                                                       12d16s17
             i10=1                                                        12d16s17
            end do                                                        12d16s17
           else                                                          6d7s22
            call ilimts(noc(is23),noc(jsb),mynprocg,mynowprog,il,ih,     6d7s22
     $         i1s,i1e,i2s,i2e)                                         12d16s17
            i10=i1s                                                       12d16s17
            i1n=noc(is23)                                                 12d19s17
            nrow=noc(isb)*noc(is23)                                       6d7s22
            ioooo=i4o(is)                                                 12d16s17
            do i2=i2s,i2e                                                 12d16s17
             if(i2.eq.i2e)i1n=i1e                                         12d16s17
             i2m=i2-1-idoub(jsb)                                         6d7s22
             i2n=i2-1                                                     12d19s17
             do i1=i10,i1n                                                12d16s17
              if(i2.le.idoub(jsb).and.i1.gt.idoub(is23))then             6d7s22
               i1m=i1-1-idoub(is23)                                       12d19s17
               iad3=iamat(isb)+iacto(isb)*i2n                             12d19s17
               do i4=0,iacto(is23)-1                                      12d16s17
                iad2=iden1(is23)+i4+iacto(is23)*i1m                       12d19s17
                fact=-2d0*bc(iad2)                                        12d16s17
                i4p=i4+idoub(is23)                                        12d16s17
                do i3=0,iacto(isb)-1                                      12d16s17
                 i3p=i3+idoub(isb)                                        12d16s17
                 iad1=ioooo+i3p+noc(isb)*i4p                              12d19s17
                 bc(iad3+i3)=bc(iad3+i3)+fact*bc(iad1)                    12d16s17
                end do                                                    12d16s17
               end do                                                     12d16s17
              end if                                                      12d16s17
              ioooo=ioooo+nrow                                            12d16s17
             end do                                                       12d16s17
             i10=1                                                        12d16s17
            end do                                                        12d16s17
           end if                                                        6d7s22
          else if((isblkp(2,is).eq.jsb.and.isblkp(4,is).eq.isb).or.      6d7s22
     $           (isblkp(2,is).eq.isb.and.isblkp(4,is).eq.jsb))then     6d7s22
           is13=isblkp(1,is)                                             12d19s17
           if(isblkp(2,is).eq.jsb.and.isblkp(4,is).eq.isb)then           6d7s22
            call ilimts(noc(is13),noc(isb),mynprocg,mynowprog,il,ih,      12d19s17
     $         i1s,i1e,i2s,i2e)                                         12d16s17
            i10=i1s                                                       12d16s17
            i1n=noc(is13)                                                 12d19s17
            nrow=noc(jsb)*noc(is13)                                       12d19s17
            ioooo=i4o(is)                                                 12d16s17
            do i2=i2s,i2e                                                 12d16s17
             if(i2.eq.i2e)i1n=i1e                                         12d16s17
             i2m=i2-1-idoub(isb)                                          12d19s17
             i2n=i2-1                                                     12d19s17
             do i1=i10,i1n                                                12d16s17
              if(i2.gt.idoub(isb).and.i1.le.idoub(is13))then              12d19s17
               i1m=i1-1                                                   12d19s17
               jden=iden1(isb)+iacto(isb)*i2m                             12d19s17
               do i4=0,idoub(jsb)-1                                      6d7s22
                iad1=ioooo+i1m+noc(is13)*i4                                12d19s17
                iad2=iamat(isb)+iacto(isb)*i4                             12d17s17
                fact=-2d0*bc(iad1)                                        12d16s17
                do ia=0,iacto(isb)-1                                      12d16s17
                 bc(iad2+ia)=bc(iad2+ia)-fact*bc(jden+ia)                 12d16s17
                end do                                                    12d16s17
               end do                                                     12d16s17
              else if(i2.le.idoub(isb).and.i1.gt.idoub(is13).and.        6d7s22
     $             ipuse.eq.1)then                                      6d7s22
               i1m=i1-1-idoub(is13)                                       12d19s17
               iad3=iamat(isb)+iacto(isb)*i2n                             12d19s17
               do i4=0,iacto(is13)-1                                      12d19s17
                iad2=iden1(is13)+i4+iacto(is13)*i1m                       12d19s17
                fact=-2d0*bc(iad2)                                        12d16s17
                i4p=i4+idoub(is13)                                        12d19s17
                do i3=0,iacto(isb)-1                                      12d16s17
                 i3p=i3+idoub(isb)                                        12d16s17
                 iad1=ioooo+i4p+noc(is13)*i3p                             12d19s17
                 bc(iad3+i3)=bc(iad3+i3)+fact*bc(iad1)                    12d16s17
                end do                                                    12d16s17
               end do                                                     12d16s17
              end if                                                      12d16s17
              ioooo=ioooo+nrow                                            12d16s17
             end do                                                       12d16s17
             i10=1                                                        12d16s17
            end do                                                        12d16s17
           else                                                          6d7s22
            call ilimts(noc(is13),noc(jsb),mynprocg,mynowprog,il,ih,      12d19s17
     $         i1s,i1e,i2s,i2e)                                         12d16s17
            i10=i1s                                                       12d16s17
            i1n=noc(is13)                                                 12d19s17
            nrow=noc(isb)*noc(is13)                                       12d19s17
            ioooo=i4o(is)                                                 12d16s17
            do i2=i2s,i2e                                                 12d16s17
             if(i2.eq.i2e)i1n=i1e                                         12d16s17
             i2m=i2-1-idoub(jsb)                                          12d19s17
             i2n=i2-1                                                     12d19s17
             do i1=i10,i1n                                                12d16s17
              if(i2.le.idoub(jsb).and.i1.gt.idoub(is13))then             6d7s22
               i1m=i1-1-idoub(is13)                                       12d19s17
               iad3=iamat(isb)+iacto(isb)*i2n                             12d19s17
               do i4=0,iacto(is13)-1                                      12d19s17
                iad2=iden1(is13)+i4+iacto(is13)*i1m                       12d19s17
                fact=-2d0*bc(iad2)                                        12d16s17
                i4p=i4+idoub(is13)                                        12d19s17
                do i3=0,iacto(isb)-1                                      12d16s17
                 i3p=i3+idoub(isb)                                        12d16s17
                 iad1=ioooo+i4p+noc(is13)*i3p                             12d19s17
                 bc(iad3+i3)=bc(iad3+i3)+fact*bc(iad1)                    12d16s17
                end do                                                    12d16s17
               end do                                                     12d16s17
              end if                                                      12d16s17
              ioooo=ioooo+nrow                                            12d16s17
             end do                                                       12d16s17
             i10=1                                                        12d16s17
            end do                                                        12d16s17
           end if                                                        6d7s22
          else if((isblkp(2,is).eq.jsb.and.isblkp(3,is).eq.isb).and.     6d7s22
     $           (isblkp(2,is).eq.isb.and.isblkp(3,is).eq.isb))then     6d7s22
           is14=isblkp(1,is)                                             12d19s17
           if(isblkp(2,is).eq.jsb.and.isblkp(3,is).eq.isb)then           6d7s22
            call ilimts(noc(isb),noc(is14),mynprocg,mynowprog,il,ih,      12d19s17
     $         i1s,i1e,i2s,i2e)                                         12d16s17
            i10=i1s                                                       12d16s17
            i1n=noc(isb)                                                  12d19s17
            nrow=noc(jsb)*noc(is14)                                       12d19s17
            ioooo=i4o(is)                                                 12d16s17
            do i2=i2s,i2e                                                 12d16s17
             if(i2.eq.i2e)i1n=i1e                                         12d16s17
             i2m=i2-1                                                     12d19s17
             do i1=i10,i1n                                                12d16s17
              if(i1.gt.idoub(isb).and.i2.le.idoub(is14))then              12d19s17
               i1m=i1-1-idoub(isb)                                        12d19s17
               jden=iden1(isb)+iacto(isb)*i1m                             12d19s17
               do i4=0,idoub(jsb)-1                                      6d7s22
                iad1=ioooo+i2m+noc(is14)*i4                               12d19s17
                iad2=iamat(isb)+iacto(isb)*i4                             12d17s17
                fact=-2d0*bc(iad1)                                        12d16s17
                do ia=0,iacto(isb)-1                                      12d16s17
                 bc(iad2+ia)=bc(iad2+ia)-fact*bc(jden+ia)                 12d16s17
                end do                                                    12d16s17
               end do                                                     12d16s17
              else if(i1.le.idoub(isb).and.i2.gt.idoub(is14).and.        6d7s22
     $              ipuse.eq.1)then                                      6d7s22
               i1m=i1-1                                                  6d7s22
               iad3=iamat(isb)+iacto(isb)*i1m                             12d19s17
               do i4=0,iacto(is14)-1                                      12d19s17
                iad2=iden1(is14)+i4+iacto(is14)*(i2m-idoub(is14))         12d19s17
                fact=-2d0*bc(iad2)                                        12d16s17
                i4p=i4+idoub(is14)                                        12d19s17
                do i3=0,iacto(isb)-1                                      12d16s17
                 i3p=i3+idoub(isb)                                        12d16s17
                 iad1=ioooo+i4p+noc(is14)*i3p                             12d19s17
                 bc(iad3+i3)=bc(iad3+i3)+fact*bc(iad1)                    12d16s17
                end do                                                    12d16s17
               end do                                                     12d16s17
              end if                                                      12d16s17
              ioooo=ioooo+nrow                                            12d16s17
             end do                                                       12d16s17
             i10=1                                                        12d16s17
            end do                                                        12d16s17
           else                                                          6d7s22
            call ilimts(noc(jsb),noc(is14),mynprocg,mynowprog,il,ih,     6d7s22
     $           i1s,i1e,i2s,i2e)                                         12d16s17
            i10=i1s                                                       12d16s17
            i1n=noc(jsb)                                                 6d7s22
            nrow=noc(isb)*noc(is14)                                      6d7s22
            ioooo=i4o(is)                                                 12d16s17
            do i2=i2s,i2e                                                 12d16s17
             if(i2.eq.i2e)i1n=i1e                                         12d16s17
             i2m=i2-1                                                     12d19s17
             do i1=i10,i1n                                                12d16s17
              if(i1.le.idoub(jsb).and.i2.gt.idoub(is14))then             6d7s22
               i1m=i1-1                                                  6d7s22
               iad3=iamat(isb)+iacto(isb)*i1m                             12d19s17
               do i4=0,iacto(is14)-1                                      12d19s17
                iad2=iden1(is14)+i4+iacto(is14)*(i2m-idoub(is14))         12d19s17
                fact=-2d0*bc(iad2)                                        12d16s17
                i4p=i4+idoub(is14)                                        12d19s17
                do i3=0,iacto(isb)-1                                      12d16s17
                 i3p=i3+idoub(isb)                                        12d16s17
                 iad1=ioooo+i4p+noc(is14)*i3p                             12d19s17
                 bc(iad3+i3)=bc(iad3+i3)+fact*bc(iad1)                    12d16s17
                end do                                                    12d16s17
               end do                                                     12d16s17
              end if                                                      12d16s17
              ioooo=ioooo+nrow                                            12d16s17
             end do                                                       12d16s17
             i10=1                                                        12d16s17
            end do                                                        12d16s17
           end if                                                        6d7s22
          end if                                                         12d16s17
         else                                                           6d16s22
c     Aad=4 sum (ad|a'a")Da'a" : l
          if(isblkp(1,is).eq.isb.and.isblkp(2,is).eq.ksb)then           6d17s22
           call ilimts(noc(isblkp(3,is)),noc(isblkp(4,is)),mynprocg,    6d17s22
     $         mynowprog,il,ih,i1s,i1e,i2s,i2e)                         6d17s22
           i10=i1s                                                      6d17s22
           i1n=noc(isblkp(3,is))                                        6d17s22
           nrow=noc(isb)*noc(ksb)                                       6d17s22
           ioooo=i4o(is)                                                6d17s22
           do i2=i2s,i2e                                                6d17s22
            if(i2.eq.i2e)i1n=i1e                                        6d17s22
            do i1=i10,i1n                                               6d17s22
             if(i1.gt.idoub(isblkp(3,is)).and.                          6d17s22
     $            i2.gt.idoub(isblkp(4,is)))then                        6d17s22
              i1m=i1-1-idoub(isblkp(3,is))                              6d17s22
              i2m=i2-1-idoub(isblkp(4,is))                              6d17s22
              iad1=iden1(isblkp(3,is))+i1m+iacto(isblkp(3,is))*i2m      6d17s22
              fact=4d0*bc(iad1)*2d0                                     6d17s22
              do i4=0,idoub(ksb)-1                                      6d17s22
               iad2=iamat(isb)+iacto(isb)*i4                            6d17s22
               iad3=ioooo+iacto(isb)*i4                                 6d17s22
               do i3=0,iacto(isb)-1                                     6d17s22
                bc(iad2+i3)=bc(iad2+i3)+fact*bc(iad3+i3)                6d17s22
               end do                                                   6d17s22
              end do                                                    6d17s22
             end if                                                     6d17s22
             ioooo=ioooo+nrow                                           6d17s22
            end do                                                      6d17s22
            i10=1                                                       6d17s22
           end do                                                       6d17s22
          end if                                                        6d16s22
c     Aad=-4 sum (da'|d'd') Daa' : m
          if(isblkp(3,is).eq.isblkp(4,is).and.isblkp(2,is).eq.ksb)then  6d17s22
           call ilimts(noc(isblkp(3,is)),noc(isblkp(4,is)),mynprocg,    6d17s22
     $         mynowprog,il,ih,i1s,i1e,i2s,i2e)                         6d17s22
           nrow=(noc(ksb)*(noc(ksb)+1))/2                               6d17s22
           i10=i1s                                                      6d17s22
           i1n=noc(isblkp(3,is))                                        6d17s22
           ioooo=i4o(is)                                                6d17s22
           do i2=i2s,i2e                                                6d17s22
            if(i2.eq.i2e)i1n=i1e                                        6d17s22
            do i1=i10,i1n                                               6d17s22
             if(i1.eq.i2.and.i1.le.idoub(isblkp(3,is)))then             6d17s22
              do i4=0,iacto(ksb)-1                                       6d17s22
               i4p=i4+idoub(ksb)                                        6d17s22
               iad1=ioooo+((i4p*(i4p+1))/2)                             6d17s22
               iad3=iden1(isb)+iacto(isb)*i4                            6d17s22
               do i3=0,idoub(ksb)-1                                     6d17s22
                iad2=iamat(isb)+iacto(isb)*i3                           6d17s22
                fact=-bc(iad1+i3)*4d0                                   6d17s22
                do ia=0,iacto(isb)-1                                    6d17s22
                 bc(iad2+ia)=bc(iad2+ia)+fact*bc(iad3+ia)               7d1s22
                end do                                                  6d17s22
               end do                                                   6d17s22
              end do                                                    6d17s22
             end if                                                     6d17s22
             ioooo=ioooo+nrow                                           6d17s22
            end do                                                      6d17s22
            i10=1                                                       6d17s22
           end do                                                       6d17s22
          end if
c     Aad= 2 sum (dd'|a'd') Daa' : n
c      ik         k   k      ik
          if(isblkp(1,is).eq.ksb.and.isblkp(3,is).eq.ksb)then           6d17s22
           call ilimts(noc(isblkp(3,is)),noc(isblkp(4,is)),mynprocg,    6d17s22
     $         mynowprog,il,ih,i1s,i1e,i2s,i2e)                         6d17s22
           if(isblkp(1,is).ne.isblkp(2,is))then                         6d17s22
            nrow=noc(isblkp(1,is))*noc(isblkp(2,is))                     6d17s22
            iswitch=1                                                   6d17s22
           else                                                         6d17s22
            nrow=(noc(ksb)*(noc(ksb)+1))/2                              6d17s22
            iswitch=0                                                   6d17s22
           end if                                                       6d17s22
           i10=i1s                                                      6d17s22
           i1n=noc(isblkp(3,is))                                        6d17s22
           ioooo=i4o(is)                                                6d17s22
           do i2=i2s,i2e                                                6d17s22
            if(i2.eq.i2e)i1n=i1e                                        6d17s22
            i2m=i2-1                                                    6d17s22
            do i1=i10,i1n                                               6d17s22
             if(i1.gt.idoub(isblkp(3,is)).and.
     $            i2.le.idoub(isblkp(4,is)))then                        6d17s22
              i1m=i1-1-idoub(ksb)                                       6d17s22
              iad3=iden1(isb)+iacto(isb)*i1m                            6d17s22
              do id=0,idoub(ksb)-1                                      6d17s22
               ix=max(id,i2m)                                           6d17s22
               in=min(id,i2m)                                           6d17s22
               itri=((ix*(ix+1))/2)+in                                  6d17s22
               iad1=ioooo+itri+iswitch*(id+idoub(ksb)*i2m-itri)         6d17s22
               fact=2d0*bc(iad1)                                        6d17s22
               iad2=iamat(isb)+iacto(isb)*id                            6d17s22
               do ia=0,iacto(isb)-1                                     6d17s22
                bc(iad2+ia)=bc(iad2+ia)+fact*bc(iad3+ia)                6d17s22
               end do                                                   6d17s22
              end do                                                    6d17s22
             end if                                                     6d17s22
             ioooo=ioooo+nrow                                           6d17s22
            end do                                                      6d17s22
            i10=1                                                       6d17s22
           end do                                                       6d17s22
          end if                                                        6d17s22
c     Aad=-2 sum (aa'|da") Da'a" : o
c      ik         i   k
          if(isblkp(1,is).eq.isb.and.isblkp(3,is).eq.ksb)then           6d17s22
           call ilimts(noc(isblkp(3,is)),noc(isblkp(4,is)),mynprocg,    6d17s22
     $         mynowprog,il,ih,i1s,i1e,i2s,i2e)                         6d17s22
           if(isblkp(1,is).eq.isblkp(2,is))then                         6d17s22
            nrow=(noc(isb)*(noc(isb)+1))/2                              6d17s22
            iswitch=0                                                   6d17s22
           else
            nrow=noc(isb)*noc(isblkp(2,is))                             6d17s22
            iswitch=1                                                   6d17s22
           end if                                                       6d17s22
           i10=i1s                                                      6d17s22
           i1n=noc(ksb)                                                 6d17s22
           ioooo=i4o(is)                                                6d17s22
           do i2=i2s,i2e                                                6d17s22
            if(i2.eq.i2e)i1n=i1e                                        6d17s22
            i2m=i2-1-idoub(isblkp(4,is))                                6d17s22
            jden=iden1(isblkp(2,is))+iacto(isblkp(2,is))*i2m            6d17s22
            do i1=i10,i1n                                               6d17s22
             if(i1.le.idoub(ksb).and.i2.gt.idoub(isblkp(4,is)))then     6d17s22
              iad1=iamat(isb)+iacto(isb)*(i1-1)                         6d17s22
              do i4=0,iacto(isblkp(2,is))-1                             6d17s22
               fact=-2d0*bc(jden+i4)                                    6d17s22
               i4p=i4+idoub(isblkp(2,is))                               6d17s22
               do i3=0,iacto(isb)-1                                     6d17s22
                i3p=i3+idoub(isb)                                       6d17s22
                ix=max(i3p,i4p)                                         6d17s22
                in=min(i3p,i4p)                                         6d17s22
                itri=((ix*(ix+1))/2)+in                                 6d17s22
                iad2=ioooo+itri+iswitch*(i3p+noc(isb)*i4p-itri)         6d17s22
                bc(iad1+i3)=bc(iad1+i3)+fact*bc(iad2)
               end do                                                   6d17s22
              end do                                                    6d17s22
             end if                                                     6d17s22
             ioooo=ioooo+nrow                                           6d17s22
            end do                                                      6d17s22
            i10=1                                                       6d17s22
           end do                                                       6d17s22
          else if(isblkp(1,is).eq.isb.and.isblkp(4,is).eq.ksb)then      6d17s22
           call ilimts(noc(isblkp(3,is)),noc(isblkp(4,is)),mynprocg,    6d17s22
     $         mynowprog,il,ih,i1s,i1e,i2s,i2e)                         6d17s22
           if(isblkp(1,is).eq.isblkp(2,is))then                         6d17s22
            nrow=(noc(isb)*(noc(isb)+1))/2                              6d17s22
            iswitch=0                                                   6d17s22
           else
            nrow=noc(isb)*noc(isblkp(2,is))                             6d17s22
            iswitch=1                                                   6d17s22
           end if                                                       6d17s22
           i10=i1s                                                      6d17s22
           i1n=noc(isblkp(3,is))                                        6d17s22
           ioooo=i4o(is)                                                6d17s22
           do i2=i2s,i2e                                                6d17s22
            if(i2.eq.i2e)i1n=i1e                                        6d17s22
            i2m=i2-1
            iad1=iamat(isb)+iacto(isb)*i2m                              6d17s22
            do i1=i10,i1n                                               6d17s22
             if(i2.le.idoub(ksb).and.i1.gt.idoub(isblkp(3,is)))then     6d17s22
              i1m=i1-1-idoub(isblkp(3,is))                              6d17s22
              jden=iden1(isblkp(2,is))+iacto(isblkp(2,is))*i1m            6d17s22
              do i4=0,iacto(isblkp(2,is))-1                             6d17s22
               fact=-2d0*bc(jden+i4)                                    6d17s22
               i4p=i4+idoub(isblkp(2,is))                               6d17s22
               do i3=0,iacto(isb)-1                                     6d17s22
                i3p=i3+idoub(isb)                                       6d17s22
                ix=max(i3p,i4p)                                         6d17s22
                in=min(i3p,i4p)                                         6d17s22
                itri=((ix*(ix+1))/2)+in                                 6d17s22
                iad2=ioooo+itri+iswitch*(i3p+noc(isb)*i4p-itri)         6d17s22
                bc(iad1+i3)=bc(iad1+i3)+fact*bc(iad2)
               end do                                                   6d17s22
              end do                                                    6d17s22
             end if                                                     6d17s22
             ioooo=ioooo+nrow                                           6d17s22
            end do                                                      6d17s22
            i10=1                                                       6d17s22
           end do                                                       6d17s22
          end if                                                        6d17s22
         end if                                                         6d16s22
        end do                                                          12d16s17
        jamat=iamat(isb)+idoub(multh(iduse,jsb))*iacto(isb)                          6d7s22
        do is=1,nsdlk1p                                                 12d20s17
         if(iduse.eq.1)then                                             6d16s22
          if(isblk1p(3,is).eq.jsb.and.isblk1p(4,is).eq.isb)then          6d8s22
           is12=isblk1p(1,is)                                            12d20s17
           call ilimts(noc(jsb),nvirtc(isb),mynprocg,mynowprog,il,ih,    6d8s22
     $         i1s,i1e,i2s,i2e)                                         12d20s17
           i10=i1s                                                       12d20s17
           i1n=noc(jsb)                                                  6d8s22
           if(ipuse.eq.1)then                                            6d7s22
            nrow=(noc(is12)*(noc(is12)+1))/2                              12d20s17
            isw=0                                                        6d7s22
           else                                                          6d7s22
            nrow=noc(is12)*noc(is12)                                     6d7s22
            isw=1                                                        6d7s22
           end if                                                        6d7s22
           i1x=ionex(is)                                                 12d20s17
           do i2=i2s,i2e                                                 12d20s17
            if(i2.eq.i2e)i1n=i1e                                         12d20s17
            i2m=i2-1                                                     12d20s17
            do i1=i10,i1n                                                12d20s17
             if(i1.gt.idoub(jsb))then                                    6d8s22
              i1m=i1-idoub(jsb)-1                                        6d8s22
              jden=iden1(jsb)+iacto(jsb)*i1m                             6d8s22
              iad3=jamat+i2m+nvirtc(isb)*idoub(jsb)                      6d8s22
              do i34=0,idoub(is12)-1                                     12d20s17
               irec=i34*(noc(is12)+1)                                    6d7s22
               itri=((i34*(i34+1))/2)+i34                                6d7s22
               itri=itri+isw*(irec-itri)                                 6d7s22
               iad1=i1x+itri                                             6d7s22
               fact=+4d0*bc(iad1)                                        12d20s17
               do ia=0,iacto(jsb)-1                                      6d8s22
                jad3=iad3+nvirtc(isb)*ia                                 12d20s17
                bc(jad3)=bc(jad3)+fact*bc(jden+ia)                       12d20s17
c     Ava= 4 sum Daa' (dd|a'v): v
c      ik         ki      i i
               end do                                                    12d20s17
              end do                                                     12d20s17
             else                                                        12d20s17
              i1m=i1-1                                                   12d20s17
              iad3=jamat+i2m+nvirtc(isb)*i1m                             12d20s17
              do i4=0,iacto(is12)-1                                      12d20s17
               i4p=i4+idoub(is12)                                        12d20s17
               jden=iden1(is12)+iacto(is12)*i4                           12d20s17
               do i3=0,iacto(is12)-1                                     12d20s17
                i3p=i3+idoub(is12)                                       12d20s17
                irec=i3p+noc(is12)*i4p                                  6d18s22
                ix=max(i3p,i4p)                                          12d20s17
                in=min(i3p,i4p)                                          12d20s17
                itri=((ix*(ix+1))/2)+in                                  6d7s22
                itri=itri+isw*(irec-itri)                                6d7s22
                iad1=i1x+itri                                            6d7s22
                bc(iad3)=bc(iad3)+4d0*bc(iad1)*bc(jden+i3)               12d20s17
c     Avd= 4 sum Da'a" (a'a"|dv): w
c      ik                    ki
               end do                                                    12d20s17
              end do                                                     12d20s17
             end if                                                      12d20s17
             i1x=i1x+nrow                                                12d20s17
            end do                                                       12d20s17
            i10=1                                                        12d20s17
           end do                                                        12d20s17
          end if                                                         12d20s17
          if(isblk1p(1,is).eq.jsb.and.isblk1p(4,is).eq.isb)then          6d8s22
           is23=isblk1p(2,is)                                            12d20s17
           call ilimts(noc(is23),nvirtc(isb),mynprocg,mynowprog,il,ih,   12d20s17
     $         i1s,i1e,i2s,i2e)                                         12d20s17
           i10=i1s                                                       12d20s17
           i1n=noc(is23)                                                 12d20s17
           if(is23.eq.jsb.and.ipuse.eq.1)then                            6d9s22
            nrow=(noc(jsb)*(noc(jsb)+1))/2                               6d9s22
            iswitch=0                                                    12d20s17
           else                                                          12d20s17
            nrow=noc(jsb)*noc(is23)                                      6d9s22
            iswitch=1                                                    12d20s17
           end if                                                        12d20s17
           i1x=ionex(is)                                                 12d21s17
           do i2=i2s,i2e                                                 12d20s17
            if(i2.eq.i2e)i1n=i1e                                         12d20s17
            i2m=i2-1                                                     12d20s17
            do i1=i10,i1n                                                12d20s17
             if(i1.le.idoub(is23))then                                   12d20s17
              i1m=i1-1                                                   12d20s17
              do ia=0,iacto(jsb)-1                                       12d20s17
               jden=iden1(jsb)+iacto(jsb)*ia                             6d7s22
               iad3=jamat+i2m+nvirtc(isb)*(idoub(jsb)+ia)                6d8s22
               do i3=0,iacto(jsb)-1                                      12d20s17
                i3p=i3+idoub(jsb)                                        6d8s22
                ix=max(i3p,i1m)                                          12d20s17
                in=min(i3p,i1m)                                          12d20s17
                ieq=((ix*(ix+1))/2)+in                                   12d20s17
                inot=i3p+noc(jsb)*i1m                                    12d20s17
                iad1=i1x+(inot-ieq)*iswitch+ieq                          12d20s17
                bc(iad3)=bc(iad3)-2d0*bc(iad1)*bc(jden+i3)               12d20s17
c     Ava= -2 sum Da'a (a'd|dv): x
               end do                                                    12d20s17
              end do                                                     12d20s17
             else                                                        12d20s17
              i1m=i1-1-idoub(is23)
              do i4=0,iacto(is23)-1                                      12d20s17
               i4p=i4+idoub(is23)                                        12d20s17
               iad2=iden1(is23)+i4+iacto(is23)*i1m                       12d20s17
               do i3=0,idoub(jsb)-1                                      6d8s22
                ix=max(i3,i4p)                                           12d20s17
                in=min(i3,i4p)                                           12d20s17
                ieq=((ix*(ix+1))/2)+in                                   12d20s17
                inot=i3+noc(jsb)*i4p                                     6d8s22
                iad1=i1x+(inot-ieq)*iswitch+ieq                          12d20s17
                iad3=jamat+i2m+nvirtc(isb)*i3                            12d20s17
                bc(iad3)=bc(iad3)-2d0*bc(iad1)*bc(iad2)                  12d20s17
c     Avd= -2 sum Da'a (da'|av): y
               end do                                                    12d20s17
              end do                                                     12d20s17
             end if                                                      12d20s17
             i1x=i1x+nrow                                                12d20s17
            end do                                                       12d20s17
            i10=1                                                        12d20s17
           end do                                                        12d20s17
          else if(isblk1p(2,is).eq.jsb.and.isblk1p(4,is).eq.isb)then     6d8s22
           is13=isblk1p(1,is)                                            12d20s17
           call ilimts(noc(is13),nvirtc(isb),mynprocg,mynowprog,il,ih,   12d20s17
     $         i1s,i1e,i2s,i2e)                                         12d20s17
           i10=i1s                                                       12d20s17
           i1n=noc(is13)                                                 12d20s17
           nrow=noc(jsb)*noc(is13)                                       6d9s22
           i1x=ionex(is)                                                 12d21s17
           do i2=i2s,i2e                                                 12d20s17
            if(i2.eq.i2e)i1n=i1e                                         12d20s17
            i2m=i2-1                                                     12d20s17
            do i1=i10,i1n                                                12d20s17
             if(i1.le.idoub(is13))then                                   12d20s17
              i1m=i1-1                                                   12d20s17
              do ia=0,iacto(jsb)-1                                       12d20s17
               jden=iden1(jsb)+iacto(jsb)*ia                             12d20s17
               iad3=jamat+i2m+nvirtc(isb)*(idoub(jsb)+ia)                12d20s17
               do i3=0,iacto(jsb)-1                                      12d20s17
                i3p=i3+idoub(jsb)                                        12d20s17
                iad1=i1x+i1m+noc(is13)*i3p                               12d20s17
                bc(iad3)=bc(iad3)-2d0*bc(iad1)*bc(jden+i3)               12d20s17
c     Ava= -2 sum Da'a (da'|dv): z
               end do                                                    12d20s17
              end do                                                     12d20s17
             else                                                        12d20s17
              i1m=i1-1-idoub(is13)
              do i4=0,iacto(is13)-1                                      12d20s17
               i4p=i4+idoub(is13)                                        12d20s17
               iad2=iden1(is13)+i4+iacto(is13)*i1m                       12d20s17
               do i3=0,idoub(jsb)-1                                      6d8s22
                iad1=i1x+i4p+noc(is13)*i3                                12d20s17
                iad3=jamat+i2m+nvirtc(isb)*i3                            12d20s17
c     Avd= -2 sum Daa' (ad|a'v): ba
                bc(iad3)=bc(iad3)-2d0*bc(iad1)*bc(iad2)                  12d20s17
               end do                                                    12d20s17
              end do                                                     12d20s17
             end if                                                      12d20s17
             i1x=i1x+nrow                                                12d20s17
            end do                                                       12d20s17
            i10=1                                                        12d20s17
           end do                                                        12d20s17
          end if                                                         12d20s17
         else                                                           6d16s22
          ksb=multh(isb,iduse)                                          6d17s22
c     Ava= 4 sum Daa' (dd|a'v): v
c      ik         ki      i i
          if(isblk1p(1,is).eq.isblk1p(2,is).and.                        6d18s22
     $         isblk1p(4,is).eq.isb)then                                6d18s22
           call ilimts(noc(isb),nvirt(isb),mynprocg,mynowprog,il,ih,    6d18s22
     $         i1s,i1e,i2s,i2e)                                         6d18s22
           nrow=(noc(isblk1p(1,is))*(noc(isblk1p(1,is))+1))/2           6d18s22
           i1x=ionex(is)                                                6d20s22
           i10=i1s                                                      6d18s22
           i1n=noc(isb)                                                 6d18s22
           do i2=i2s,i2e                                                6d18s22
            if(i2.eq.i2e)i1n=i1e                                        6d18s22
            i2m=i2-1                                                    6d18s22
            do i1=i10,i1n                                               6d18s22
             if(i1.gt.idoub(isb))then                                   6d18s22
              i1m=i1-1-idoub(isb)                                       6d18s22
              iadd=iden1(ksb)+iacto(ksb)*i1m                            6d18s22
              do id=1,idoub(isblk1p(1,is))                              6d18s22
               iadi=i1x+((id*(id+1))/2)-1                               6d19s22
               ff=bc(iadi)*4d0                                          6d18s22
               do ia=0,iacto(ksb)-1                                      6d18s22
                iada=jamat+i2m+(ia+idoub(ksb))*nvirt(isb)               6d19s22
                bc(iada)=bc(iada)+ff*bc(iadd+ia)                        6d18s22
               end do                                                   6d18s22
              end do                                                    6d18s22
             end if                                                     6d18s22
             i1x=i1x+nrow                                               6d18s22
            end do                                                      6d18s22
            i10=1                                                       6d18s22
           end do                                                       6d18s22
          end if                                                        6d18s22
c     Avd= 4 sum Da'a" (a'a"|dv): w
c      ik               k i  ki
c      ik               i k  ki
          if(isblk1p(3,is).eq.ksb.and.isblk1p(4,is).eq.isb)then         6d18s22
           call ilimts(noc(ksb),nvirt(isb),mynprocg,mynowprog,il,ih,    6d18s22
     $         i1s,i1e,i2s,i2e)                                         6d18s22
           nrow=noc(isblk1p(1,is))*noc(isblk1p(2,is))                   6d18s22
           i1x=ionex(is)                                                6d18s22
           i10=i1s                                                      6d18s22
           i1n=noc(ksb)                                                 6d18s22
           do i2=i2s,i2e                                                6d18s22
            if(i2.eq.i2e)i1n=i1e                                        6d18s22
            i2m=i2-1                                                    6d18s22
            do i1=i10,i1n                                               6d18s22
             if(i1.le.idoub(ksb))then                                   6d18s22
              i1m=i1-1                                                  6d18s22
              iada=jamat+i2m+nvirt(isb)*i1m                             6d18s22
              do iapp=0,iacto(isblk1p(2,is))-1                          6d18s22
               iadi=i1x+idoub(isblk1p(1,is))+noc(isblk1p(1,is))         6d18s22
     $              *(iapp+idoub(isblk1p(2,is)))                        6d18s22
               iadd=iden1(isblk1p(1,is))+iacto(isblk1p(1,is))*iapp      6d18s22
               do iap=0,iacto(isblk1p(1,is))-1                          6d18s22
                bc(iada)=bc(iada)+8d0*bc(iadi+iap)*bc(iadd+iap)         6d20s22
               end do                                                   6d18s22
              end do                                                    6d18s22
             end if                                                     6d18s22
             i1x=i1x+nrow                                               6d18s22
            end do                                                      6d18s22
            i10=1                                                       6d18s22
           end do                                                       6d18s22
          end if                                                        6d20s22
c     Avd= -2 sum Da'a (da'|av): y
c      ik          i k  ki  ki
c      ik          k i  kk  ii
          if(isblk1p(1,is).eq.ksb.and.isblk1p(4,is).eq.isb)then         6d18s22
           call ilimts(noc(isblk1p(3,is)),nvirt(isb),mynprocg,mynowprog,6d20s22
     $          il,ih,i1s,i1e,i2s,i2e)                                  6d20s22
           if(ksb.eq.isblk1p(2,is))then                                 6d20s22
            nrow=(noc(ksb)*(noc(ksb)+1))/2                              6d20s22
            iswitch=0                                                   6d20s22
           else                                                         6d20s22
            nrow=noc(ksb)*noc(isblk1p(2,is))                            6d20s22
            iswitch=1                                                   6d20s22
           end if                                                       6d20s22
           i1x=ionex(is)                                                6d18s22
           i10=i1s                                                      6d18s22
           i1n=noc(isblk1p(3,is))                                       6d20s22
           do i2=i2s,i2e                                                6d18s22
            if(i2.eq.i2e)i1n=i1e                                        6d18s22
            i2m=i2-1                                                    6d18s22
            do i1=i10,i1n                                               6d18s22
             if(i1.gt.idoub(isblk1p(3,is)))then                         6d20s22
              i1m=i1-1-idoub(isblk1p(3,is))                             6d20s22
              do iap=0,iacto(isblk1p(2,is))-1                           6d18s22
               iapx=iap+idoub(isblk1p(2,is))                            6d18s22
               iadd=iden1(isblk1p(2,is))+iap+iacto(isblk1p(2,is))*i1m   6d18s22
               ff=-2d0*bc(iadd)                                         6d18s22
               do id=0,idoub(ksb)-1                                     6d18s22
                iada=jamat+i2m+nvirt(isb)*id                            6d18s22
                ix=max(iapx,id)                                         6d20s22
                in=min(iapx,id)                                         6d20s22
                itri=((ix*(ix+1))/2)+in                                 6d20s22
                iadi=i1x+itri+iswitch*(id+noc(ksb)*iapx-itri)           6d20s22
                bc(iada)=bc(iada)+ff*bc(iadi)                           6d20s22
               end do                                                   6d18s22
              end do                                                    6d18s22
             end if                                                     6d18s22
             i1x=i1x+nrow                                               6d18s22
            end do                                                      6d18s22
            i10=1                                                       6d18s22
           end do                                                       6d18s22
c     Avd= -2 sum Daa' (ad|a'v): ba
c      ik          ik    k   i
c      ik          ki    k   i
          else if(isblk1p(2,is).eq.ksb.and.isblk1p(4,is).eq.isb)then    6d20s22
           call ilimts(noc(isblk1p(3,is)),nvirt(isb),mynprocg,mynowprog,6d18s22
     $         il,ih,i1s,i1e,i2s,i2e)                                   6d18s22
           nrow=noc(isblk1p(1,is))*noc(ksb)                             6d18s22
           i1x=ionex(is)                                                6d18s22
           i10=i1s                                                      6d18s22
           i1n=noc(isblk1p(3,is))                                       6d18s22
           do i2=i2s,i2e                                                6d18s22
            if(i2.eq.i2e)i1n=i1e                                        6d18s22
            i2m=i2-1                                                    6d18s22
            do i1=i10,i1n                                               6d18s22
             if(i1.gt.idoub(isblk1p(3,is)))then                         6d18s22
              i1m=i1-1-idoub(isblk1p(3,is))                             6d18s22
              do id=0,idoub(ksb)-1                                      6d18s22
               iada=jamat+i2m+nvirt(isb)*id                             6d18s22
               iadd=iden1(isblk1p(1,is))+iacto(isblk1p(1,is))*i1m       6d18s22
               iadi=i1x+idoub(isblk1p(1,is))+noc(isblk1p(1,is))*id      6d18s22
               do ia=0,iacto(isblk1p(1,is))-1                           6d18s22
                bc(iada)=bc(iada)-2d0*bc(iadd+ia)*bc(iadi+ia)           6d18s22
               end do                                                   6d18s22
              end do                                                    6d18s22
             end if                                                     6d18s22
             i1x=i1x+nrow                                               6d18s22
            end do                                                      6d18s22
            i10=1                                                       6d18s22
           end do                                                       6d18s22
          end if                                                        6d18s22
c     Ava= -2 sum Da'a (a'd|dv): x
c      ik          i k  i    i
          if(isblk1p(1,is).eq.isb.and.isblk1p(4,is).eq.isb)then         6d18s22
           call ilimts(noc(isblk1p(3,is)),nvirt(isb),mynprocg,mynowprog,6d18s22
     $         il,ih,i1s,i1e,i2s,i2e)                                   6d18s22
           if(isblk1p(2,is).eq.isb)then                                 6d18s22
            nrow=(noc(isb)*(noc(isb)+1))/2                              6d18s22
            iswitch=0                                                   6d18s22
           else                                                         6d18s22
            nrow=noc(isblk1p(2,is))*noc(isb)                            6d18s22
            iswitch=1                                                   6d18s22
           end if                                                       6d18s22
           i1x=ionex(is)                                                6d18s22
           i10=i1s                                                      6d18s22
           i1n=noc(isblk1p(3,is))                                       6d18s22
           do i2=i2s,i2e                                                6d18s22
            if(i2.eq.i2e)i1n=i1e                                        6d18s22
            i2m=i2-1                                                    6d18s22
            do i1=i10,i1n                                               6d18s22
             if(i1.le.idoub(isblk1p(3,is)))then                         6d18s22
              i1m=i1-1                                                  6d18s22
              do iap=0,iacto(isb)-1                                     6d18s22
               iapx=iap+idoub(isb)                                      6d18s22
               ix=max(i1m,iapx)                                         6d18s22
               in=min(i1m,iapx)                                         6d18s22
               itri=((ix*(ix+1))/2)+in                                  6d18s22
               iadi=i1x+itri+iswitch*(iapx+noc(isb)*i1m-itri)           6d18s22
               ff=-2d0*bc(iadi)                                         6d18s22
               do ia=0,iacto(ksb)-1                                     6d18s22
                iadd=iden1(isb)+iap+iacto(isb)*ia                       6d18s22
                iada=jamat+i2m+nvirt(isb)*(ia+idoub(ksb))               6d19s22
                bc(iada)=bc(iada)+ff*bc(iadd)                           6d18s22
               end do                                                   6d18s22
              end do                                                    6d18s22
             end if                                                     6d18s22
             i1x=i1x+nrow                                               6d18s22
            end do                                                      6d18s22
            i10=1                                                       6d18s22
           end do                                                       6d18s22
c     Ava= -2 sum Da'a (da'|dv): z
c      ik          i k   i   i
          else if(isblk1p(2,is).eq.isb.and.isblk1p(4,is).eq.isb)then    6d19s22
           call ilimts(noc(isblk1p(3,is)),nvirt(isb),mynprocg,mynowprog,6d18s22
     $         il,ih,i1s,i1e,i2s,i2e)                                   6d18s22
           if(isblk1p(1,is).eq.isb)then                                 6d18s22
            nrow=(noc(isb)*(noc(isb)+1))/2                              6d18s22
            iswitch=0                                                   6d18s22
           else                                                         6d18s22
            nrow=noc(isblk1p(1,is))*noc(isb)                            6d18s22
            iswitch=1                                                   6d18s22
           end if                                                       6d18s22
           i1x=ionex(is)                                                6d18s22
           i10=i1s                                                      6d18s22
           i1n=noc(isblk1p(3,is))                                       6d18s22
           do i2=i2s,i2e                                                6d18s22
            if(i2.eq.i2e)i1n=i1e                                        6d18s22
            i2m=i2-1                                                    6d18s22
            do i1=i10,i1n                                               6d18s22
             if(i1.le.idoub(isblk1p(3,is)))then                         6d18s22
              i1m=i1-1                                                  6d18s22
              do ia=0,iacto(ksb)-1                                      6d18s22
               iada=jamat+i2m+nvirt(isb)*(ia+idoub(ksb))                6d19s22
               iadd=iden1(isb)+iacto(isb)*ia                            6d18s22
               do iap=0,iacto(isb)-1                                    6d18s22
                iapx=iap+idoub(isb)                                     6d18s22
                ix=max(iapx,i1m)                                        6d18s22
                in=min(iapx,i1m)                                        6d18s22
                itri=((ix*(ix+1))/2)+in                                 6d18s22
                iadi=i1x+itri+iswitch*(i1m+noc(isblk1p(1,is))*iapx-itri)6d18s22
                bc(iada)=bc(iada)-2d0*bc(iadd+iap)*bc(iadi)             6d18s22
               end do                                                   6d18s22
              end do                                                    6d18s22
             end if                                                     6d18s22
             i1x=i1x+nrow                                               6d18s22
            end do                                                      6d18s22
            i10=1                                                       6d18s22
           end do                                                       6d18s22
          end if                                                        6d18s22
         end if                                                         6d16s22
        end do                                                          12d20s17
       end do                                                           12d16s17
      end if                                                            12d16s17
      if(idoit(5).ne.0)then                                             12d21s17
c
c     aaaa part
c
       do isb=1,nsymb                                                   12d21s17
        jsb=multh(isb,ipuse)                                            6d17s22
        if(iduse.eq.1)then                                              6d16s22
         do isd=1,nsdlk                                                  12d21s17
          if(j2den(isd).gt.0.and.                                       6d29s22
     $         (isblk(1,isd).eq.isb.or.                                 6d30s22
     $         (ipuse.ne.1.and.isblk(2,isd).eq.isb)))then               6d30s22
           i2hit=0
           if(isblk(1,isd).eq.isblk(2,isd))then                          12d21s17
            nrowd=(iacto(isblk(1,isd))*(iacto(isblk(1,isd))+1))/2        12d21s17
            jswitch=0                                                    12d21s17
            fuse=-4d0                                                    6d9s22
           else                                                          12d21s17
            nrowd=iacto(isblk(1,isd))*iacto(isblk(2,isd))                12d21s17
            jswitch=1                                                    12d21s17
            fuse=-2d0                                                    6d9s22
           end if
           if(isblk(3,isd).eq.isblk(4,isd))then                          6d9s22
            fmul=4d0                                                     6d9s22
           else                                                          6d9s22
            fmul=2d0                                                     6d9s22
           end if                                                        6d9s22
           do is=1,nsdlkp                                                12d21s17
            if(isblkp(1,is).eq.jsb.and.isblkp(2,is).eq.isblk(2,isd).and. 6d8s22     Aad=sum Gaa'a''a
     $        isblkp(3,is).eq.isblk(3,isd).and.                         12d21s17
     $          isblkp(4,is).eq.isblk(4,isd))then                       12d21s17
             i2hit=1
             call ilimts(noc(isblkp(3,is)),noc(isblkp(4,is)),mynprocg,   12d21s17
     $          mynowprog,il,ih,i1s,i1e,i2s,i2e)                        12d21s17
             i10=i1s                                                     12d21s17
             i1n=noc(isblkp(3,is))                                       12d21s17
             if(isblkp(1,is).eq.isblkp(2,is).and.ipuse.eq.1)then         6d7s22
              nrowi=(noc(isblkp(1,is))*(noc(isblkp(1,is))+1))/2          12d21s17
              iswitch=0                                                  12d21s17
             else                                                        12d21s17
              nrowi=noc(isblkp(1,is))*noc(isblkp(2,is))                  12d21s17
              iswitch=1                                                  12d21s17
             end if                                                      12d21s17
             ioooo=i4o(is)                                               12d21s17
             do i2=i2s,i2e                                               12d21s17
              if(i2.eq.i2e)i1n=i1e                                       12d21s17
              i2m=i2-1-idoub(isblkp(4,is))                               12d21s17
              do i1=i10,i1n                                              12d21s17
               if(i1.gt.idoub(isblkp(3,is)).and.                         12d21s17
     $             i2.gt.idoub(isblkp(4,is)))then                       12d21s17
                i1m=i1-1-idoub(isblkp(3,is))                             12d21s17
                ix=max(i2m,i1m)                                          12d21s17
                in=min(i2m,i1m)                                          12d21s17
                ieq=((ix*(ix+1))/2)+in                                   12d21s17
                inot=i1m+iacto(isblk(3,isd))*i2m                         12d21s17
                j2=j2den(isd)+nrowd*((inot-ieq)*jswitch+ieq)             12d21s17
                do i5=0,idoub(jsb)-1                                     6d7s22
                 do i4=0,iacto(isblkp(2,is))-1                            12d21s17
                  i4p=i4+idoub(isblkp(2,is))                             12d21s17
                  ix=max(i4p,i5)                                         12d21s17
                  in=min(i4p,i5)                                         12d21s17
                  ieq=((ix*(ix+1))/2)+in                                 12d21s17
                  inot=i5+noc(jsb)*i4p                                   6d8s22
                  iad1=ioooo+(inot-ieq)*iswitch+ieq                      12d21s17
                  do i3=0,iacto(isb)-1                                   6d9s22
                   ix=max(i3,i4)                                         12d21s17
                   in=min(i3,i4)                                         12d21s17
                   ieq=((ix*(ix+1))/2)+in                                12d21s17
                   inot=i3+iacto(isb)*i4                                 12d21s17
                   iad2=j2+(inot-ieq)*jswitch+ieq                        12d21s17
                   iad3=iamat(isb)+i3+iacto(isb)*i5                      12d21s17
c      35   34 1  2     54  1    2
c     Aad= Gaa'a''a''' (da'|a''a''')
                   bc(iad3)=bc(iad3)+fuse*bc(iad2)*bc(iad1)              12d21s17
                  end do                                                 12d21s17
                 end do                                                  12d21s17
                end do                                                   12d21s17
               end if                                                    12d21s17
               ioooo=ioooo+nrowi                                         12d21s17
              end do                                                     12d21s17
              i10=1                                                      12d21s17
             end do                                                      12d21s17
            else if(isblkp(1,is).eq.jsb.and.                            6d29s22
     $            isblkp(2,is).eq.isblk(1,isd).and.                     6d29s22
     $        isblkp(3,is).eq.isblk(4,isd).and.                         6d29s22
     $          isblkp(4,is).eq.isblk(3,isd))then                       6d29s22
             i2hit=1
             call ilimts(noc(isblkp(3,is)),noc(isblkp(4,is)),mynprocg,   12d21s17
     $          mynowprog,il,ih,i1s,i1e,i2s,i2e)                        12d21s17
             i10=i1s                                                     12d21s17
             i1n=noc(isblkp(3,is))                                       12d21s17
             if(isblkp(1,is).eq.isblkp(2,is).and.ipuse.eq.1)then         6d7s22
              nrowi=(noc(isblkp(1,is))*(noc(isblkp(1,is))+1))/2          12d21s17
              iswitch=0                                                  12d21s17
             else                                                        12d21s17
              nrowi=noc(isblkp(1,is))*noc(isblkp(2,is))                  12d21s17
              iswitch=1                                                  12d21s17
             end if                                                      12d21s17
             ioooo=i4o(is)                                               12d21s17
             do i2=i2s,i2e                                               12d21s17
              if(i2.eq.i2e)i1n=i1e                                       12d21s17
              i2m=i2-1-idoub(isblkp(4,is))                               12d21s17
              do i1=i10,i1n                                              12d21s17
               if(i1.gt.idoub(isblkp(3,is)).and.                         12d21s17
     $             i2.gt.idoub(isblkp(4,is)))then                       12d21s17
                i1m=i1-1-idoub(isblkp(3,is))                             12d21s17
                inotr=i2m+iacto(isblkp(4,is))*i1m                       6d30s22
                j2=j2den(isd)+nrowd*inotr                                6d29s22
                do i5=0,idoub(jsb)-1                                     6d7s22
                 do i4=0,iacto(isblkp(2,is))-1                            12d21s17
                  i4p=i4+idoub(isblkp(2,is))                             12d21s17
                  ix=max(i4p,i5)                                         12d21s17
                  in=min(i4p,i5)                                         12d21s17
                  ieq=((ix*(ix+1))/2)+in                                 12d21s17
                  inot=i5+noc(jsb)*i4p                                   6d8s22
                  iad1=ioooo+(inot-ieq)*iswitch+ieq                      12d21s17
                  do i3=0,iacto(isb)-1                                   6d9s22
                   inot=i4+iacto(isblkp(2,is))*i3                       6d29s22
                   iad2=j2+inot                                         6d29s22
                   iad3=iamat(isb)+i3+iacto(isb)*i5                      12d21s17
c      35   43 2  1     54  1    2
c     Aad= Ga'aa'''a'' (da'|a''a''')
                   bc(iad3)=bc(iad3)+fuse*bc(iad2)*bc(iad1)              12d21s17
                  end do                                                 12d21s17
                 end do                                                  12d21s17
                end do                                                   12d21s17
               end if                                                    12d21s17
               ioooo=ioooo+nrowi                                         12d21s17
              end do                                                     12d21s17
              i10=1                                                      12d21s17
             end do                                                      12d21s17
            else if(isblkp(2,is).eq.jsb.and.isblkp(1,is).eq.isblk(1,isd)8d4s22
     $           .and.isblkp(3,is).eq.isblk(3,isd).and.                 6d8s22
     $          isblkp(4,is).eq.isblk(4,isd))then                       12d21s17
             i2hit=1
             call ilimts(noc(isblkp(3,is)),noc(isblkp(4,is)),mynprocg,   12d21s17
     $          mynowprog,il,ih,i1s,i1e,i2s,i2e)                        12d21s17
             i10=i1s                                                     12d21s17
             i1n=noc(isblkp(3,is))                                       12d21s17
             if(isblkp(2,is).eq.isblkp(1,is).and.ipuse.eq.1)then         6d7s22
              nrowi=(noc(isblkp(2,is))*(noc(isblkp(2,is))+1))/2          12d21s17
              iswitch=0                                                  12d21s17
             else                                                        12d21s17
              nrowi=noc(isblkp(2,is))*noc(isblkp(1,is))                  12d21s17
              iswitch=1                                                  12d21s17
             end if                                                      12d21s17
             ioooo=i4o(is)                                               12d21s17
             do i2=i2s,i2e                                               12d21s17
              if(i2.eq.i2e)i1n=i1e                                       12d21s17
              i2m=i2-1-idoub(isblkp(4,is))                               12d21s17
              do i1=i10,i1n                                              12d21s17
               if(i1.gt.idoub(isblkp(3,is)).and.                         12d21s17
     $             i2.gt.idoub(isblkp(4,is)))then                       12d21s17
                i1m=i1-1-idoub(isblkp(3,is))                             12d21s17
                ix=max(i2m,i1m)                                          12d21s17
                in=min(i2m,i1m)                                          12d21s17
                ieq=((ix*(ix+1))/2)+in                                   12d21s17
                inot=i1m+iacto(isblk(3,isd))*i2m                         12d21s17
                j2=j2den(isd)+nrowd*((inot-ieq)*jswitch+ieq)             12d21s17
                do i5=0,idoub(jsb)-1                                     6d7s22
                 do i4=0,iacto(isblkp(1,is))-1                            12d21s17
                  i4p=i4+idoub(isblkp(1,is))                             12d21s17
                  ix=max(i4p,i5)                                         12d21s17
                  in=min(i4p,i5)                                         12d21s17
                  ieq=((ix*(ix+1))/2)+in                                 12d21s17
                  inot=i4p+noc(isblkp(1,is))*i5                         6d29s22
                  iad1=ioooo+(inot-ieq)*iswitch+ieq                      12d21s17
                  do i3=0,iacto(isb)-1                                   6d9s22
                   ix=max(i3,i4)                                         12d21s17
                   in=min(i3,i4)                                         12d21s17
                   ieq=((ix*(ix+1))/2)+in                                12d21s17
                   inot=i4+iacto(isblkp(1,is))*i3                       7d8s22
                   iad2=j2+(inot-ieq)*jswitch+ieq                        12d21s17
                   iad3=iamat(isb)+i3+iacto(isb)*i5                      12d21s17
c      35   43 1  2     4 5 1    2
c     Aad= Ga'aa''a''' (a'd|a''a''')
                   bc(iad3)=bc(iad3)+fuse*bc(iad2)*bc(iad1)              12d21s17
                  end do                                                 12d21s17
                 end do                                                  12d21s17
                end do                                                   12d21s17
               end if                                                    12d21s17
               ioooo=ioooo+nrowi                                         12d21s17
              end do                                                     12d21s17
              i10=1                                                      12d21s17
             end do                                                      12d21s17
            end if                                                       12d21s17
           end do                                                        12d21s17
           jamat=iamat(isb)+iacto(isb)*idoub(jsb)                        6d7s22
           jamatj=iamat(jsb)+iacto(jsb)*idoub(isb)                       6d9s22
           do is=1,nsdlk1p                                               12d22s17
            if(isblk1p(4,is).eq.jsb.and.isblk1p(3,is).eq.isblk(2,isd).   6d9s22    Avd=sum Gaa'a''a'
     $          and.isblk1p(1,is).eq.isblk(3,isd).and.isblk1p(2,is).eq. 12d22s17   ji      i2 3  4
     $          isblk(4,isd))then                                       12d22s17           21 2  1
c     Ava=fmul sum Gaa'''a'a''(a'a''|a'''v): bc
c      ji           i                    j
             call ilimts(noc(isblk1p(3,is)),nvirtc(jsb),mynprocg,        6d9s22
     $          mynowprog,il,ih,i1s,i1e,i2s,i2e)                        12d22s17
             i10=i1s                                                     12d22s17
             i1n=noc(isblk1p(3,is))                                      12d22s17
             if(isblk1p(1,is).eq.isblk1p(2,is).and.ipuse.eq.1)then       6d9s22
              nrowi=(noc(isblk1p(1,is))*(noc(isblk1p(1,is))+1))/2        12d22s17
              iswitch=0                                                  12d22s17
             else                                                        12d22s17
              nrowi=noc(isblk1p(1,is))*noc(isblk1p(2,is))                12d22s17
              iswitch=1                                                  12d22s17
             end if                                                      12d22s17
             i1x=ionex(is)                                               12d22s17
             do i2=i2s,i2e                                               12d22s17
              if(i2.eq.i2e)i1n=i1e                                       12d22s17
              i2m=i2-1                                                   12d22s17
c     a'''
              do i1=i10,i1n                                              12d22s17
               if(i1.gt.idoub(isblk1p(3,is)))then                        12d22s17
                i1m=i1-1-idoub(isblk1p(3,is))                            12d22s17
c     a''
                do i4=0,iacto(isblk1p(2,is))-1                           12d22s17
                 i4p=i4+idoub(isblk1p(2,is))                             12d22s17
c     a'
                 do i3=0,iacto(isblk1p(1,is))-1                          12d22s17
                  i3p=i3+idoub(isblk1p(1,is))                            12d22s17
                  ix=max(i3p,i4p)                                        12d22s17
                  in=min(i3p,i4p)                                        12d22s17
                  ieq=((ix*(ix+1))/2)+in                                 12d22s17
                  inot=i3p+noc(isblk1p(1,is))*i4p                        12d22s17
                  iad1=i1x+(inot-ieq)*iswitch+ieq                        12d22s17
                  ix=max(i3,i4)                                          12d22s17
                  in=min(i3,i4)                                          12d22s17
                  ieq=((ix*(ix+1))/2)+in                                 12d22s17
                  inot=i3+iacto(isblk(3,isd))*i4                         12d22s17
                  j2=j2den(isd)+nrowd*((inot-ieq)*jswitch+ieq)           12d22s17
c     a
                  do i5=0,iacto(isb)-1                                   12d22s17
                   iad3=jamatj+i2m+nvirtc(jsb)*(i5+idoub(isb))           6d9s22
                   ix=max(i5,i1m)                                        12d22s17
                   in=min(i5,i1m)                                        12d22s17
                   ieq=((ix*(ix+1))/2)+in                                12d22s17
                   inot=i5+iacto(isb)*i1m                                6d9s22
                   iad2=j2+(inot-ieq)*jswitch+ieq                        12d22s17
c     Ava=fmul sum Gaa'''a'a''(a'a''|a'''v): bc
c      ji           i                    j
                   bc(iad3)=bc(iad3)+fmul*bc(iad2)*bc(iad1)              12d22s17
                  end do                                                 12d22s17
                 end do                                                  12d22s17
                end do                                                   12d22s17
               end if                                                    12d22s17
               i1x=i1x+nrowi                                             12d22s17
              end do                                                     12d22s17
              i10=1                                                      12d22s17
             end do                                                      12d22s17
            else if(ipuse.ne.1.and.isblk1p(4,is).eq.jsb.and.            6d30s22
     $            isblk1p(3,is).eq.isblk(1,isd).                        7d1s22
     $          and.isblk1p(1,is).eq.isblk(3,isd).and.isblk1p(2,is).eq. 6d30s22
     $          isblk(4,isd).and.isblk(1,isd).ne.isblk(2,isd))then      6d30s22
c     Ava=fmul sum Ga'''aa'a''(a'a''|a'''v): bc
c      ji               i                j
             call ilimts(noc(isblk1p(3,is)),nvirtc(jsb),mynprocg,        6d9s22
     $          mynowprog,il,ih,i1s,i1e,i2s,i2e)                        12d22s17
             i10=i1s                                                     12d22s17
             i1n=noc(isblk1p(3,is))                                      12d22s17
             if(isblk1p(1,is).eq.isblk1p(2,is).and.ipuse.eq.1)then       6d9s22
              nrowi=(noc(isblk1p(1,is))*(noc(isblk1p(1,is))+1))/2        12d22s17
              iswitch=0                                                  12d22s17
             else                                                        12d22s17
              nrowi=noc(isblk1p(1,is))*noc(isblk1p(2,is))                12d22s17
              iswitch=1                                                  12d22s17
             end if                                                      12d22s17
             i1x=ionex(is)                                               12d22s17
             do i2=i2s,i2e                                               12d22s17
              if(i2.eq.i2e)i1n=i1e                                       12d22s17
              i2m=i2-1                                                   12d22s17
c     a'''
              do i1=i10,i1n                                              12d22s17
               if(i1.gt.idoub(isblk1p(3,is)))then                        12d22s17
                i1m=i1-1-idoub(isblk1p(3,is))                            12d22s17
c     a''
                do i4=0,iacto(isblk1p(2,is))-1                           12d22s17
                 i4p=i4+idoub(isblk1p(2,is))                             12d22s17
c     a'
                 do i3=0,iacto(isblk1p(1,is))-1                          12d22s17
                  i3p=i3+idoub(isblk1p(1,is))                            12d22s17
                  ix=max(i3p,i4p)                                        12d22s17
                  in=min(i3p,i4p)                                        12d22s17
                  ieq=((ix*(ix+1))/2)+in                                 12d22s17
                  inot=i3p+noc(isblk1p(1,is))*i4p                        12d22s17
                  iad1=i1x+(inot-ieq)*iswitch+ieq                        12d22s17
                  ix=max(i3,i4)                                          12d22s17
                  in=min(i3,i4)                                          12d22s17
                  ieq=((ix*(ix+1))/2)+in                                 12d22s17
                  inot=i3+iacto(isblk(3,isd))*i4                         12d22s17
                  j2=j2den(isd)+nrowd*((inot-ieq)*jswitch+ieq)           12d22s17
c     a
                  do i5=0,iacto(isb)-1                                   12d22s17
                   iad3=jamatj+i2m+nvirtc(jsb)*(i5+idoub(isb))           6d9s22
                   ix=max(i5,i1m)                                        12d22s17
                   in=min(i5,i1m)                                        12d22s17
                   ieq=((ix*(ix+1))/2)+in                                12d22s17
                   inot=i1m+iacto(isblk1p(3,is))*i5                     7d1s22
                   iad2=j2+(inot-ieq)*jswitch+ieq                        12d22s17
c     Ava=fmul sum Ga'''aa'a''(a'a''|a'''v): bc
c      ji               i                j
                   bc(iad3)=bc(iad3)+fmul*bc(iad2)*bc(iad1)              12d22s17
                  end do                                                 12d22s17
                 end do                                                  12d22s17
                end do                                                   12d22s17
               end if                                                    12d22s17
               i1x=i1x+nrowi                                             12d22s17
              end do                                                     12d22s17
              i10=1                                                      12d22s17
             end do                                                      12d22s17
            end if                                                       12d22s17
           end do                                                        12d22s17
          else if(j2den(isd).gt.0.and.isblk(2,isd).eq.isb)then           12d21s17
           nrowd=iacto(isblk(1,isd))*iacto(isblk(2,isd))                 12d21s17
           do is=1,nsdlkp                                                12d21s17
            if(isblkp(1,is).eq.isblk(1,isd).and.isblkp(2,is).eq.jsb.and. 6d8s22
     $        isblkp(3,is).eq.isblk(3,isd).and.                         12d21s17
     $          isblkp(4,is).eq.isblk(4,isd))then                       12d21s17
             call ilimts(noc(isblkp(3,is)),noc(isblkp(4,is)),mynprocg,   12d21s17
     $          mynowprog,il,ih,i1s,i1e,i2s,i2e)                        12d21s17
             i10=i1s                                                     12d21s17
             i1n=noc(isblkp(3,is))                                       12d21s17
             nrowi=noc(isblkp(1,is))*noc(isblkp(2,is))                   12d21s17
             ioooo=i4o(is)                                               12d21s17
             do i2=i2s,i2e                                               12d21s17
              if(i2.eq.i2e)i1n=i1e                                       12d21s17
              i2m=i2-1-idoub(isblkp(4,is))                               12d21s17
              do i1=i10,i1n                                              12d21s17
               if(i1.gt.idoub(isblkp(3,is)).and.                         12d21s17
     $             i2.gt.idoub(isblkp(4,is)))then                       12d21s17
                i1m=i1-1-idoub(isblkp(3,is))                             12d21s17
                j2=j2den(isd)+nrowd*(i1m+iacto(isblk(3,isd))*i2m)        12d21s17
                do i5=0,idoub(jsb)-1                                     6d7s22
                 do i4=0,iacto(isblkp(1,is))-1                            12d21s17
                  i4p=i4+idoub(isblkp(1,is))                             12d21s17
                  iad1=ioooo+i4p+noc(isblkp(1,is))*i5                    12d21s17
                  do i3=0,iacto(isb)-1                                   12d21s17
                   iad2=j2+i4+iacto(isblkp(1,is))*i3                     12d21s17
                   iad3=iamat(isb)+i3+iacto(isb)*i5                      12d21s17
                   bc(iad3)=bc(iad3)-2d0*bc(iad2)*bc(iad1)               12d21s17
                  end do                                                 12d21s17
                 end do                                                  12d21s17
                end do                                                   12d21s17
               end if                                                    12d21s17
               ioooo=ioooo+nrowi                                         12d21s17
              end do                                                     12d21s17
              i10=1                                                      12d21s17
             end do                                                      12d21s17
            end if                                                       12d21s17
           end do                                                        12d21s17
           jamat=iamat(isb)+iacto(isb)*idoub(jsb)                        6d7s22
           jamatj=iamat(jsb)+iacto(jsb)*idoub(isb)                       6d9s22
           do is=1,nsdlk1p                                               12d22s17
            if(isblk1p(4,is).eq.jsb.and.isblk1p(3,is).eq.isblk(1,isd).   6d9s22
     $          and.isblk1p(1,is).eq.isblk(3,isd).and.isblk1p(2,is).eq. 12d22s17   ji      1i 3  4
     $          isblk(4,isd))then                                       12d22s17
             call ilimts(noc(isblk1p(3,is)),nvirtc(jsb),mynprocg,        6d9s22
     $          mynowprog,il,ih,i1s,i1e,i2s,i2e)                        12d22s17
             i10=i1s                                                     12d22s17
             i1n=noc(isblk1p(3,is))                                      12d22s17
             nrowi=noc(isblk1p(1,is))*noc(isblk1p(2,is))                 12d22s17
             i1x=ionex(is)                                               12d22s17
             do i2=i2s,i2e                                               12d22s17
              if(i2.eq.i2e)i1n=i1e                                       12d22s17
              i2m=i2-1                                                   12d22s17
              do i1=i10,i1n                                              12d22s17
               if(i1.gt.idoub(isblk1p(3,is)))then                        12d22s17
                i1m=i1-1-idoub(isblk1p(3,is))                            12d22s17
                do i4=0,iacto(isblk1p(2,is))-1                           12d22s17
                 i4p=i4+idoub(isblk1p(2,is))                             12d22s17
                 do i3=0,iacto(isblk1p(1,is))-1                          12d22s17
                  i3p=i3+idoub(isblk1p(1,is))                            12d22s17
                  iad1=i1x+i3p+noc(isblk1p(1,is))*i4p                    12d22s17
                  j2=j2den(isd)+nrowd*(i3+iacto(isblk(3,isd))*i4)        12d22s17
                  do i5=0,iacto(isb)-1                                   6d9s22
                   iad3=jamatj+i2m+nvirtc(jsb)*(i5+idoub(isb))           6d9s22
                   iad2=j2+i1m+iacto(isblk(1,isd))*i5                    12d22s17
                   bc(iad3)=bc(iad3)+2d0*bc(iad2)*bc(iad1)               12d22s17
                  end do                                                 12d22s17
                 end do                                                  12d22s17
                end do                                                   12d22s17
               end if                                                    12d22s17
               i1x=i1x+nrowi                                             12d22s17
              end do                                                     12d22s17
              i10=1                                                      12d22s17
             end do                                                      12d22s17
            else if(isblk1p(4,is).eq.jsb.and.                            6d9s22
     $           isblk1p(3,is).eq.isblk(1,isd).                         12d22s17
     $          and.isblk1p(2,is).eq.isblk(3,isd).and.isblk1p(1,is).eq. 12d22s17
     $          isblk(4,isd))then                                       12d22s17
             call ilimts(noc(isblk1p(3,is)),nvirtc(jsb),mynprocg,        6d9s22
     $          mynowprog,il,ih,i1s,i1e,i2s,i2e)                        12d22s17
             i10=i1s                                                     12d22s17
             i1n=noc(isblk1p(3,is))                                      12d22s17
             nrowi=noc(isblk1p(1,is))*noc(isblk1p(2,is))                 12d22s17
             i1x=ionex(is)                                               12d22s17
             do i2=i2s,i2e                                               12d22s17
              if(i2.eq.i2e)i1n=i1e                                       12d22s17
              i2m=i2-1                                                   12d22s17
              do i1=i10,i1n                                              12d22s17
               if(i1.gt.idoub(isblk1p(3,is)))then                        12d22s17
                i1m=i1-1-idoub(isblk1p(3,is))                            12d22s17
                do i4=0,iacto(isblk1p(2,is))-1                           12d22s17
                 i4p=i4+idoub(isblk1p(2,is))                             12d22s17
                 do i3=0,iacto(isblk1p(1,is))-1                          12d22s17
                  i3p=i3+idoub(isblk1p(1,is))                            12d22s17
                  iad1=i1x+i3p+noc(isblk1p(1,is))*i4p                    12d22s17
                  j2=j2den(isd)+nrowd*(i3+iacto(isblk(3,isd))*i4)        12d22s17
                  do i5=0,iacto(isb)-1                                   6d9s22
                   iad3=jamat+i2m+nvirtc(jsb)*(i5+idoub(isb))            6d9s22
                   iad2=j2+i5+iacto(isblk(1,isd))*i1m                    12d22s17
                   bc(iad3)=bc(iad3)+2d0*bc(iad2)*bc(iad1)               12d22s17
                  end do                                                 12d22s17
                 end do                                                  12d22s17
                end do                                                   12d22s17
               end if                                                    12d22s17
               i1x=i1x+nrowi                                             12d22s17
              end do                                                     12d22s17
              i10=1                                                      12d22s17
             end do                                                      12d22s17
            end if                                                       12d22s17
           end do                                                        12d22s17
          end if
         end do                                                          12d21s17
        else                                                            6d16s22
         ksb=multh(isb,iduse)                                           6d17s22
         do isd=1,nsblkder                                              6d16s22
          if(j2den(isd).gt.0.and.                                       6d17s22
     $         (isblkder(1,isd).eq.isb.or.isblkder(3,isd).eq.isb.or.    6d19s22
     $         isb.eq.isblkder(2,isd).or.isb.eq.isblkder(4,isd).or.     7d1s22
     $         isblkder(1,isd).eq.ksb.or.isblkder(2,isd).eq.ksb.or.     7d1s22
     $         isblkder(3,isd).eq.ksb.or.isblkder(4,isd).eq.ksb))then   7d1s22
           i2hit=0
           nrowd=iacto(isblkder(1,isd))*iacto(isblkder(2,isd))          6d16s22
           do is=1,nsdlkp                                                12d21s17
            if(isblkp(1,is).eq.isblkp(2,is))then                        6d17s22
             fuse=-1d0                                                  6d17s22
            else                                                        6d17s22
             fuse=-2d0                                                  6d17s22
            end if                                                      6d17s22
c
c     Aad=sum Gaa'a''a''' (da'|a''a''')
c      ik      i           k
            if(isblkp(1,is).eq.ksb.and.isblkp(2,is).eq.isblkder(2,isd)  6d17s22
     $           .and.isblkp(3,is).eq.isblkder(3,isd).and.              6d16s22
     $          isblkp(4,is).eq.isblkder(4,isd))then                    6d16s22
C     Aad=fuse sum Gaa'a''a''' (da'|a''a''')
             i2hit=1                                                    6d16s22
             call ilimts(noc(isblkp(3,is)),noc(isblkp(4,is)),mynprocg,   12d21s17
     $          mynowprog,il,ih,i1s,i1e,i2s,i2e)                        12d21s17
             i10=i1s                                                     12d21s17
             i1n=noc(isblkp(3,is))                                       12d21s17
             if(isblkp(1,is).eq.isblkp(2,is))then                       6d16s22
              nrowi=(noc(isblkp(1,is))*(noc(isblkp(1,is))+1))/2          12d21s17
              iswitch=0                                                  12d21s17
             else                                                        12d21s17
              nrowi=noc(isblkp(1,is))*noc(isblkp(2,is))                  12d21s17
              iswitch=1                                                  12d21s17
             end if                                                      12d21s17
             ioooo=i4o(is)                                               12d21s17
             fuseq=fuse                                                 7d1s22
             if(isblkder(1,isd).ne.isblkder(2,isd).and.                 7d1s22
     $          isblkder(3,isd).ne.isblkder(4,isd))fuseq=fuseq*0.5d0    7d1s22
             do i2=i2s,i2e                                               12d21s17
              if(i2.eq.i2e)i1n=i1e                                       12d21s17
              i2m=i2-1-idoub(isblkp(4,is))                               12d21s17
              do i1=i10,i1n                                              12d21s17
               if(i1.gt.idoub(isblkp(3,is)).and.                         12d21s17
     $             i2.gt.idoub(isblkp(4,is)))then                       12d21s17
                i1m=i1-1-idoub(isblkp(3,is))                             12d21s17
                inot=i1m+iacto(isblkder(3,isd))*i2m                     6d16s22
                j2=j2den(isd)+nrowd*inot                                6d16s22
                do i5=0,idoub(ksb)-1                                    6d17s22
                 do i4=0,iacto(isblkp(2,is))-1                            12d21s17
                  i4p=i4+idoub(isblkp(2,is))                             12d21s17
                  ix=max(i4p,i5)                                         12d21s17
                  in=min(i4p,i5)                                         12d21s17
                  ieq=((ix*(ix+1))/2)+in                                 12d21s17
                  inot=i5+noc(ksb)*i4p                                  6d17s22
                  iad1=ioooo+(inot-ieq)*iswitch+ieq                      12d21s17
                  do i3=0,iacto(isb)-1                                   6d9s22
                   iad2=j2+i3+iacto(isb)*i4                                 12d21s17
                   iad3=iamat(isb)+i3+iacto(isb)*i5                      12d21s17
C     Aad=fuse sum Gaa'a''a''' (da'|a''a''')
                   bc(iad3)=bc(iad3)+fuseq*bc(iad2)*bc(iad1)              6d17s22
                  end do                                                 12d21s17
                 end do                                                  12d21s17
                end do                                                   12d21s17
               end if                                                    12d21s17
               ioooo=ioooo+nrowi                                         12d21s17
              end do                                                     12d21s17
              i10=1                                                      12d21s17
             end do                                                      12d21s17
            else if(isblkp(2,is).eq.ksb.and.                            7d1s22
     $            isblkp(1,is).eq.isblkder(2,isd)                       7d1s22
     $           .and.isblkp(3,is).eq.isblkder(3,isd).and.              6d16s22
     $          isblkp(4,is).eq.isblkder(4,isd))then                    6d16s22
C     Aad=fuse sum Gaa'a''a''' (a'd|a''a''')
             i2hit=1
             call ilimts(noc(isblkp(3,is)),noc(isblkp(4,is)),mynprocg,   12d21s17
     $          mynowprog,il,ih,i1s,i1e,i2s,i2e)                        12d21s17
             i10=i1s                                                     12d21s17
             i1n=noc(isblkp(3,is))                                       12d21s17
             if(isblkp(2,is).eq.isblkp(1,is))then                       6d16s22
              nrowi=(noc(isblkp(2,is))*(noc(isblkp(2,is))+1))/2          12d21s17
              iswitch=0                                                  12d21s17
             else                                                        12d21s17
              nrowi=noc(isblkp(2,is))*noc(isblkp(1,is))                  12d21s17
              iswitch=1                                                  12d21s17
             end if                                                      12d21s17
             ioooo=i4o(is)                                               12d21s17
             fuseq=fuse                                                 7d1s22
             if(isblkder(1,isd).ne.isblkder(2,isd).and.                 7d1s22
     $          isblkder(3,isd).ne.isblkder(4,isd))fuseq=fuseq*0.5d0    7d1s22
             do i2=i2s,i2e                                               12d21s17
              if(i2.eq.i2e)i1n=i1e                                       12d21s17
              i2m=i2-1-idoub(isblkp(4,is))                               12d21s17
              do i1=i10,i1n                                              12d21s17
               if(i1.gt.idoub(isblkp(3,is)).and.                         12d21s17
     $             i2.gt.idoub(isblkp(4,is)))then                       12d21s17
                i1m=i1-1-idoub(isblkp(3,is))                             12d21s17
                inot=i1m+iacto(isblkder(3,isd))*i2m                         12d21s17
                j2=j2den(isd)+nrowd*inot                                6d16s22
                do i5=0,idoub(ksb)-1                                    6d17s22
                 do i4=0,iacto(isblkp(1,is))-1                          7d1s22
                  i4p=i4+idoub(isblkp(1,is))                            7d1s22
                  ix=max(i4p,i5)                                         12d21s17
                  in=min(i4p,i5)                                         12d21s17
                  ieq=((ix*(ix+1))/2)+in                                 12d21s17
                  inot=i4p+noc(isblkp(1,is))*i5                         7d1s22
                  iad1=ioooo+(inot-ieq)*iswitch+ieq                      12d21s17
                  do i3=0,iacto(isb)-1                                   6d9s22
                   iad2=j2+i3+iacto(isb)*i4                             7d1s22
                   iad3=iamat(isb)+i3+iacto(isb)*i5                      12d21s17
C     Aad=fuse sum Gaa'a''a''' (a'd|a''a''')
                   bc(iad3)=bc(iad3)+fuseq*bc(iad2)*bc(iad1)             6d17s22
                  end do                                                 12d21s17
                 end do                                                  12d21s17
                end do                                                   12d21s17
               end if                                                    12d21s17
               ioooo=ioooo+nrowi                                         12d21s17
              end do                                                     12d21s17
              i10=1                                                      12d21s17
             end do                                                      12d21s17
            else if(isblkp(2,is).eq.ksb.and.                            7d1s22
     $            isblkp(1,is).eq.isblkder(1,isd)                       7d1s22
     $           .and.isblkp(3,is).eq.isblkder(3,isd).and.              6d16s22
     $          isblkp(4,is).eq.isblkder(4,isd))then                    6d16s22
C     Aad=fuse sum Ga'aa''a''' (a'd|a''a''') BBB2
c      ik             i           k
             call ilimts(noc(isblkp(3,is)),noc(isblkp(4,is)),mynprocg,   12d21s17
     $          mynowprog,il,ih,i1s,i1e,i2s,i2e)                        12d21s17
             i10=i1s                                                     12d21s17
             i1n=noc(isblkp(3,is))                                       12d21s17
             if(isblkp(2,is).eq.isblkp(1,is))then                       6d16s22
              nrowi=(noc(isblkp(2,is))*(noc(isblkp(2,is))+1))/2          12d21s17
              iswitch=0                                                  12d21s17
             else                                                        12d21s17
              nrowi=noc(isblkp(2,is))*noc(isblkp(1,is))                  12d21s17
              iswitch=1                                                  12d21s17
             end if                                                      12d21s17
             fuseq=fuse                                                 7d1s22
             if(isblkder(1,isd).ne.isblkder(2,isd).and.                 7d1s22
     $          isblkder(3,isd).ne.isblkder(4,isd))fuseq=fuseq*0.5d0    7d1s22
             ioooo=i4o(is)                                               12d21s17
             do i2=i2s,i2e                                               12d21s17
              if(i2.eq.i2e)i1n=i1e                                       12d21s17
              i2m=i2-1-idoub(isblkp(4,is))                               12d21s17
              do i1=i10,i1n                                              12d21s17
               if(i1.gt.idoub(isblkp(3,is)).and.                         12d21s17
     $             i2.gt.idoub(isblkp(4,is)))then                       12d21s17
                i1m=i1-1-idoub(isblkp(3,is))                             12d21s17
                inot=i1m+iacto(isblkder(3,isd))*i2m                         12d21s17
                j2=j2den(isd)+nrowd*inot                                6d16s22
                do i5=0,idoub(ksb)-1                                    6d17s22
                 do i4=0,iacto(isblkp(1,is))-1                          7d1s22
                  i4p=i4+idoub(isblkp(1,is))                            7d1s22
                  ix=max(i4p,i5)                                         12d21s17
                  in=min(i4p,i5)                                         12d21s17
                  ieq=((ix*(ix+1))/2)+in                                 12d21s17
                  inot=i4p+noc(isblkp(1,is))*i5                         7d1s22
                  iad1=ioooo+(inot-ieq)*iswitch+ieq                      12d21s17
                  do i3=0,iacto(isb)-1                                   6d9s22
                   iad2=j2+i4+iacto(isblkp(1,is))*i3                    7d1s22
                   iad3=iamat(isb)+i3+iacto(isb)*i5                      12d21s17
c     Aad=fuse sum Ga'aa''a''' (a'd|a''a''') BBB2
                   bc(iad3)=bc(iad3)+fuseq*bc(iad2)*bc(iad1)             6d17s22
                  end do                                                 12d21s17
                 end do                                                  12d21s17
                end do                                                   12d21s17
               end if                                                    12d21s17
               ioooo=ioooo+nrowi                                         12d21s17
              end do                                                     12d21s17
              i10=1                                                      12d21s17
             end do                                                      12d21s17
            else if(isblkp(1,is).eq.ksb.and.                            7d1s22
     $            isblkp(2,is).eq.isblkder(4,isd)                       7d1s22
     $           .and.isblkp(3,is).eq.isblkder(1,isd).and.              6d16s22
     $          isblkp(4,is).eq.isblkder(2,isd))then                    6d16s22
c     Aad=sum Ga''a'''aa' (da'|a''a''') BB2
c      ik             i    k
             call ilimts(noc(isblkp(3,is)),noc(isblkp(4,is)),mynprocg,   12d21s17
     $          mynowprog,il,ih,i1s,i1e,i2s,i2e)                        12d21s17
             i10=i1s                                                     12d21s17
             i1n=noc(isblkp(3,is))                                       12d21s17
             if(isblkp(2,is).eq.isblkp(1,is))then                       6d16s22
              nrowi=(noc(isblkp(2,is))*(noc(isblkp(2,is))+1))/2          12d21s17
              iswitch=0                                                  12d21s17
             else                                                        12d21s17
              nrowi=noc(isblkp(2,is))*noc(isblkp(1,is))                  12d21s17
              iswitch=1                                                  12d21s17
             end if                                                      12d21s17
             fuseq=fuse                                                 7d1s22
             if(isblkder(1,isd).ne.isblkder(2,isd).and.                 7d1s22
     $          isblkder(3,isd).ne.isblkder(4,isd))fuseq=fuseq*0.5d0    7d1s22
             ioooo=i4o(is)                                               12d21s17
             do i2=i2s,i2e                                               12d21s17
              if(i2.eq.i2e)i1n=i1e                                       12d21s17
              i2m=i2-1-idoub(isblkp(4,is))                               12d21s17
              do i1=i10,i1n                                              12d21s17
               if(i1.gt.idoub(isblkp(3,is)).and.                         12d21s17
     $             i2.gt.idoub(isblkp(4,is)))then                       12d21s17
                i1m=i1-1-idoub(isblkp(3,is))                             12d21s17
                inot=i1m+iacto(isblkder(1,isd))*i2m                     7d1s22
                j2=j2den(isd)+inot                                      7d1s22
                do i5=0,idoub(ksb)-1                                    6d17s22
                 do i4=0,iacto(isblkp(1,is))-1                          7d1s22
                  i4p=i4+idoub(isblkp(1,is))                            7d1s22
                  ix=max(i4p,i5)                                         12d21s17
                  in=min(i4p,i5)                                         12d21s17
                  ieq=((ix*(ix+1))/2)+in                                 12d21s17
                  inot=i4p+noc(isblkp(1,is))*i5                         7d1s22
                  iad1=ioooo+(inot-ieq)*iswitch+ieq                      12d21s17
                  do i3=0,iacto(isb)-1                                   6d9s22
                   iad2=j2+nrowd*(i3+iacto(isb)*i4)                     7d1s22
                   iad3=iamat(isb)+i3+iacto(isb)*i5                      12d21s17
c     Aad=sum Ga''a'''aa' (da'|a''a''') BB2
                   bc(iad3)=bc(iad3)+fuseq*bc(iad2)*bc(iad1)             6d17s22
                  end do                                                 12d21s17
                 end do                                                  12d21s17
                end do                                                   12d21s17
               end if                                                    12d21s17
               ioooo=ioooo+nrowi                                         12d21s17
              end do                                                     12d21s17
              i10=1                                                      12d21s17
             end do                                                      12d21s17
            end if                                                      6d17s22
           end do                                                        12d21s17
           jamat=iamat(isb)+iacto(isb)*idoub(ksb)                        6d7s22
           do is=1,nsdlk1p                                               12d22s17
            if(isblk1p(4,is).eq.isb.and.                                6d16s22
     $          isblkder(1,isd).eq.ksb.and.                             6d19s22
     $           isblk1p(3,is).eq.isblkder(2,isd).                      6d16s22    Avd=sum Gaa'a''a'
     $        and.isblk1p(1,is).eq.isblkder(3,isd).and.isblk1p(2,is).eq. 12d22s17   ji      i2 3  4
     $          isblkder(4,isd))then                                     6d16s22            21 2  1
c     Ava=fmul sum Gaa'''a'a''(a'a''|a'''v): bc
c      ik           k                    i
             call ilimts(noc(isblk1p(3,is)),nvirtc(isb),mynprocg,        6d9s22
     $          mynowprog,il,ih,i1s,i1e,i2s,i2e)                        12d22s17
             i10=i1s                                                     12d22s17
             i1n=noc(isblk1p(3,is))                                      12d22s17
             if(isblk1p(1,is).eq.isblk1p(2,is))then                     6d16s22
              nrowi=(noc(isblk1p(1,is))*(noc(isblk1p(1,is))+1))/2        12d22s17
              iswitch=0                                                  12d22s17
              fmul=1d0                                                  7d1s22
             else                                                        12d22s17
              nrowi=noc(isblk1p(1,is))*noc(isblk1p(2,is))                12d22s17
              iswitch=1                                                  12d22s17
              fmul=2d0                                                  7d1s22
             end if                                                      12d22s17
             if(isblkder(1,isd).ne.isblkder(2,isd).and.                 7d1s22
     $          isblkder(3,isd).ne.isblkder(4,isd))fmul=fmul*0.5d0      7d1s22
             i1x=ionex(is)                                               12d22s17
             do i2=i2s,i2e                                               12d22s17
              if(i2.eq.i2e)i1n=i1e                                       12d22s17
              i2m=i2-1                                                   12d22s17
              do i1=i10,i1n                                              12d22s17
               if(i1.gt.idoub(isblk1p(3,is)))then                        12d22s17
                i1m=i1-1-idoub(isblk1p(3,is))                            12d22s17
                do i4=0,iacto(isblk1p(2,is))-1                           12d22s17
                 i4p=i4+idoub(isblk1p(2,is))                             12d22s17
                 do i3=0,iacto(isblk1p(1,is))-1                          12d22s17
                  i3p=i3+idoub(isblk1p(1,is))                            12d22s17
                  ix=max(i3p,i4p)                                        12d22s17
                  in=min(i3p,i4p)                                        12d22s17
                  ieq=((ix*(ix+1))/2)+in                                 12d22s17
                  inot=i3p+noc(isblk1p(1,is))*i4p                        12d22s17
                  iad1=i1x+(inot-ieq)*iswitch+ieq                        12d22s17
                  ff=bc(iad1)*fmul                                      6d19s22
                  inot=i3+iacto(isblkder(3,isd))*i4                     6d19s22
                  j2=j2den(isd)+nrowd*inot                              6d16s22
                  do i5=0,iacto(ksb)-1                                   12d22s17
                   iad3=jamat+i2m+nvirtc(isb)*(i5+idoub(ksb))           6d19s22
                   iad2=j2+i5+iacto(ksb)*i1m                            6d19s22
                   bc(iad3)=bc(iad3)+bc(iad2)*ff                        6d19s22
                  end do                                                 12d22s17
                 end do                                                  12d22s17
                end do                                                   12d22s17
               end if                                                    12d22s17
               i1x=i1x+nrowi                                             12d22s17
              end do                                                     12d22s17
              i10=1                                                      12d22s17
             end do                                                      12d22s17
            else if(isblk1p(4,is).eq.isb.and.                                6d16s22
     $          isblkder(2,isd).eq.ksb.and.                             6d19s22
     $           isblk1p(3,is).eq.isblkder(1,isd).                      6d16s22    Avd=sum Gaa'a''a'
     $        and.isblk1p(1,is).eq.isblkder(3,isd).and.isblk1p(2,is).eq. 12d22s17   ji      i2 3  4
     $          isblkder(4,isd))then                                     6d16s22            21 2  1
c     Ava=fmul sum Ga'''aa'a''(a'a''|a'''v): bcp
c      ik               k                i
             call ilimts(noc(isblk1p(3,is)),nvirtc(isb),mynprocg,        6d9s22
     $          mynowprog,il,ih,i1s,i1e,i2s,i2e)                        12d22s17
             i10=i1s                                                     12d22s17
             i1n=noc(isblk1p(3,is))                                      12d22s17
             if(isblk1p(1,is).eq.isblk1p(2,is))then                     7d1s22
              nrowi=(noc(isblk1p(1,is))*(noc(isblk1p(1,is))+1))/2       7d1s22
              fmul=1d0                                                  7d1s22
              iswitch=0                                                 7d1s22
             else                                                       7d1s22
              nrowi=noc(isblk1p(1,is))*noc(isblk1p(2,is))                12d22s17
              fmul=2d0                                                   7d
              iswitch=1
             end if                                                     7d1s22
             if(isblkder(1,isd).ne.isblkder(2,isd).and.                 7d1s22
     $          isblkder(3,isd).ne.isblkder(4,isd))fmul=fmul*0.5d0      7d1s22
             i1x=ionex(is)                                               12d22s17
             do i2=i2s,i2e                                               12d22s17
              if(i2.eq.i2e)i1n=i1e                                       12d22s17
              i2m=i2-1                                                   12d22s17
              do i1=i10,i1n                                              12d22s17
               if(i1.gt.idoub(isblk1p(3,is)))then                        12d22s17
                i1m=i1-1-idoub(isblk1p(3,is))                            12d22s17
                do i4=0,iacto(isblk1p(2,is))-1                           12d22s17
                 i4p=i4+idoub(isblk1p(2,is))                             12d22s17
                 do i3=0,iacto(isblk1p(1,is))-1                          12d22s17
                  i3p=i3+idoub(isblk1p(1,is))                            12d22s17
                  ix=max(i3p,i4p)                                        12d22s17
                  in=min(i3p,i4p)                                        12d22s17
                  ieq=((ix*(ix+1))/2)+in                                 12d22s17
                  inot=i3p+noc(isblk1p(1,is))*i4p                        12d22s17
                  iad1=i1x+(inot-ieq)*iswitch+ieq                        12d22s17
                  ff=bc(iad1)*fmul                                      6d19s22
                  inot=i3+iacto(isblkder(3,isd))*i4                     6d19s22
                  j2=j2den(isd)+nrowd*inot                              6d16s22
                  do i5=0,iacto(ksb)-1                                   12d22s17
                   iad3=jamat+i2m+nvirtc(isb)*(i5+idoub(ksb))           6d19s22
                   iad2=j2+i1m+iacto(isblk1p(3,is))*i5                  6d20s22
                   bc(iad3)=bc(iad3)+bc(iad2)*ff                        6d19s22
                  end do                                                 12d22s17
                 end do                                                  12d22s17
                end do                                                   12d22s17
               end if                                                    12d22s17
               i1x=i1x+nrowi                                             12d22s17
              end do                                                     12d22s17
              i10=1                                                      12d22s17
             end do                                                      12d22s17
            else if(isblk1p(4,is).eq.isb.and.                                6d16s22
     $          isblkder(4,isd).eq.ksb.and.                             6d19s22
     $           isblk1p(3,is).eq.isblkder(3,isd).                      6d16s22
     $        and.isblk1p(1,is).eq.isblkder(1,isd).and.isblk1p(2,is).eq. 12d22s17
     $          isblkder(2,isd))then                                     6d16s22
c     Ava=fmul sum Ga'a''a'''a(a'a''|a'''v): bcp2
c      ik                    k           i
             call ilimts(noc(isblk1p(3,is)),nvirtc(isb),mynprocg,        6d9s22
     $          mynowprog,il,ih,i1s,i1e,i2s,i2e)                        12d22s17
             i10=i1s                                                     12d22s17
             i1n=noc(isblk1p(3,is))                                      12d22s17
             if(isblk1p(1,is).eq.isblk1p(2,is))then                     7d1s22
              nrowi=(noc(isblk1p(1,is))*(noc(isblk1p(1,is))+1))/2       7d1s22
              fmul=1d0                                                  7d1s22
              iswitch=0                                                 7d1s22
             else                                                       7d1s22
              nrowi=noc(isblk1p(1,is))*noc(isblk1p(2,is))                12d22s17
              fmul=2d0                                                   7d
              iswitch=1                                                 7d1s22
             end if                                                     7d1s22
             if(isblkder(1,isd).ne.isblkder(2,isd).and.                 7d1s22
     $          isblkder(3,isd).ne.isblkder(4,isd))fmul=fmul*0.5d0      7d1s22
             i1x=ionex(is)                                               12d22s17
             do i2=i2s,i2e                                               12d22s17
              if(i2.eq.i2e)i1n=i1e                                       12d22s17
              i2m=i2-1                                                   12d22s17
              do i1=i10,i1n                                              12d22s17
               if(i1.gt.idoub(isblk1p(3,is)))then                        12d22s17
                i1m=i1-1-idoub(isblk1p(3,is))                            12d22s17
                do i4=0,iacto(isblk1p(2,is))-1                           12d22s17
                 i4p=i4+idoub(isblk1p(2,is))                             12d22s17
                 do i3=0,iacto(isblk1p(1,is))-1                          12d22s17
                  i3p=i3+idoub(isblk1p(1,is))                            12d22s17
                  inot=i3p+noc(isblk1p(1,is))*i4p                        12d22s17
                  ix=max(i3p,i4p)                                       7d1s22
                  in=min(i3p,i4p)                                       7d1s22
                  itri=((ix*(ix+1))/2)+in                               7d1s22
                  iad1=i1x+iswitch*(inot-itri)+itri                     7d1s22
                  ff=bc(iad1)*fmul                                      6d19s22
                  inot=i3+iacto(isblkder(1,isd))*i4                     7d1s22
                  j2=j2den(isd)+inot                                    7d1s22
                  do i5=0,iacto(ksb)-1                                   12d22s17
                   iad3=jamat+i2m+nvirtc(isb)*(i5+idoub(ksb))           6d19s22
                   iad2=j2+nrowd*(i1m+iacto(isblk1p(3,is))*i5)          7d1s22
                   bc(iad3)=bc(iad3)+bc(iad2)*ff                        6d19s22
                  end do                                                 12d22s17
                 end do                                                  12d22s17
                end do                                                   12d22s17
               end if                                                    12d22s17
               i1x=i1x+nrowi                                             12d22s17
              end do                                                     12d22s17
              i10=1                                                      12d22s17
             end do                                                      12d22s17
            else if(isblk1p(4,is).eq.isb.and.                                6d16s22
     $          isblkder(3,isd).eq.ksb.and.                             7d1s22
     $           isblk1p(3,is).eq.isblkder(4,isd).                      7d1s22
     $        and.isblk1p(1,is).eq.isblkder(1,isd).and.isblk1p(2,is).eq. 12d22s17
     $          isblkder(2,isd))then                                     6d16s22
c     Ava=fmul sum Ga'a''aa'''(a'a''|a'''v): bcp3
c      ik                k               i
             call ilimts(noc(isblk1p(3,is)),nvirtc(isb),mynprocg,        6d9s22
     $          mynowprog,il,ih,i1s,i1e,i2s,i2e)                        12d22s17
             i10=i1s                                                     12d22s17
             i1n=noc(isblk1p(3,is))                                      12d22s17
             if(isblk1p(1,is).eq.isblk1p(2,is))then                     7d1s22
              nrowi=(noc(isblk1p(1,is))*(noc(isblk1p(1,is))+1))/2       7d1s22
              fmul=1d0                                                  7d1s22
              iswitch=0                                                 7d1s22
             else                                                       7d1s22
              nrowi=noc(isblk1p(1,is))*noc(isblk1p(2,is))                12d22s17
              fmul=2d0                                                   7d
              iswitch=1                                                 7d1s22
             end if                                                     7d1s22
             if(isblkder(1,isd).ne.isblkder(2,isd).and.                 7d1s22
     $          isblkder(3,isd).ne.isblkder(4,isd))fmul=fmul*0.5d0      7d1s22
             i1x=ionex(is)                                               12d22s17
             do i2=i2s,i2e                                               12d22s17
              if(i2.eq.i2e)i1n=i1e                                       12d22s17
              i2m=i2-1                                                   12d22s17
              do i1=i10,i1n                                              12d22s17
               if(i1.gt.idoub(isblk1p(3,is)))then                        12d22s17
                i1m=i1-1-idoub(isblk1p(3,is))                            12d22s17
                do i4=0,iacto(isblk1p(2,is))-1                           12d22s17
                 i4p=i4+idoub(isblk1p(2,is))                             12d22s17
                 do i3=0,iacto(isblk1p(1,is))-1                          12d22s17
                  i3p=i3+idoub(isblk1p(1,is))                            12d22s17
                  inot=i3p+noc(isblk1p(1,is))*i4p                        12d22s17
                  ix=max(i3p,i4p)                                       7d1s22
                  in=min(i3p,i4p)                                       7d1s22
                  itri=((ix*(ix+1))/2)+in                               7d1s22
                  iad1=i1x+iswitch*(inot-itri)+itri                     7d1s22
                  ff=bc(iad1)*fmul                                      6d19s22
                  inot=i3+iacto(isblkder(1,isd))*i4                     7d1s22
                  j2=j2den(isd)+inot                                    7d1s22
                  do i5=0,iacto(ksb)-1                                   12d22s17
                   iad3=jamat+i2m+nvirtc(isb)*(i5+idoub(ksb))           6d19s22
                   iad2=j2+nrowd*(i5+iacto(ksb)*i1m)                    7d1s22
                   bc(iad3)=bc(iad3)+bc(iad2)*ff                        6d19s22
                  end do                                                 12d22s17
                 end do                                                  12d22s17
                end do                                                   12d22s17
               end if                                                    12d22s17
               i1x=i1x+nrowi                                             12d22s17
              end do                                                     12d22s17
              i10=1                                                      12d22s17
             end do                                                      12d22s17
            end if                                                       12d22s17
           end do                                                        12d22s17
          end if                                                        6d16s22
         end do                                                         6d16s22
        end if                                                          6d16s22
       end do                                                           12d21s17
      end if                                                            12d21s17
      if(idoit(3).ne.0.and.mynowprog.eq.mynprocg-1)then                 12d15s17
c
c     aa part
c
       jh0=ih0
       do isb=1,nsymb                                                   12d13s17
        jsb=multh(isb,ipuse)                                            6d16s22
        jsd=multh(isb,iduse)                                            6d16s22
        jamat=iamat(jsd)
        nb=noc(isb)+nvirtc(isb)
        nj=noc(jsb)+nvirtc(jsb)                                         6d8s22
        nd=noc(jsd)+nvirtc(jsd)                                         6d16s22
        kh0=jh0+idoub(isb)
        do id=0,idoub(jsb)-1                                            6d7s22
         jden1=iden1(isb)
         do ia=0,iacto(jsd)-1                                           6d16s22
          do iap=0,iacto(isb)-1                                         12d13s17
           bc(jamat+ia)=bc(jamat+ia)-2d0*bc(kh0+iap)*bc(jden1+iap)      12d23s17      di         i i
          end do                                                        12d13s17
          jden1=jden1+iacto(isb)
         end do
         kh0=kh0+nb
         jamat=jamat+iacto(jsd)                                         6d16s22
        end do                                                          6d16s22
        if(iduse.eq.1)then                                              6d17s22
         jamat=iamat(isb)+idoub(jsb)*iacto(isb)
         jamat=jamat+idoub(jsb)*nvirtc(isb)                              6d16s22
        else                                                            6d17s22
         jamat=iamat(isb)+idoub(jsd)*iacto(isb)
         jamat=jamat+idoub(jsd)*nvirtc(isb)                              6d16s22
        end if                                                          6d17s22
        jden=iden1(jsb)                                                 6d9s22 Ava=2 sum Hva' Da'a
        do io=0,iacto(multh(iduse,jsb))-1                               6d7s22  ij        ij   j j
         kh0=jh0+noc(isb)+nb*idoub(jsb)                                 6d9s22  id        ii   i d
         do ia=0,iacto(jsb)-1                                           6d9s22
          do iv=0,nvirtc(isb)-1
           bc(jamat+iv)=bc(jamat+iv)+bc(kh0+iv)*bc(jden)*2d0
          end do
          kh0=kh0+nb
          jden=jden+1
         end do
         jamat=jamat+nvirtc(isb)                                        12d15s17
        end do
        jh0=jh0+nb*nj                                                   6d8s22
       end do                                                           12d13s17
      end if
      call dws_gsumf(bc(iamat(1)),natot)                                12d15s17
      do isb=1,nsymb                                                    4d2s18
       jsb=multh(isb,multh(iduse,ipuse))                                6d16s22
       rmssz(isb)=0d0                                                   4d2s18
       nda=idoub(jsb)*iacto(isb)                                        6d7s22
       nov=noc(jsb)*nvirtc(isb)                                         6d7s22
       if(nda+nov.gt.0)then                                             4d2s18
        do i=0,nda+nov-1                                                 4d2s18
         rmssz(isb)=rmssz(isb)+bc(iamat(isb)+i)**2                       4d2s18
        end do                                                           4d2s18
        rmssz(isb)=sqrt(rmssz(isb)/dfloat(nda+nov))                     4d2s18
       end if                                                           4d2s18
      end do                                                            4d2s18
      if(idwsdeb.ge.10)then
      write(6,*)('rms amat size in buildcasgrad: '),(rmssz(isb),isb=1,  4d2s18
     $     nsymb),nsdlk1p                                                       4d2s18
       write(6,*)('amat from buildcasgrad: '),nsdlk1p
       write(6,*)('goal '),bc(igoal)
       do isb=1,nsymb
        jsb=multh(isb,multh(iduse,ipuse))                               6d16s22
        write(6,*)('for symmetry block '),isb,jsb,nsdlk1p               6d7s22
        if(iacto(isb)*idoub(jsb).gt.0)then
         write(6,*)('active-doub part '),iamat(isb)
         call prntm2(bc(iamat(isb)),iacto(isb),idoub(jsb),iacto(isb))
         if(nsymb.eq.1)then
          itimes2=ibcoff                                                1d2s23
          ibcoff=itimes2+iacto(isb)*idoub(jsb)                          1d2s23
          do ix2x=0,iacto(isb)*idoub(jsb)-1                             1d2s23
           bc(itimes2+ix2x)=bc(iamat(isb)+ix2x)*2d0                     1d2s23
          end do                                                        1d2s23
          call printa(bc(itimes2),iacto,idoub,1,0,idoub,0,1,0,
     $         bc(ibcoff))
          ibcoff=itimes2                                                1d2s23
         end if
        end if
        jamat=iamat(isb)+iacto(isb)*idoub(jsb)                          6d7s22
        if(noc(jsb)*nvirtc(isb).gt.0)then                               6d7s22
         write(6,*)('occupied-virt part '),jamat
         call prntm2(bc(jamat),nvirtc(isb),noc(jsb),nvirtc(isb))        6d7s22
         if(nsymb.eq.1)then
          itimes2=ibcoff
          ibcoff=itimes2+nvirtc(isb)*noc(jsb)
          write(6,*)('times2 ')
          do ix2x=0,nvirtc(isb)*noc(jsb)-1
           bc(itimes2+ix2x)=bc(jamat+ix2x)*2d0
          end do
          call printa(bc(itimes2),nvirtc,noc,1,0,noc,0,1,0,bc(ibcoff))
          ibcoff=itimes2
         end if
        end if
       end do
      end if
      iamatb(1)=ibcoff                                                  2d16s18
      ibcoff=iamatb(1)+natot                                            2d16s18
      call enough('buildcasgrad.  2',bc,ibc)
      do i=0,natot-1                                                    2d16s18
       bc(iamatb(1)+i)=bc(iamat(1)+i)                                   2d16s18
      end do                                                            2d16s18
      return
      end
