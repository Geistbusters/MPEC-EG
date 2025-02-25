c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine buildhesscas(ioooo,ionex,jmats,kmats,h0mo,noc,ihessa,  4d25s17
     $     ihessb,ihessc,ihessd,                                        12d4s17
     $     nvirtc,ipsym,multh,iden1,jdenpt,idoub,iacto,idwsdeb,bc,ibc)  11d9s22
      implicit real*8 (a-h,o-z)
c
c     build orbital rotation hessian matrix
c     ihessa: rotations between doubly occupied and active,
c             both bra and ket
c     ihessb: bra: rotations between doubly occupied and active, ket:
c             rotations between occupied and virtual
c             dimensions: active, doub, occ, virt
c     ihessc: rotations between occupied and virtual, both bra and ket
c     ihessd: diagonal of hessian                                       12d4s17
c
      logical lprint,lsym                                               11d30s17
      parameter (idgoalx=2)
      dimension h0mo(1),ioooo(1),ionex(1),jmats(1),kmats(1),noc(1),
     $     ihessc(8,8),nvirtc(1),multh(8,8),ipt(8),iden1(8),jdenpt(*),   4d11s17
     $     idoub(8),iacto(8),ihessa(8,8),ihessb(8,8),idoit(7),          12d4s17
     $     ihessd(8)                                                    12d4s17
      dimension igoalx(idgoalx)                                         10d10s17
      common/lowersymcm/nsymbgx,iptno(8),ipts(8),nhsz(8),ipao(8)        4d25s18
      include "common.hf"
      include "common.store"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      common/unitcm/iunit                                               11d9s17
      data idoit/7*1/                                                   9d25s17
      idoit(4)=1
      iunit=6                                                           11d9s17
      itest=0                                                           6d19s17
      ioffh=0
      do isb=1,nsymb
       ipt(isb)=ioffh
       nh=noc(isb)+nvirtc(isb)
       ioffh=ioffh+nh*nh
      end do
      do isb=1,nsymb
       do isa=1,nsymb                                                   11d30s17
        lsym=isa.le.isb                                                 11d30s17
c
c     the symmetry label isb,isa correspond to the bra and ket
c     occupied orbital symmetries. note that I've coded up things so the
c     order of occupied orbital indicies is bra,ket but the order of
c     virtual orbital indicies is ket,bra.
c     because of the distribution across processors, we need
c     isb,isa to be the virtual labels, which is only different when
c     the perturbation is not totally symmetric.                        3d29s16
c
        ihessa(isa,isb)=-1                                              4d20s17
        ihessb(isa,isb)=-1                                              4d20s17
        ihessc(isa,isb)=-1                                               12d5s16
        isav=isa                                                        3d29s16
        isbv=isb                                                        3d29s16
        isao=multh(isa,ipsym)                                           3d28s16
        isbo=multh(isb,ipsym)                                           3d28s16
        if(noc(isbo)*noc(isao).gt.0)then                                  3d28s16
         if(isa.eq.isb)then                                             12d26s17
          ihessd(isa)=ibcoff                                            12d4s17
          nhessd=idoub(isao)*iacto(isa)+noc(isao)*nvirtc(isa)           12d4s17
          ibcoff=ihessd(isa)+nhessd                                     12d4s17
          call enough('buildhesscas.  1',bc,ibc)
          do i=0,nhessd-1                                               12d4s17
           bc(ihessd(isa)+i)=0d0                                        12d4s17
          end do                                                        12d4s17
         end if                                                         12d4s17
         if(lsym)then                                                   12d26s17
          ihessc(isa,isb)=ibcoff                                           3d28s16
          call ilimts(nvirtc(isb),nvirtc(isa),mynprocg,mynowprog,ilh,    3d29s16
     $        ihh,                                                      3d28s16
     $        i1sh,i1eh,i2sh,i2eh)                                      3d28s16
          nhere=ihh+1-ilh                                                3d28s16
          nrowh=noc(isao)*noc(isbo)                                        3d28s16
          ibcoff=ihessc(isa,isb)+nhere*noc(isao)*noc(isbo)                  3d28s16
          nrowa=idoub(isao)*idoub(isbo)                                  4d20s17
          ncola=iacto(isb)*iacto(isa)                                    4d20s17
          ihessa(isa,isb)=ibcoff                                         4d20s17
          ibcoff=ihessa(isa,isb)+nrowa*ncola                             4d20s17
          do i1234=0,nhere*nrowh+nrowa*ncola-1                           11d30s17
           bc(ihessc(isa,isb)+i1234)=0d0                                  4d11s17
          end do                                                         4d11s17
         end if
         call ilimts(noc(isbo),nvirtc(isb),mynprocg,mynowprog,ilb,ihb,  4d25s17
     $        i1sb,i1eb,i2sb,i2eb)                                      4d25s17
         nhereb=ihb+1-ilb                                               4d25s17
         nrowb=idoub(isao)*iacto(isa)                                   4d25s17
         ihessb(isa,isb)=ibcoff                                         4d25s17
         ibcoff=ihessb(isa,isb)+nrowb*nhereb                            4d25s17
         call enough('buildhesscas.  2',bc,ibc)
         do i1234=0,nrowb*nhereb-1                                      11d30s17
          bc(ihessb(isa,isb)+i1234)=0d0                                  4d11s17
         end do                                                         4d11s17
        end if                                                          12d26s17
       end do                                                           12d26s17
      end do                                                            12d26s17
      do isb=1,nsymb
       do isa=1,nsymb                                                   11d30s17
        lsym=isa.le.isb                                                 11d30s17
        if(idwsdeb.gt.10)write(6,*)('for isa,isb = '),isa,isb,lsym      11d30s17
c
c     the symmetry label isb,isa correspond to the bra and ket
c     occupied orbital symmetries. note that I've coded up things so the
c     order of occupied orbital indicies is bra,ket but the order of
c     virtual orbital indicies is ket,bra.
c     because of the distribution across processors, we need
c     isb,isa to be the virtual labels, which is only different when
c     the perturbation is not totally symmetric.                        3d29s16
c
        isav=isa                                                        3d29s16
        isbv=isb                                                        3d29s16
        isao=multh(isa,ipsym)                                           3d28s16
        isbo=multh(isb,ipsym)                                           3d28s16
        if(noc(isbo)*noc(isao).gt.0)then                                  3d28s16
         if(lsym)then                                                   11d30s17
          if(isa.eq.isb)then                                            12d4s17
           nhessd=idoub(isao)*iacto(isa)+noc(isao)*nvirtc(isa)          12d4s17
          end if                                                        12d4s17
          call ilimts(nvirtc(isb),nvirtc(isa),mynprocg,mynowprog,ilh,    3d29s16
     $        ihh,                                                      3d28s16
     $        i1sh,i1eh,i2sh,i2eh)                                      3d28s16
          nhere=ihh+1-ilh                                                3d28s16
          nrowh=noc(isao)*noc(isbo)                                        3d28s16
          ig1=-1
          ig2=-1
          if(isb.eq.1.and.isa.eq.1)then
          ig1=ihessc(isa,isb)
          end if
          nrowa=idoub(isao)*idoub(isbo)                                  4d20s17
          ncola=iacto(isb)*iacto(isa)                                    4d20s17
         end if                                                         11d30s17
         call ilimts(noc(isbo),nvirtc(isb),mynprocg,mynowprog,ilb,ihb,  4d25s17
     $        i1sb,i1eb,i2sb,i2eb)                                      4d25s17
         nhereb=ihb+1-ilb                                               4d25s17
         nrowb=idoub(isao)*iacto(isa)                                   4d25s17
         if(idwsdeb.gt.10)
     $        write(6,*)('space for ihessc: '),noc(isao)*noc(isbo),nhere
         ibtop=ibcoff                                                   3d29s16
         call enough('buildhesscas.  3',bc,ibc)
         igotj=0                                                        4d11s17
c
c     for rotations between doub&act vs doub&virt
c     because of distribution on onex, we need to work on full matrices
c     here for J part.
c
         ihessbt=ibcoff                                                 4d25s17
         nwds=iacto(isa)*idoub(isao)*nvirtc(isb)*noc(isbo)              9d19s17
         itest=itest+1
         if(nwds.gt.0)then                                              4d25s17
          ibcoff=ihessbt+nwds                                           4d25s17
          call enough('buildhesscas.  4',bc,ibc)
          do i=0,nwds-1                                                 4d25s17
           bc(ihessbt+i)=0d0                                            4d25s17
          end do                                                        4d25s17
          do is=1,nsdlk1                                                 4d25s17
c
c     ddaa: nd,v* part e
c
           if(idoit(4).ne.0)then                                        9d25s17
           if(isblk1(3,is).eq.isbo.and.isblk1(4,is).eq.isbv)then        9d25s17
            if(isblk1(1,is).eq.isao)then                                9d18s17
             iswitch=iabs(isao-isav)                                    9d18s17
             iswitch=iswitch/max(1,iswitch)                             9d18s17
             i10=i1sb                                                   9d18s17
             i1n=noc(isbo)                                              9d18s17
             if(isao.eq.isav)then                                       9d18s17
              nrow=(noc(isao)*(noc(isao)+1))/2                          9d18s17
             else                                                       9d18s17
              nrow=noc(isao)*noc(isav)                                  9d18s17
             end if                                                     9d18s17
             i1x=ionex(is)                                              9d18s17
             iadh=ihessb(isa,isb)                                       9d18s17
             do i2=i2sb,i2eb                                            9d18s17
              if(i2.eq.i2eb)i1n=i1eb                                    9d18s17
              i2m=i2-1
              do i1=i10,i1n                                             9d18s17
               if(i1.le.idoub(isbo))then                                9d18s17
                do i4=0,idoub(isao)-1                                   9d18s17
                 do i3=0,iacto(isav)-1                                  9d18s17
                  i3p=i3+idoub(isav)                                    9d18s17
                  in=min(i3p,i4)                                        9d18s17
                  ix=max(i3p,i4)                                        9d18s17
                  ieq=((ix*(ix+1))/2)+in                                9d18s17
                  inot=i4+noc(isao)*i3p                                 9d18s17
                  iadi=i1x+(inot-ieq)*iswitch+ieq                        9d18s17
                  jadh=iadh+iacto(isav)*i4                              9d18s17
                  iden=iden1(isav)+iacto(isav)*i3                       9d18s17
                  fact=-8d0*bc(iadi)                                    9d22s17
                  do i5=0,iacto(isav)-1                                 9d18s17
                   bc(jadh+i5)=bc(jadh+i5)+fact*bc(iden+i5)              9d18s17
                  end do
                 end do                                                 9d18s17
                end do                                                  9d18s17
               else                                                     9d18s17
                i1m=i1-1                                                9d19s17
                jhessbt=ihessbt+nrowb*(idoub(isbo)+noc(isbo)*i2m)       9d20s17
                i1m=i1m-idoub(isbo)                                     9d19s17
                iden=iden1(isbo)+iacto(isbo)*i1m                        9d18s17
                do i4=0,idoub(isao)-1                                   9d18s17
                 do i3=0,iacto(isav)-1                                  9d18s17
                  i3p=i3+idoub(isav)                                    9d18s17
                  in=min(i4,i3p)                                        9d18s17
                  ix=max(i4,i3p)                                        9d18s17
                  ieq=((ix*(ix+1))/2)+in                                9d18s17
                  inot=i4+noc(isao)*i3p                                 11d9s17
                  iadi=i1x+(inot-ieq)*iswitch+ieq                       9d18s17
                  fact=bc(iadi)*8d0                                     9d18s17
                  jadh=jhessbt+i3+iacto(isav)*i4                        9d19s17
                  do i5=0,iacto(isbo)-1                                 9d18s17
                   bc(jadh+i5*nrowb)=bc(jadh+i5*nrowb)+fact*bc(iden+i5) 9d20s17
                  end do                                                9d18s17
                 end do                                                 9d18s17
                end do                                                  9d18s17
               end if                                                   9d18s17
               i1x=i1x+nrow                                             9d18s17
               iadh=iadh+nrowb                                          9d18s17
              end do                                                    9d18s17
              i10=1                                                     9d18s17
             end do                                                     9d18s17
            else if(isblk1(2,is).eq.isao)then                           9d18s17
             i10=i1sb                                                   9d18s17
             i1n=noc(isbo)                                              9d18s17
             nrow=noc(isao)*noc(isav)                                   9d18s17
             i1x=ionex(is)                                              9d18s17
             iadh=ihessb(isa,isb)                                       9d18s17
             do i2=i2sb,i2eb                                            9d18s17
              i2m=i2-1                                                  9d19s17
              if(i2.eq.i2eb)i1n=i1eb                                    9d18s17
              do i1=i10,i1n                                             9d18s17
               if(i1.le.idoub(isbo))then                                9d18s17
                do i4=0,idoub(isao)-1                                   9d18s17
                 do i3=0,iacto(isav)-1                                  9d18s17
                  i3p=i3+idoub(isav)                                    9d18s17
                  inot=i3p+noc(isav)*i4                                  9d18s17
                  iadi=i1x+inot                                         9d18s17
                  jadh=iadh+iacto(isav)*i4                              9d18s17
                  iden=iden1(isav)+iacto(isav)*i3                       9d18s17
                  fact=-8d0*bc(iadi)                                    9d22s17
                  do i5=0,iacto(isav)-1                                 9d18s17
                   bc(jadh+i5)=bc(jadh+i5)+fact*bc(iden+i5)              9d18s17
                  end do
                 end do                                                 9d18s17
                end do                                                  9d18s17
               else                                                     9d18s17
                i1m=i1-1                                                9d19s17
                jhessbt=ihessbt+nrowb*(idoub(isbo)+noc(isbo)*i2m)       9d20s17
                i1m=i1m-idoub(isbo)                                     9d19s17
                iden=iden1(isbo)+iacto(isbo)*i1m                        9d18s17
                do i4=0,idoub(isao)-1                                   9d18s17
                 do i3=0,iacto(isav)-1                                  9d18s17
                  i3p=i3+idoub(isav)                                    9d18s17
                  inot=i3p+noc(isav)*i4                                 11d9s17
                  iadi=i1x+inot                                         9d18s17
                  fact=bc(iadi)*8d0                                     9d18s17
                  jadh=jhessbt+i3+iacto(isav)*i4                        9d19s17
                  do i5=0,iacto(isbo)-1                                 9d18s17
                   bc(jadh+i5*nrowb)=bc(jadh+i5*nrowb)+fact*bc(iden+i5) 9d20s17
                  end do                                                9d18s17
                 end do                                                 9d18s17
                end do                                                  9d18s17
               end if                                                   9d18s17
               i1x=i1x+nrow                                             9d18s17
               iadh=iadh+nrowb                                          9d18s17
              end do                                                    9d18s17
              i10=1                                                     9d18s17
             end do                                                     9d18s17
            end if                                                      9d18s17
           end if                                                       9d18s17
c
c     ddaa: nd,v* part f
c
           if(isblk1(3,is).eq.isao.and.isblk1(4,is).eq.isbv)then        9d18s17
            call ilimts(noc(isao),nvirtc(isbv),mynprocg,mynowprog,il,   9d20s17
     $          ih,i1s,i1e,i2s,i2e)                                     9d20s17
            if(isblk1(1,is).eq.isav)then                                9d18s17
             iswitch=iabs(isbo-isav)                                    9d18s17
             iswitch=iswitch/max(1,iswitch)                             9d18s17
             i10=i1s                                                    9d20s17
             i1n=noc(isao)                                              9d18s17
             if(isbo.eq.isav)then                                       9d18s17
              nrow=(noc(isbo)*(noc(isbo)+1))/2                          11d9s17
             else                                                       9d18s17
              nrow=noc(isbo)*noc(isav)                                  9d18s17
             end if                                                     9d18s17
             i1x=ionex(is)                                              9d18s17
             iadh=ihessb(isa,isb)                                       9d18s17
             do i2=i2s,i2e                                              9d20s17
              if(i2.eq.i2e)i1n=i1e                                      9d20s17
              i2m=i2-1
              do i1=i10,i1n                                             9d18s17
               if(i1.le.idoub(isao))then                                9d18s17
                i1m=i1-1
                jhessbt=ihessbt+iacto(isav)*(i1m                        9d20s17
     $               +idoub(isao)*noc(isbo)*i2m)                          9d20s17
                do i4=0,idoub(isbo)-1                                   9d18s17
                 jadh=jhessbt+nrowb*i4                                  9d20s17
                 do i3=0,iacto(isav)-1                                  9d18s17
                  i3p=i3+idoub(isav)                                    9d18s17
                  in=min(i3p,i4)                                        9d18s17
                  ix=max(i3p,i4)                                        9d18s17
                  ieq=((ix*(ix+1))/2)+in                                9d18s17
                  inot=i3p+noc(isav)*i4                                 9d20s17
                  iadi=i1x+(inot-ieq)*iswitch+ieq                        9d18s17
                  iden=iden1(isav)+iacto(isav)*i3                       9d18s17
                  fact=2d0*bc(iadi)                                     9d22s17
                  do i5=0,iacto(isav)-1                                 9d18s17
                   bc(jadh+i5)=bc(jadh+i5)+fact*bc(iden+i5)              9d18s17
                  end do
                 end do                                                 9d18s17
                end do                                                  9d18s17
                jhessbt=ihessbt+iacto(isav)*(i1m+idoub(isao)*           10d24s17
     $               (idoub(isbo)+noc(isbo)*i2m))                       9d20s17
                do i4=0,iacto(isav)-1                                   9d18s17
                 i4p=i4+idoub(isav)                                     9d20s17
                 jadh=jhessbt+i4                                        9d20s17
                 do i3=0,iacto(isbo)-1                                  9d20s17
                  iden=iden1(isbo)+iacto(isbo)*i3                       9d20s17
                  i3p=i3+idoub(isbo)                                    9d20s17
                  in=min(i4p,i3p)                                        9d18s17
                  ix=max(i4p,i3p)                                        9d18s17
                  ieq=((ix*(ix+1))/2)+in                                9d18s17
                  inot=i4p+noc(isav)*i3p                                 9d20s17
                  iadi=i1x+(inot-ieq)*iswitch+ieq                       9d18s17
                  fact=-bc(iadi)*2d0                                     9d18s17
                  do i5=0,iacto(isbo)-1                                 9d18s17
                   iarg=jadh+i5*nrowb                                   9d20s17
                   bc(iarg)=bc(iarg)+fact*bc(iden+i5)                   9d20s17
                  end do                                                9d18s17
                 end do                                                 9d18s17
                end do                                                  9d18s17
               end if                                                   9d18s17
               i1x=i1x+nrow                                             9d18s17
               iadh=iadh+nrowb                                          9d18s17
              end do                                                    9d18s17
              i10=1                                                     9d18s17
             end do                                                     9d18s17
            else if(isblk1(2,is).eq.isav)then                           9d18s17
             i10=i1s                                                    9d20s17
             i1n=noc(isao)                                              9d20s17
             nrow=noc(isbo)*noc(isav)                                   9d20s17
             i1x=ionex(is)                                              9d18s17
             iadh=ihessb(isa,isb)                                       9d18s17
             do i2=i2s,i2e                                              9d20s17
              i2m=i2-1                                                  9d19s17
              if(i2.eq.i2e)i1n=i1e                                      9d20s17
              do i1=i10,i1n                                             9d18s17
               if(i1.le.idoub(isao))then                                9d18s17
                i1m=i1-1                                                9d20s17
                jhessbt=ihessbt+iacto(isav)*(i1m                        9d20s17
     $               +idoub(isao)*noc(isbo)*i2m)                        9d20s17
                do i4=0,idoub(isbo)-1                                   9d22s17
                 jadh=jhessbt+nrowb*i4                                  9d20s17
                 do i3=0,iacto(isav)-1                                  9d18s17
                  i3p=i3+idoub(isav)                                    9d18s17
                  inot=i4+noc(isbo)*i3p                                 11d8s17
                  iadi=i1x+inot                                         9d18s17
                  iden=iden1(isav)+iacto(isav)*i3                       9d18s17
                  fact=2d0*bc(iadi)                                     9d22s17
                  do i5=0,iacto(isav)-1                                 9d18s17
                   bc(jadh+i5)=bc(jadh+i5)+fact*bc(iden+i5)              9d18s17
                  end do
                 end do                                                 9d18s17
                end do                                                  9d18s17
                jhessbt=ihessbt+iacto(isav)*(i1m+idoub(isao)*           9d20s17
     $               (idoub(isbo)+noc(isbo)*i2m))                       9d20s17
                do i4=0,iacto(isav)-1                                   9d18s17
                 i4p=i4+idoub(isav)                                     9d20s17
                 jadh=jhessbt+i4                                        9d20s17
                 do i3=0,iacto(isbo)-1                                  9d20s17
                  iden=iden1(isbo)+iacto(isbo)*i3                       9d20s17
                  i3p=i3+idoub(isbo)                                    9d20s17
                  inot=i3p+noc(isbo)*i4p                                 9d20s17
                  iadi=i1x+inot                                         9d20s17
                  fact=-bc(iadi)*2d0                                     9d18s17
                  do i5=0,iacto(isbo)-1                                 9d18s17
                   bc(jadh+i5*nrowb)=bc(jadh+i5*nrowb)+fact*bc(iden+i5) 9d20s17
                  end do                                                9d18s17
                 end do                                                 9d18s17
                end do                                                  9d18s17
               end if                                                   9d18s17
               i1x=i1x+nrow                                             9d18s17
               iadh=iadh+nrowb                                          9d18s17
              end do                                                    9d18s17
              i10=1                                                     9d18s17
             end do                                                     9d18s17
            end if                                                      9d18s17
           end if                                                       9d18s17
c
c     ddaa: nd,v* part g
c
           if(isblk1(3,is).eq.isav.and.isblk1(4,is).eq.isbv)then        9d18s17
            call ilimts(noc(isav),nvirtc(isbv),mynprocg,mynowprog,il,   9d20s17
     $          ih,i1s,i1e,i2s,i2e)                                     9d20s17
            if(isblk1(1,is).eq.isao)then                                9d18s17
             iswitch=iabs(isbo-isao)                                    9d18s17
             iswitch=iswitch/max(1,iswitch)                             9d18s17
             i10=i1s                                                    9d20s17
             i1n=noc(isav)                                              9d18s17
             if(isbo.eq.isao)then                                       9d18s17
              nrow=(noc(isao)*(noc(isao)+1))/2                          9d18s17
             else                                                       9d18s17
              nrow=noc(isbo)*noc(isao)                                  9d18s17
             end if                                                     9d18s17
             i1x=ionex(is)                                              9d18s17
             iadh=ihessb(isa,isb)                                       9d18s17
             do i2=i2s,i2e                                              9d20s17
              if(i2.eq.i2e)i1n=i1e                                      9d20s17
              i2m=i2-1
              do i1=i10,i1n                                             9d18s17
               if(i1.gt.idoub(isav))then                                9d18s17
                i1m=i1-1-idoub(isav)                                    9d20s17
                jhessbt=ihessbt+nrowb*noc(isbo)*i2m                     9d20s17
                iden=iden1(isav)+iacto(isav)*i1m                        9d20s17
                do i4=0,idoub(isbo)-1                                   9d18s17
                 jadh=jhessbt+nrowb*i4                                  9d20s17
                 do i3=0,idoub(isao)-1                                  9d20s17
                  kadh=jadh+i3*iacto(isav)                              9d20s17
                  in=min(i3,i4)                                         9d20s17
                  ix=max(i3,i4)                                         9d20s17
                  ieq=((ix*(ix+1))/2)+in                                9d18s17
                  inot=i3+noc(isao)*i4                                  9d20s17
                  iadi=i1x+(inot-ieq)*iswitch+ieq                        9d18s17
                  fact=bc(iadi)*2d0                                     9d22s17
                  do i5=0,iacto(isav)-1                                 9d18s17
                   bc(kadh+i5)=bc(kadh+i5)+fact*bc(iden+i5)              9d18s17
                  end do
                 end do                                                 9d18s17
                end do                                                  9d18s17
                jhessbt=ihessbt+i1m+nrowb*(idoub(isbo)+noc(isbo)*i2m)   11d9s17
                do i4=0,idoub(isao)-1                                   9d22s17
                 jadh=jhessbt+iacto(isav)*i4                            9d22s17
                 do i3=0,iacto(isbo)-1                                  9d20s17
                  iden=iden1(isbo)+iacto(isbo)*i3                       9d20s17
                  i3p=i3+idoub(isbo)                                    9d20s17
                  in=min(i4,i3p)                                        9d22s17
                  ix=max(i4,i3p)                                        9d22s17
                  ieq=((ix*(ix+1))/2)+in                                9d18s17
                  inot=i4+noc(isao)*i3p                                 9d22s17
                  iadi=i1x+(inot-ieq)*iswitch+ieq                       9d18s17
                  fact=-bc(iadi)*2d0                                     9d18s17
                  do i5=0,iacto(isbo)-1                                 9d18s17
                   iarg=jadh+i5*nrowb                                   9d20s17
                   bc(iarg)=bc(iarg)+fact*bc(iden+i5)                   9d20s17
                  end do                                                9d18s17
                 end do                                                 9d18s17
                end do                                                  9d18s17
               end if                                                   9d18s17
               i1x=i1x+nrow                                             9d18s17
               iadh=iadh+nrowb                                          9d18s17
              end do                                                    9d18s17
              i10=1                                                     9d18s17
             end do                                                     9d18s17
            else if(isblk1(2,is).eq.isao)then                           9d18s17
             i10=i1s                                                    9d20s17
             i1n=noc(isav)                                              9d20s17
             nrow=noc(isbo)*noc(isao)                                   9d20s17
             i1x=ionex(is)                                              9d18s17
             iadh=ihessb(isa,isb)                                       9d18s17
             do i2=i2s,i2e                                              9d20s17
              i2m=i2-1                                                  9d19s17
              if(i2.eq.i2e)i1n=i1e                                      9d20s17
              do i1=i10,i1n                                             9d18s17
               if(i1.gt.idoub(isav))then                                9d18s17
                i1m=i1-1-idoub(isav)                                    9d20s17
                jhessbt=ihessbt+nrowb*noc(isbo)*i2m                     9d20s17
                iden=iden1(isav)+iacto(isav)*i1m                        9d20s17
                do i4=0,idoub(isbo)-1                                   9d18s17
                 jadh=jhessbt+nrowb*i4                                  9d20s17
                 do i3=0,idoub(isao)-1                                  9d20s17
                  kadh=jadh+i3*iacto(isav)                              9d20s17
                  inot=i4+noc(isbo)*i3                                  9d20s17
                  iadi=i1x+inot                                         9d20s17
                  fact=bc(iadi)*2d0                                     9d22s17
                  do i5=0,iacto(isav)-1                                 9d18s17
                   bc(kadh+i5)=bc(kadh+i5)+fact*bc(iden+i5)              9d18s17
                  end do
                 end do                                                 9d18s17
                end do                                                  9d18s17
                jhessbt=ihessbt+i1m+nrowb*(idoub(isbo)+noc(isbo)*i2m)   9d22s17
                do i4=0,idoub(isao)-1                                   9d22s17
                 jadh=jhessbt+iacto(isav)*i4                            9d22s17
                 do i3=0,iacto(isbo)-1                                  9d20s17
                  iden=iden1(isbo)+iacto(isbo)*i3                       9d20s17
                  i3p=i3+idoub(isbo)                                    9d20s17
                  inot=i3p+noc(isbo)*i4                                 9d22s17
                  iadi=i1x+inot                                         9d22s17
                  fact=-bc(iadi)*2d0                                     9d18s17
                  do i5=0,iacto(isbo)-1                                 9d18s17
                   iarg=jadh+i5*nrowb                                   9d20s17
                   bc(iarg)=bc(iarg)+fact*bc(iden+i5)                   9d20s17
                  end do                                                9d18s17
                 end do                                                 9d18s17
                end do                                                  9d18s17
               end if                                                   9d18s17
               i1x=i1x+nrow                                             9d18s17
               iadh=iadh+nrowb                                          9d18s17
              end do                                                    9d18s17
              i10=1                                                     9d18s17
             end do                                                     9d18s17
            end if                                                      9d18s17
            end if                                                      9d25s17
           end if                                                       9d18s17
c
c     dddd: "K" part
c
           if(idoit(3).ne.0)then                                        9d25s17
           igot=0
           if(isblk1(1,is).eq.isa.and.isblk1(2,is).eq.isao.and.
     $          isblk1(3,is).eq.isbo.and.isblk1(4,is).eq.isb)then        4d25s17
            igot=1
            i10=i1sb                                                      4d25s17
            i1n=noc(isbo)                                                4d25s17
            if(isa.eq.isao)then                                          4d25s17
             nrow=(noc(isa)*(noc(isa)+1))/2                              4d25s17
            else
             nrow=noc(isa)*noc(isao)                                     4d25s17
            end if                                                       4d25s17
            iox=ionex(is)                                                4d25s17
            ihb=ihessb(isa,isb)                                          4d25s17
            do i2=i2sb,i2eb                                              4d25s17
             if(i2.eq.i2eb)i1n=i1eb                                      4d25s17
             do i1=i10,i1n                                               4d25s17
              if(i1.le.idoub(isbo))then                                  4d25s17
               if(isa.eq.isao)then                                       4d25s17
                do i4=0,idoub(isao)-1                                    4d25s17
                 do i3=0,iacto(isa)-1                                    4d25s17
                  i3p=i3+idoub(isa)                                     5d4s17
                  ix=max(i4,i3p)                                         4d25s17
                  in=min(i4,i3p)                                         4d25s17
                  iadb=ihb+i3+iacto(isa)*i4                              4d25s17
                  iadx=iox+((ix*(ix+1))/2)+in                            4d25s17
                  bc(iadb)=bc(iadb)+16d0*bc(iadx)                       7d19s17
                 end do                                                  4d25s17
                end do                                                   4d25s17
               else                                                      4d25s17
                do i4=0,idoub(isao)-1                                     4d25s17
                 do i3=0,iacto(isa)-1                                     4d25s17
                  i3p=i3+idoub(isa)                                     5d4s17
                  iadb=ihb+i3+iacto(isa)*i4                               4d25s17
                  iadx=iox+i3p+noc(isa)*i4                                4d25s17
                  bc(iadb)=bc(iadb)+16d0*bc(iadx)                       7d19s17
                 end do                                                   4d25s17
                end do                                                    4d25s17
               end if                                                    4d25s17
              end if                                                     4d25s17
              iox=iox+nrow                                               4d25s17
              ihb=ihb+nrowb                                              4d25s17
             end do                                                      4d25s17
             i10=1                                                       4d25s17
            end do                                                       4d25s17
           end if                                                        4d25s17
           if(isblk1(2,is).eq.isa.and.isblk1(1,is).eq.isao.and.          4d25s17
     $        isblk1(3,is).eq.isbo.and.isblk1(4,is).eq.isb.and.         4d25s17
     $        igot.eq.0)then                                            4d25s17
            igot=1
            i10=i1sb                                                      4d25s17
            i1n=noc(isbo)                                                4d25s17
            nrow=noc(isa)*noc(isao)                                      4d25s17
            iox=ionex(is)                                                4d25s17
            ihb=ihessb(isa,isb)                                          4d25s17
            do i2=i2sb,i2eb                                                4d25s17
             if(i2.eq.i2eb)i1n=i1eb                                        4d25s17
             do i1=i10,i1n                                               4d25s17
              if(i1.le.idoub(isbo))then                                 5d4s17
               do i4=0,idoub(isao)-1                                     4d25s17
                do i3=0,iacto(isa)-1                                     4d25s17
                 i3p=i3+idoub(isa)                                       5d4s17
                 iadb=ihb+i3+iacto(isa)*i4                               4d25s17
                 iadx=iox+i4+noc(isao)*i3p                               4d25s17
                 bc(iadb)=bc(iadb)+16d0*bc(iadx)                        7d19s17
                end do                                                   4d25s17
               end do                                                    4d25s17
              end if                                                    5d4s17
              iox=iox+nrow                                               4d25s17
              ihb=ihb+nrowb                                              4d25s17
             end do                                                      4d25s17
             i10=1                                                       4d25s17
            end do                                                       4d25s17
           end if                                                        4d25s17
c
c     dddd: "J" part
c
           igot=0                                                        4d25s17
           if(isblk1(1,is).eq.isao.and.isblk1(2,is).eq.isbo.and.         4d25s17
     $         isblk1(3,is).eq.isa.and.isblk1(4,is).eq.isb)then           4d25s17
            igot=1                                                       4d25s17
            call ilimts(noc(isa),nvirtc(isb),mynprocg,mynowprog,il,ih,  4d25s17
     $           i1s,i1e,i2s,i2e)                                       4d25s17
            i10=i1s                                                     6d19s17
            i1n=noc(isa)                                                4d25s17
            iox=ionex(is)                                               4d25s17
            ihb=ihessbt                                                 4d25s17
            if(isao.eq.isbo)then                                        4d25s17
             nrow=(noc(isao)*(noc(isao)+1))/2                           4d25s17
            else                                                        4d25s17
             nrow=noc(isao)*noc(isbo)                                   4d25s17
            end if                                                      4d25s17
            do i2=i2s,i2e                                               4d25s17
             i2m=i2-1                                                   4d25s17
             if(i2.eq.i2e)i1n=i1e                                       4d25s17
             do i1=i10,i1n                                              4d25s17
              i1m=i1-1-idoub(isa)
              if(i1.gt.idoub(isa))then                                  4d25s17
               if(isao.eq.isbo)then                                     4d25s17
                do i4=0,idoub(isbo)-1                                    4d25s17
                 do i3=0,idoub(isao)-1                                   4d25s17
                  ix=max(i3,i4)                                         4d25s17
                  in=min(i3,i4)                                         4d25s17
                  iadx=iox+((ix*(ix+1))/2)+in                           4d25s17
                  iadb=ihessbt+i1m+iacto(isa)*(i3+idoub(isao)*          6d19s17
     $                 (i4+noc(isbo)*i2m))                              9d19s17
                  bc(iadb)=bc(iadb)-4d0*bc(iadx)                        7d19s17
                 end do                                                  4d25s17
                end do                                                   4d25s17
               else                                                     4d25s17
                do i4=0,idoub(isbo)-1                                    4d25s17
                 do i3=0,idoub(isao)-1                                   4d25s17
                  iadx=iox+i3+noc(isao)*i4                              4d25s17
                  iadb=ihessbt+i1m+iacto(isa)*(i3+idoub(isao)*          6d19s17
     $                 (i4+noc(isbo)*i2m))                              9d19s17
                  bc(iadb)=bc(iadb)-4d0*bc(iadx)                        7d19s17
                 end do                                                  4d25s17
                end do                                                   4d25s17
               end if                                                   4d25s17
              end if                                                    4d25s17
              iox=iox+nrow                                              4d25s17
             end do                                                     4d25s17
             i10=1                                                      4d25s17
            end do                                                      4d25s17
           end if                                                        4d25s17
           if(isblk1(2,is).eq.isao.and.isblk1(1,is).eq.isbo.and.         4d25s17
     $         isblk1(3,is).eq.isa.and.isblk1(4,is).eq.isb.and.         4d25s17
     $          igot.eq.0)then                                          4d25s17
            igot=1                                                       4d25s17
            call ilimts(noc(isa),nvirtc(isb),mynprocg,mynowprog,il,ih,  4d25s17
     $           i1s,i1e,i2s,i2e)                                       4d25s17
            i10=i1s                                                     6d19s17
            i1n=noc(isa)                                                4d25s17
            iox=ionex(is)                                               4d25s17
            nrow=noc(isao)*noc(isbo)                                    4d25s17
            do i2=i2s,i2e                                               4d25s17
             i2m=i2-1                                                   4d25s17
             if(i2.eq.i2e)i1n=i1e                                       4d25s17
             do i1=i10,i1n                                              4d25s17
              i1m=i1-1-idoub(isa)
              if(i1.gt.idoub(isa))then                                  4d25s17
                call dumb(i1,i2)
               do i4=0,idoub(isbo)-1                                    4d25s17
                do i3=0,idoub(isao)-1                                   4d25s17
                 iadx=iox+i4+noc(isbo)*i3                               4d25s17
                 iadb=ihessbt+i1m+iacto(isa)*(i3+idoub(isao)*           6d19s17
     $                 (i4+noc(isbo)*i2m))                              9d19s17
                 bc(iadb)=bc(iadb)-4d0*bc(iadx)                         11d8s17
                end do                                                  4d25s17
               end do                                                   4d25s17
              end if                                                    4d25s17
              iox=iox+nrow                                              4d25s17
             end do                                                     4d25s17
             i10=1                                                      4d25s17
            end do                                                      4d25s17
           end if                                                        4d25s17
           igot=0                                                        4d25s17
           if(isblk1(1,is).eq.isa.and.isblk1(2,is).eq.isbo.and.         7d13s17
     $         isblk1(3,is).eq.isao.and.isblk1(4,is).eq.isb)then        7d13s17
            igot=1                                                       4d25s17
            call ilimts(noc(isao),nvirtc(isb),mynprocg,mynowprog,il,ih, 7d13s17
     $           i1s,i1e,i2s,i2e)                                       4d25s17
            i10=i1s                                                     6d19s17
            i1n=noc(isao)                                               7d13s17
            iox=ionex(is)                                               4d25s17
            ihb=ihessbt                                                 4d25s17
            if(isa.eq.isbo)then                                         7d13s17
             nrow=(noc(isa)*(noc(isa)+1))/2                             7d13s17
            else                                                        4d25s17
             nrow=noc(isa)*noc(isbo)                                    7d13s17
            end if                                                      4d25s17
            do i2=i2s,i2e                                               4d25s17
             i2m=i2-1                                                   4d25s17
             if(i2.eq.i2e)i1n=i1e                                       4d25s17
             do i1=i10,i1n                                              4d25s17
              i1m=i1-1                                                  7d13s17
              if(i1.le.idoub(isao))then                                 7d13s17
               if(isa.eq.isbo)then                                      7d13s17
                do i4=0,idoub(isbo)-1                                    4d25s17
                 do i3=0,iacto(isa)-1                                   7d13s17
                  i3p=i3+idoub(isa)                                     7d13s17
                  ix=max(i3p,i4)                                        7d13s17
                  in=min(i3p,i4)                                        7d13s17
                  iadx=iox+((ix*(ix+1))/2)+in                           4d25s17
                  iadb=ihessbt+i3+iacto(isa)*(i1m+idoub(isao)*          7d13s17
     $                 (i4+noc(isbo)*i2m))                              9d19s17
                  bc(iadb)=bc(iadb)-4d0*bc(iadx)                        7d19s17
                 end do                                                  4d25s17
                end do                                                   4d25s17
               else                                                     4d25s17
                do i4=0,idoub(isbo)-1                                    4d25s17
                 do i3=0,iacto(isa)-1                                   7d13s17
                  i3p=i3+idoub(isa)                                     7d13s17
                  iadx=iox+i3p+noc(isa)*i4                              7d13s17
                  iadb=ihessbt+i3+iacto(isa)*(i1m+idoub(isao)*          7d13s17
     $                 (i4+noc(isbo)*i2m))                              9d19s17
                  bc(iadb)=bc(iadb)-4d0*bc(iadx)                        7d19s17
                 end do                                                  4d25s17
                end do                                                   4d25s17
               end if                                                   4d25s17
              end if                                                    4d25s17
              iox=iox+nrow                                              4d25s17
             end do                                                     4d25s17
             i10=1                                                      4d25s17
            end do                                                      4d25s17
           end if                                                        4d25s17
           if(isblk1(2,is).eq.isa.and.isblk1(1,is).eq.isbo.and.         7d13s17
     $         isblk1(3,is).eq.isao.and.isblk1(4,is).eq.isb.and.        7d13s17
     $          igot.eq.0)then                                          4d25s17
            igot=1                                                       4d25s17
            call ilimts(noc(isao),nvirtc(isb),mynprocg,mynowprog,il,ih, 7d13s17
     $           i1s,i1e,i2s,i2e)                                       4d25s17
            i10=i1s                                                     6d19s17
            i1n=noc(isao)                                               7d13s17
            iox=ionex(is)                                               4d25s17
            nrow=noc(isa)*noc(isbo)                                     7d13s17
            do i2=i2s,i2e                                               4d25s17
             i2m=i2-1                                                   4d25s17
             if(i2.eq.i2e)i1n=i1e                                       4d25s17
             do i1=i10,i1n                                              4d25s17
              i1m=i1-1                                                  7d13s17
              if(i1.le.idoub(isao))then                                 7d13s17
                call dumb(i1,i2)
               do i4=0,idoub(isbo)-1                                    4d25s17
                do i3=0,iacto(isa)-1                                    7d13s17
                 i3p=i3+idoub(isa)                                      7d13s17
                 iadx=iox+i4+noc(isbo)*i3p                              7d13s17
                 iadb=ihessbt+i3+iacto(isa)*(i1m+idoub(isao)*           7d13s17
     $                 (i4+noc(isbo)*i2m))                              9d19s17
                 bc(iadb)=bc(iadb)-4d0*bc(iadx)                         11d8s17
                end do                                                  4d25s17
               end do                                                   4d25s17
              end if                                                    4d25s17
              iox=iox+nrow                                              4d25s17
             end do                                                     4d25s17
             i10=1                                                      4d25s17
            end do                                                      4d25s17
           end if                                                        4d25s17
           end if                                                       9d25s17
c
c     aaaa: nd vn part a
c
           if(idoit(5).ne.0)then                                        9d25s17
           if(isblk1(3,is).eq.isao.and.isblk1(4,is).eq.isbv)then        9d25s17
            is1=isblk1(1,is)                                            9d25s17
            is2=isblk1(2,is)                                            9d25s17
            iswitch=iabs(is1-is2)                                       9d25s17
            iswitch=iswitch/max(1,iswitch)                              9d25s17
            call ilimts(noc(isao),nvirtc(isbv),mynprocg,mynowprog,il,ih,9d25s17
     $           i1s,i1e,i2s,i2e)                                       4d25s17
            if(is1.eq.is2)then                                          9d25s17
             nrow=(noc(is1)*(noc(is1)+1))/2                             9d25s17
            else                                                        9d25s17
             nrow=noc(is1)*noc(is2)                                     9d25s17
            end if                                                      9d25s17
            do isd=1,nsdlk                                              9d25s17
             if(isblk(1,isd).eq.isav.and.isblk(2,isd).eq.isbo.and.      9d25s17
     $            isblk(3,isd).eq.is1.and.jdenpt(isd).gt.0)then         9d25s17
              if(isav.eq.isbo)then                                      9d25s17
               fmult=-4d0                                               9d25s17
              else                                                      9d25s17
               fmult=-2d0                                               9d25s17
              end if                                                    9d25s17
              jswitch=iabs(isav-isbo)                                   9d25s17
              jswitch=jswitch/max(1,jswitch)                            9d25s17
              if(isav.eq.isbo)then                                      9d25s17
               nrowd=(iacto(isav)*(iacto(isav)+1))/2                    9d25s17
              else                                                      9d25s17
               nrowd=iacto(isav)*iacto(isbo)                            9d25s17
              end if                                                    9d25s17
              i10=i1s                                                     6d19s17
              i1n=noc(isao)                                               7d13s17
              iox=ionex(is)                                               4d25s17
              do i2=i2s,i2e                                             9d25s17
               if(i2.eq.i2e)i1n=i1e                                     9d25s17
               i2m=i2-1                                                 9d25s17
               do i1=i10,i1n                                            9d25s17
                if(i1.le.idoub(isao))then                                9d25s17
                 i1m=i1-1                                               9d25s17
                 jadh=ihessbt+iacto(isav)*(i1m+idoub(isao)*(idoub(isbo) 9d25s17
     $                +noc(isbo)*i2m))                                  9d25s17
                 do i4=0,iacto(is2)-1                                   9d25s17
                  i4p=i4+idoub(is2)                                     9d25s17
                  do i3=0,iacto(is1)-1                                  9d25s17
                   i3p=i3+idoub(is1)                                    9d25s17
                   ix=max(i3p,i4p)                                      9d25s17
                   in=min(i3p,i4p)                                      9d25s17
                   ieq=((ix*(ix+1))/2)+in                               9d25s17
                   inot=i3p+noc(is1)*i4p                                9d25s17
                   iadi=iox+(inot-ieq)*iswitch+ieq                      9d25s17
                   fact=fmult*bc(iadi)                                  9d25s17
                   ix=max(i3,i4)                                        9d25s17
                   in=min(i3,i4)                                        9d25s17
                   ieq=((ix*(ix+1))/2)+in                               9d25s17
                   inot=i3+iacto(is1)*i4                                9d25s17
                   icol=(inot-ieq)*iswitch+ieq                          9d25s17
                   iden=jdenpt(isd)+nrowd*icol                          9d25s17
                   do i7=0,iacto(isbo)-1                                9d25s17
                    iadh=jadh+nrowb*i7                                  9d25s17
                    do i6=0,iacto(isav)-1                               9d25s17
                     ix=max(i7,i6)                                      9d25s17
                     in=min(i7,i6)                                      9d25s17
                     ieq=((ix*(ix+1))/2)+in                             9d25s17
                     inot=i6+iacto(isav)*i7                             9d25s17
                     jden=iden+(inot-ieq)*jswitch+ieq                   9d25s17
                     bc(iadh+i6)=bc(iadh+i6)+fact*bc(jden)              9d25s17
                    end do                                              9d25s17
                   end do                                               9d25s17
                  end do                                                9d25s17
                 end do                                                 9d25s17
                end if                                                  9d25s17
                iox=iox+nrow                                            9d25s17
               end do                                                   9d25s17
               i10=1                                                    9d25s17
              end do                                                    9d25s17
             else if(isblk(2,isd).eq.isav.and.isblk(1,isd).eq.isbo.and. 9d25s17
     $            isblk(3,isd).eq.is1.and.jdenpt(isd).gt.0)then         9d25s17
              nrowd=iacto(isav)*iacto(isbo)                             9d25s17
              i10=i1s                                                     6d19s17
              i1n=noc(isao)                                               7d13s17
              iox=ionex(is)                                               4d25s17
              do i2=i2s,i2e                                             9d25s17
               if(i2.eq.i2e)i1n=i1e                                     9d25s17
               i2m=i2-1                                                 9d25s17
               do i1=i10,i1n                                            9d25s17
                if(i1.le.idoub(isao))then                                9d25s17
                 i1m=i1-1                                               9d25s17
                 jadh=ihessbt+iacto(isav)*(i1m+idoub(isao)*(idoub(isbo) 9d25s17
     $                +noc(isbo)*i2m))                                  9d25s17
                 do i4=0,iacto(is2)-1                                   9d25s17
                  i4p=i4+idoub(is2)                                     9d25s17
                  do i3=0,iacto(is1)-1                                  9d25s17
                   i3p=i3+idoub(is1)                                    9d25s17
                   iadi=iox+i3p+noc(is1)*i4p                            9d25s17
                   fact=-2d0*bc(iadi)                                   9d25s17
                   inot=i3+iacto(is1)*i4                                9d25s17
                   iden=jdenpt(isd)+nrowd*inot                          9d25s17
                   do i7=0,iacto(isbo)-1                                9d25s17
                    iadh=jadh+nrowb*i7                                  9d25s17
                    do i6=0,iacto(isav)-1                               9d25s17
                     jden=iden+i7+iacto(isbo)*i6                        9d25s17
                     bc(iadh+i6)=bc(iadh+i6)+fact*bc(jden)              9d25s17
                    end do                                              9d25s17
                   end do                                               9d25s17
                  end do                                                9d25s17
                 end do                                                 9d25s17
                end if                                                  9d25s17
                iox=iox+nrow                                            9d25s17
               end do                                                   9d25s17
               i10=1                                                    9d25s17
              end do                                                    9d25s17
             else if(isblk(2,isd).eq.isav.and.isblk(1,isd).eq.isbo.and. 9d25s17
     $            isblk(4,isd).eq.is1.and.jdenpt(isd).gt.0)then         9d25s17
              nrowd=iacto(isav)*iacto(isbo)                             9d25s17
              i10=i1s                                                     6d19s17
              i1n=noc(isao)                                               7d13s17
              iox=ionex(is)                                               4d25s17
              do i2=i2s,i2e                                             9d25s17
               if(i2.eq.i2e)i1n=i1e                                     9d25s17
               i2m=i2-1                                                 9d25s17
               do i1=i10,i1n                                            9d25s17
                if(i1.le.idoub(isao))then                                9d25s17
                 i1m=i1-1                                               9d25s17
                 jadh=ihessbt+iacto(isav)*(i1m+idoub(isao)*(idoub(isbo) 9d25s17
     $                +noc(isbo)*i2m))                                  9d25s17
                 do i4=0,iacto(is2)-1                                   9d25s17
                  i4p=i4+idoub(is2)                                     9d25s17
                  do i3=0,iacto(is1)-1                                  9d25s17
                   i3p=i3+idoub(is1)                                    9d25s17
                   iadi=iox+i3p+noc(is1)*i4p                            9d25s17
                   fact=-2d0*bc(iadi)                                   9d25s17
                   inot=i4+iacto(is2)*i3                                9d25s17
                   iden=jdenpt(isd)+nrowd*inot                          9d25s17
                   do i7=0,iacto(isbo)-1                                9d25s17
                    iadh=jadh+nrowb*i7                                  9d25s17
                    do i6=0,iacto(isav)-1                               9d25s17
                     jden=iden+i7+iacto(isbo)*i6                        9d25s17
                     bc(iadh+i6)=bc(iadh+i6)+fact*bc(jden)              9d25s17
                    end do                                              9d25s17
                   end do                                               9d25s17
                  end do                                                9d25s17
                 end do                                                 9d25s17
                end if                                                  9d25s17
                iox=iox+nrow                                            9d25s17
               end do                                                   9d25s17
               i10=1                                                    9d25s17
              end do                                                    9d25s17
             else if(isblk(1,isd).eq.isav.and.isblk(2,isd).eq.isbo.and. 9d25s17
     $            isblk(4,isd).eq.is1.and.jdenpt(isd).gt.0)then         9d25s17
              nrowd=iacto(isav)*iacto(isbo)                             9d25s17
              i10=i1s                                                     6d19s17
              i1n=noc(isao)                                               7d13s17
              iox=ionex(is)                                               4d25s17
              do i2=i2s,i2e                                             9d25s17
               if(i2.eq.i2e)i1n=i1e                                     9d25s17
               i2m=i2-1                                                 9d25s17
               do i1=i10,i1n                                            9d25s17
                if(i1.le.idoub(isao))then                                9d25s17
                 i1m=i1-1                                               9d25s17
                 jadh=ihessbt+iacto(isav)*(i1m+idoub(isao)*(idoub(isbo) 9d25s17
     $                +noc(isbo)*i2m))                                  9d25s17
                 do i4=0,iacto(is2)-1                                   9d25s17
                  i4p=i4+idoub(is2)                                     9d25s17
                  do i3=0,iacto(is1)-1                                  9d25s17
                   i3p=i3+idoub(is1)                                    9d25s17
                   iadi=iox+i3p+noc(is1)*i4p                            9d25s17
                   fact=-2d0*bc(iadi)                                   9d25s17
                   inot=i4+iacto(is2)*i3                                9d25s17
                   iden=jdenpt(isd)+nrowd*inot                          9d25s17
                   do i7=0,iacto(isbo)-1                                9d25s17
                    iadh=jadh+nrowb*i7                                  9d25s17
                    do i6=0,iacto(isav)-1                               9d25s17
                     jden=iden+i6+iacto(isav)*i7                        9d25s17
                     bc(iadh+i6)=bc(iadh+i6)+fact*bc(jden)              9d25s17
                    end do                                              9d25s17
                   end do                                               9d25s17
                  end do                                                9d25s17
                 end do                                                 9d25s17
                end if                                                  9d25s17
                iox=iox+nrow                                            9d25s17
               end do                                                   9d25s17
               i10=1                                                    9d25s17
              end do                                                    9d25s17
             end if                                                     9d25s17
            end do                                                      9d25s17
           end if                                                       9d25s17
c
c     aaaa nd vn part b
c
           if(isblk1(4,is).eq.isbv.and.isblk1(1,is).eq.isao)then        9d25s17
            is2=isblk1(2,is)                                            9d25s17
            is3=isblk1(3,is)                                            9d25s17
            iswitch=iabs(isao-is2)                                      9d25s17
            iswitch=iswitch/max(1,iswitch)                              9d25s17
            call ilimts(noc(is3),nvirtc(isbv),mynprocg,mynowprog,il,ih, 9d25s17
     $           i1s,i1e,i2s,i2e)                                       4d25s17
            if(isao.eq.is2)then                                         9d25s17
             nrow=(noc(isao)*(noc(isao)+1))/2                           9d25s17
            else                                                        9d25s17
             nrow=noc(isao)*noc(is2)                                    9d25s17
            end if                                                      9d25s17
            do isd=1,nsdlk                                              9d25s17
             if(isblk(1,isd).eq.isav.and.isblk(2,isd).eq.is2.and.       9d25s17
     $            isblk(3,isd).eq.isbo.and.jdenpt(isd).gt.0)then        9d25s17
              jswitch=iabs(isav-is2)                                    9d25s17
              jswitch=jswitch/max(1,jswitch)                            9d25s17
              kswitch=iabs(isbo-is3)                                    9d25s17
              kswitch=kswitch/max(1,kswitch)                            9d25s17
              if(isav.eq.is2)then                                       9d25s17
               fmult=-8d0                                               9d25s17
               nrowd=(iacto(isav)*(iacto(isav)+1))/2                    9d25s17
              else                                                      9d25s17
               fmult=-2d0                                               9d25s17
               nrowd=iacto(isav)*iacto(is2)                             9d25s17
              end if                                                    9d25s17
              i10=i1s                                                   9d25s17
              i1n=noc(is3)                                              9d25s17
              iox=ionex(is)                                             9d25s17
              do i2=i2s,i2e                                             9d25s17
               if(i2.eq.i2e)i1n=i1e                                     9d25s17
               i2m=i2-1                                                 9d25s17
               do i1=i10,i1n                                            9d25s17
                if(i1.gt.idoub(is3))then                                9d25s17
                 i1m=i1-1-idoub(is3)                                    9d25s17
                 do i4=0,iacto(is2)-1                                   9d25s17
                  i4p=i4+idoub(is2)                                     9d25s17
                  do i3=0,idoub(isao)-1                                 9d25s17
                   ix=max(i3,i4p)                                       9d25s17
                   in=min(i3,i4p)                                       9d25s17
                   ieq=((ix*(ix+1))/2)+in                               9d25s17
                   inot=i3+noc(isao)*i4p                                9d25s17
                   iadi=iox+(inot-ieq)*iswitch+ieq                      9d25s17
                   fact=fmult*bc(iadi)                                  9d25s17
                   iadh=ihessbt+iacto(isav)*(i3+idoub(isao)             9d25s17
     $                  *(idoub(isbo)+noc(isbo)*i2m))                   9d25s17
                   do i7=0,iacto(isbo)-1                                9d25s17
                    ix=max(i7,i1m)                                      9d25s17
                    in=min(i7,i1m)                                      9d25s17
                    ieq=((ix*(ix+1))/2)+in                              9d25s17
                    inot=i7+iacto(isbo)*i1m                             9d25s17
                    icol=(inot-ieq)*kswitch+ieq                         9d25s17
                    iden=jdenpt(isd)+nrowd*icol                         9d25s17
                    jadh=iadh+nrowb*i7                                  9d25s17
                    do i6=0,iacto(isav)-1                               9d25s17
                     ix=max(i6,i4)                                      9d25s17
                     in=min(i6,i4)                                      9d25s17
                     ieq=((ix*(ix+1))/2)+in                             9d25s17
                     inot=i6+iacto(isav)*i4                             9d25s17
                     jden=iden+(inot-ieq)*jswitch+ieq                   9d25s17
                     bc(jadh+i6)=bc(jadh+i6)+fact*bc(jden)              9d25s17
                    end do                                              9d25s17
                   end do                                               9d25s17
                  end do                                                9d25s17
                 end do                                                 9d25s17
                end if                                                  9d25s17
                iox=iox+nrow                                            9d25s17
               end do                                                   9d25s17
               i10=1                                                    9d25s17
              end do                                                    9d25s17
             else if(isblk(2,isd).eq.isav.and.isblk(1,isd).eq.is2.and.  9d25s17
     $            isblk(3,isd).eq.isbo.and.jdenpt(isd).gt.0)then        9d25s17
              fmult=-2d0                                                9d25s17
              nrowd=iacto(isav)*iacto(is2)                              9d25s17
              i10=i1s                                                   9d25s17
              i1n=noc(is3)                                              9d25s17
              iox=ionex(is)                                             9d25s17
              do i2=i2s,i2e                                             9d25s17
               if(i2.eq.i2e)i1n=i1e                                     9d25s17
               i2m=i2-1                                                 9d25s17
               do i1=i10,i1n                                            9d25s17
                if(i1.gt.idoub(is3))then                                9d25s17
                 i1m=i1-1-idoub(is3)                                    9d25s17
                 do i4=0,iacto(is2)-1                                   9d25s17
                  i4p=i4+idoub(is2)                                     9d25s17
                  do i3=0,idoub(isao)-1                                 9d25s17
                   ix=max(i3,i4p)                                       9d25s17
                   in=min(i3,i4p)                                       9d25s17
                   ieq=((ix*(ix+1))/2)+in                               9d25s17
                   inot=i3+noc(isao)*i4p                                9d25s17
                   iadi=iox+(inot-ieq)*iswitch+ieq                      9d25s17
                   fact=fmult*bc(iadi)                                  9d25s17
                   iadh=ihessbt+iacto(isav)*(i3+idoub(isao)             9d25s17
     $                  *(idoub(isbo)+noc(isbo)*i2m))                   9d25s17
                   do i7=0,iacto(isbo)-1                                9d25s17
                    inot=i7+iacto(isbo)*i1m                             9d25s17
                    iden=jdenpt(isd)+nrowd*inot                         9d25s17
                    jadh=iadh+nrowb*i7                                  9d25s17
                    do i6=0,iacto(isav)-1                               9d25s17
                     jden=iden+i4+iacto(is2)*i6                         9d25s17
                     bc(jadh+i6)=bc(jadh+i6)+fact*bc(jden)              9d25s17
                    end do                                              9d25s17
                   end do                                               9d25s17
                  end do                                                9d25s17
                 end do                                                 9d25s17
                end if                                                  9d25s17
                iox=iox+nrow                                            9d25s17
               end do                                                   9d25s17
               i10=1                                                    9d25s17
              end do                                                    9d25s17
             else if(isblk(2,isd).eq.isav.and.isblk(1,isd).eq.is2.and.  9d25s17
     $            isblk(4,isd).eq.isbo.and.jdenpt(isd).gt.0)then        9d25s17
              fmult=-2d0                                                9d25s17
              nrowd=iacto(isav)*iacto(is2)                              9d25s17
              i10=i1s                                                   9d25s17
              i1n=noc(is3)                                              9d25s17
              iox=ionex(is)                                             9d25s17
              do i2=i2s,i2e                                             9d25s17
               if(i2.eq.i2e)i1n=i1e                                     9d25s17
               i2m=i2-1                                                 9d25s17
               do i1=i10,i1n                                            9d25s17
                if(i1.gt.idoub(is3))then                                9d25s17
                 i1m=i1-1-idoub(is3)                                    9d25s17
                 do i4=0,iacto(is2)-1                                   9d25s17
                  i4p=i4+idoub(is2)                                     9d25s17
                  do i3=0,idoub(isao)-1                                 9d25s17
                   ix=max(i3,i4p)                                       9d25s17
                   in=min(i3,i4p)                                       9d25s17
                   ieq=((ix*(ix+1))/2)+in                               9d25s17
                   inot=i3+noc(isao)*i4p                                9d25s17
                   iadi=iox+(inot-ieq)*iswitch+ieq                      9d25s17
                   fact=fmult*bc(iadi)                                  9d25s17
                   iadh=ihessbt+iacto(isav)*(i3+idoub(isao)             9d25s17
     $                  *(idoub(isbo)+noc(isbo)*i2m))                   9d25s17
                   do i7=0,iacto(isbo)-1                                9d25s17
                    inot=i1m+iacto(is3)*i7                              9d25s17
                    iden=jdenpt(isd)+nrowd*inot                         9d25s17
                    jadh=iadh+nrowb*i7                                  9d25s17
                    do i6=0,iacto(isav)-1                               9d25s17
                     jden=iden+i4+iacto(is2)*i6                         9d25s17
                     bc(jadh+i6)=bc(jadh+i6)+fact*bc(jden)              9d25s17
                    end do                                              9d25s17
                   end do                                               9d25s17
                  end do                                                9d25s17
                 end do                                                 9d25s17
                end if                                                  9d25s17
                iox=iox+nrow                                            9d25s17
               end do                                                   9d25s17
               i10=1                                                    9d25s17
              end do                                                    9d25s17
             else if(isblk(1,isd).eq.isav.and.isblk(2,isd).eq.is2.and.  9d25s17
     $            isblk(4,isd).eq.isbo.and.jdenpt(isd).gt.0)then        9d25s17
              fmult=-2d0                                                9d25s17
              nrowd=iacto(isav)*iacto(is2)                              9d25s17
              i10=i1s                                                   9d25s17
              i1n=noc(is3)                                              9d25s17
              iox=ionex(is)                                             9d25s17
              do i2=i2s,i2e                                             9d25s17
               if(i2.eq.i2e)i1n=i1e                                     9d25s17
               i2m=i2-1                                                 9d25s17
               do i1=i10,i1n                                            9d25s17
                if(i1.gt.idoub(is3))then                                9d25s17
                 i1m=i1-1-idoub(is3)                                    9d25s17
                 do i4=0,iacto(is2)-1                                   9d25s17
                  i4p=i4+idoub(is2)                                     9d25s17
                  do i3=0,idoub(isao)-1                                 9d25s17
                   ix=max(i3,i4p)                                       9d25s17
                   in=min(i3,i4p)                                       9d25s17
                   ieq=((ix*(ix+1))/2)+in                               9d25s17
                   inot=i3+noc(isao)*i4p                                9d25s17
                   iadi=iox+(inot-ieq)*iswitch+ieq                      9d25s17
                   fact=fmult*bc(iadi)                                  9d25s17
                   iadh=ihessbt+iacto(isav)*(i3+idoub(isao)             9d25s17
     $                  *(idoub(isbo)+noc(isbo)*i2m))                   9d25s17
                   do i7=0,iacto(isbo)-1                                9d25s17
                    inot=i1m+iacto(is3)*i7                              9d25s17
                    iden=jdenpt(isd)+nrowd*inot                         9d25s17
                    jadh=iadh+nrowb*i7                                  9d25s17
                    do i6=0,iacto(isav)-1                               9d25s17
                     jden=iden+i6+iacto(isav)*i4                        9d25s17
                     bc(jadh+i6)=bc(jadh+i6)+fact*bc(jden)              9d25s17
                    end do                                              9d25s17
                   end do                                               9d25s17
                  end do                                                9d25s17
                 end do                                                 9d25s17
                end if                                                  9d25s17
                iox=iox+nrow                                            9d25s17
               end do                                                   9d25s17
               i10=1                                                    9d25s17
              end do                                                    9d25s17
             end if                                                     9d25s17
            end do                                                      9d25s17
           else if(isblk1(4,is).eq.isbv.and.isblk1(2,is).eq.isao)then   9d25s17
            is1=isblk1(1,is)                                            9d25s17
            is3=isblk1(3,is)                                            9d25s17
            call ilimts(noc(is3),nvirtc(isbv),mynprocg,mynowprog,il,ih, 9d25s17
     $           i1s,i1e,i2s,i2e)                                       4d25s17
            nrow=noc(isao)*noc(is1)                                     9d25s17
            do isd=1,nsdlk                                              9d25s17
             if(isblk(1,isd).eq.isav.and.isblk(2,isd).eq.is1.and.       9d25s17
     $            isblk(3,isd).eq.isbo.and.jdenpt(isd).gt.0)then        9d25s17
              jswitch=iabs(isav-is1)                                    9d25s17
              jswitch=jswitch/max(1,jswitch)                            9d25s17
              kswitch=iabs(isbo-is3)                                    9d25s17
              kswitch=kswitch/max(1,kswitch)                            9d25s17
              if(isav.eq.is1)then                                       9d25s17
               fmult=-8d0                                               9d25s17
               nrowd=(iacto(isav)*(iacto(isav)+1))/2                    9d25s17
              else                                                      9d25s17
               fmult=-2d0                                               9d25s17
               nrowd=iacto(isav)*iacto(is1)                             9d25s17
              end if                                                    9d25s17
              i10=i1s                                                   9d25s17
              i1n=noc(is3)                                              9d25s17
              iox=ionex(is)                                             9d25s17
              do i2=i2s,i2e                                             9d25s17
               if(i2.eq.i2e)i1n=i1e                                     9d25s17
               i2m=i2-1                                                 9d25s17
               do i1=i10,i1n                                            9d25s17
                if(i1.gt.idoub(is3))then                                9d25s17
                 i1m=i1-1-idoub(is3)                                    9d25s17
                 do i4=0,idoub(isao)-1                                   9d25s17
                  do i3=0,iacto(is1)-1                                  9d25s17
                   i3p=i3+idoub(is1)                                    9d25s17
                   iadi=iox+i3p+noc(is1)*i4                             9d25s17
                   fact=fmult*bc(iadi)                                  9d25s17
                   iadh=ihessbt+iacto(isav)*(i4+idoub(isao)             9d25s17
     $                  *(idoub(isbo)+noc(isbo)*i2m))                   9d25s17
                   do i7=0,iacto(isbo)-1                                9d25s17
                    ix=max(i7,i1m)                                      9d25s17
                    in=min(i7,i1m)                                      9d25s17
                    ieq=((ix*(ix+1))/2)+in                              9d25s17
                    inot=i7+iacto(isbo)*i1m                             9d25s17
                    icol=(inot-ieq)*kswitch+ieq                         9d25s17
                    iden=jdenpt(isd)+nrowd*icol                         9d25s17
                    jadh=iadh+nrowb*i7                                  9d25s17
                    do i6=0,iacto(isav)-1                               9d25s17
                     ix=max(i6,i3)                                      9d25s17
                     in=min(i6,i3)                                      9d25s17
                     ieq=((ix*(ix+1))/2)+in                             9d25s17
                     inot=i6+iacto(isav)*i3                             9d25s17
                     jden=iden+(inot-ieq)*jswitch+ieq                   9d25s17
                     bc(jadh+i6)=bc(jadh+i6)+fact*bc(jden)              9d25s17
                    end do                                              9d25s17
                   end do                                               9d25s17
                  end do                                                9d25s17
                 end do                                                 9d25s17
                end if                                                  9d25s17
                iox=iox+nrow                                            9d25s17
               end do                                                   9d25s17
               i10=1                                                    9d25s17
              end do                                                    9d25s17
             else if(isblk(2,isd).eq.isav.and.isblk(1,isd).eq.is1.and.  9d25s17
     $            isblk(3,isd).eq.isbo.and.jdenpt(isd).gt.0)then        9d25s17
              fmult=-2d0                                                9d25s17
              nrowd=iacto(isav)*iacto(is1)                              9d25s17
              i10=i1s                                                   9d25s17
              i1n=noc(is3)                                              9d25s17
              iox=ionex(is)                                             9d25s17
              do i2=i2s,i2e                                             9d25s17
               if(i2.eq.i2e)i1n=i1e                                     9d25s17
               i2m=i2-1                                                 9d25s17
               do i1=i10,i1n                                            9d25s17
                if(i1.gt.idoub(is3))then                                9d25s17
                 i1m=i1-1-idoub(is3)                                    9d25s17
                 do i4=0,idoub(isao)-1                                  9d25s17
                  do i3=0,iacto(is1)-1                                  9d25s17
                   i3p=i3+idoub(is1)                                    9d25s17
                   iadi=iox+i3p+noc(is1)*i4                             9d25s17
                   fact=fmult*bc(iadi)                                  9d25s17
                   iadh=ihessbt+iacto(isav)*(i4+idoub(isao)             9d25s17
     $                  *(idoub(isbo)+noc(isbo)*i2m))                   9d25s17
                   do i7=0,iacto(isbo)-1                                9d25s17
                    inot=i7+iacto(isbo)*i1m                             9d25s17
                    iden=jdenpt(isd)+nrowd*inot                         9d25s17
                    jadh=iadh+nrowb*i7                                  9d25s17
                    do i6=0,iacto(isav)-1                               9d25s17
                     jden=iden+i3+iacto(is1)*i6                         9d25s17
                     bc(jadh+i6)=bc(jadh+i6)+fact*bc(jden)              9d25s17
                    end do                                              9d25s17
                   end do                                               9d25s17
                  end do                                                9d25s17
                 end do                                                 9d25s17
                end if                                                  9d25s17
                iox=iox+nrow                                            9d25s17
               end do                                                   9d25s17
               i10=1                                                    9d25s17
              end do                                                    9d25s17
             else if(isblk(2,isd).eq.isav.and.isblk(1,isd).eq.is1.and.  9d25s17
     $            isblk(4,isd).eq.isbo.and.jdenpt(isd).gt.0)then        9d25s17
              fmult=-2d0                                                9d25s17
              nrowd=iacto(isav)*iacto(is1)                              9d25s17
              i10=i1s                                                   9d25s17
              i1n=noc(is3)                                              9d25s17
              iox=ionex(is)                                             9d25s17
              do i2=i2s,i2e                                             9d25s17
               if(i2.eq.i2e)i1n=i1e                                     9d25s17
               i2m=i2-1                                                 9d25s17
               do i1=i10,i1n                                            9d25s17
                if(i1.gt.idoub(is3))then                                9d25s17
                 i1m=i1-1-idoub(is3)                                    9d25s17
                 do i4=0,idoub(isao)-1                                  9d25s17
                  do i3=0,iacto(is1)-1                                  9d25s17
                   i3p=i3+idoub(is1)                                    9d25s17
                   iadi=iox+i3p+noc(is1)*i4                             9d25s17
                   fact=fmult*bc(iadi)                                  9d25s17
                   iadh=ihessbt+iacto(isav)*(i4+idoub(isao)             9d25s17
     $                  *(idoub(isbo)+noc(isbo)*i2m))                   9d25s17
                   do i7=0,iacto(isbo)-1                                9d25s17
                    inot=i1m+iacto(is3)*i7                              9d25s17
                    iden=jdenpt(isd)+nrowd*inot                         9d25s17
                    jadh=iadh+nrowb*i7                                  9d25s17
                    do i6=0,iacto(isav)-1                               9d25s17
                     jden=iden+i3+iacto(is1)*i6                         9d25s17
                     bc(jadh+i6)=bc(jadh+i6)+fact*bc(jden)              9d25s17
                    end do                                              9d25s17
                   end do                                               9d25s17
                  end do                                                9d25s17
                 end do                                                 9d25s17
                end if                                                  9d25s17
                iox=iox+nrow                                            9d25s17
               end do                                                   9d25s17
               i10=1                                                    9d25s17
              end do                                                    9d25s17
             else if(isblk(1,isd).eq.isav.and.isblk(2,isd).eq.is1.and.  9d25s17
     $            isblk(4,isd).eq.isbo.and.jdenpt(isd).gt.0)then        9d25s17
              fmult=-2d0                                                9d25s17
              nrowd=iacto(isav)*iacto(is1)                              9d25s17
              i10=i1s                                                   9d25s17
              i1n=noc(is3)                                              9d25s17
              iox=ionex(is)                                             9d25s17
              do i2=i2s,i2e                                             9d25s17
               if(i2.eq.i2e)i1n=i1e                                     9d25s17
               i2m=i2-1                                                 9d25s17
               do i1=i10,i1n                                            9d25s17
                if(i1.gt.idoub(is3))then                                9d25s17
                 i1m=i1-1-idoub(is3)                                    9d25s17
                 do i4=0,idoub(isao)-1                                  9d25s17
                  do i3=0,iacto(is1)-1                                  9d25s17
                   i3p=i3+idoub(is1)                                    9d25s17
                   iadi=iox+i3p+noc(is1)*i4
                   fact=fmult*bc(iadi)                                  9d25s17
                   iadh=ihessbt+iacto(isav)*(i4+idoub(isao)             9d25s17
     $                  *(idoub(isbo)+noc(isbo)*i2m))                   9d25s17
                   do i7=0,iacto(isbo)-1                                9d25s17
                    inot=i1m+iacto(is3)*i7                              9d25s17
                    iden=jdenpt(isd)+nrowd*inot                         9d25s17
                    jadh=iadh+nrowb*i7                                  9d25s17
                    do i6=0,iacto(isav)-1                               9d25s17
                     jden=iden+i6+iacto(isav)*i3                        9d25s17
                     bc(jadh+i6)=bc(jadh+i6)+fact*bc(jden)              9d25s17
                    end do                                              9d25s17
                   end do                                               9d25s17
                  end do                                                9d25s17
                 end do                                                 9d25s17
                end if                                                  9d25s17
                iox=iox+nrow                                            9d25s17
               end do                                                   9d25s17
               i10=1                                                    9d25s17
              end do                                                    9d25s17
             end if                                                     9d25s17
            end do                                                      9d25s17
           end if                                                       9d25s17
           end if                                                       9d25s17
          end do                                                         4d25s17
          call dws_gsumf(bc(ihessbt),nwds)                              4d25s17
          i10=i1sb                                                      4d25s17
          i1n=noc(isbo)                                                 4d25s17
          ihb=ihessb(isa,isb)                                           4d25s17
          do i2=i2sb,i2eb                                               4d25s17
           i2m=i2-1                                                     4d25s17
           if(i2.eq.i2eb)i1n=i1eb                                       4d25s17
           do i1=i10,i1n                                                4d25s17
            i1m=i1-1                                                    4d25s17
             iadt=ihessbt+nrowb*(i1m+noc(isbo)*i2m)                     9d19s17
c
c     dddd
c
             do i34=0,idoub(isao)*iacto(isa)-1                          4d25s17
              bc(ihb+i34)=bc(ihb+i34)+bc(iadt+i34)                      5d4s17
             end do                                                     4d25s17
            ihb=ihb+nrowb                                                4d25s17
           end do                                                       4d25s17
           i10=1                                                        4d25s17
          end do                                                        4d25s17
          ibcoff=ihessbt                                                4d25s17
         end if                                                         4d25s17
         do is=1,nsdlk
          if(min(idoub(isao),idoub(isbo),noc(isb),noc(isa)).ne.0        7d11s22
     $         .and.lsym)then                                           7d11s22
c
c     for rotations between double and active orbitals
c
c     dddd for shift in energy due to doubly occupied orbials                7d20s17
c     K part
c
           if(idoit(3).ne.0)then                                        9d25s17
           igot=0                                                       4d20s17
           if(isblk(1,is).eq.isa.and.isblk(2,is).eq.isao.and.
     $          isblk(3,is).eq.isb.and.isblk(4,is).eq.isbo)then          4d20s17
            igot=1                                                      4d20s17
            call ilimts(noc(isb),noc(isbo),mynprocg,mynowprog,ilo,iho,    4d20s17
     $         i1so,i1eo,i2so,i2eo)                                     4d20s17
            if(isa.eq.isao)then
             nrowo=(noc(isa)*(noc(isa)+1))/2                             4d20s17
            else                                                         4d20s17
             nrowo=noc(isa)*noc(isao)                                    4d20s17
            end if                                                       4d20s17
            i10=i1so                                                     4d20s17
            i1n=noc(isb)                                                 4d20s17
            i4o=ioooo(is)                                                4d20s17
            do i2=i2so,i2eo                                                4d20s17
             i2m=i2-1                                                   4d20s17
             if(i2.eq.i2eo)i1n=i1eo                                       4d20s17
             do i1=i10,i1n                                               4d20s17
              i1m=i1-1-idoub(isb)                                       4d20s17
              if(i2.le.idoub(isbo).and.i1.gt.idoub(isb))then            4d20s17
               if(isa.eq.isao)then
                do i4=0,idoub(isao)-1                                   4d20s17
                 do i3=0,iacto(isa)-1                                   4d20s17
                  i3p=i3+idoub(isa)                                     4d20s17
                  ix=max(i4,i3p)                                        4d20s17
                  in=min(i4,i3p)                                        4d20s17
                  iado=i4o+((ix*(ix+1))/2)+in                           4d20s17
                  iadh=ihessa(isa,isb)+i2m+idoub(isbo)*(i4+idoub(isao)*  4d20s17
     $                 (i3+iacto(isa)*i1m))                             4d20s17
                  bc(iadh)=bc(iadh)+16d0*bc(iado)                       4d20s17
                 end do                                                 4d20s17
                end do                                                  4d20s17
               else
                do i4=0,idoub(isao)-1                                   4d20s17
                 do i3=0,iacto(isa)-1                                   4d20s17
                  i3p=i3+idoub(isa)                                     4d20s17
                  iado=i4o+i3p+noc(isa)*i4                              4d20s17
                  iadh=ihessa(isa,isb)+i2m+idoub(isbo)*(i4+idoub(isao)*  4d20s17
     $                 (i3+iacto(isa)*i1m))                             4d20s17
                  bc(iadh)=bc(iadh)+16d0*bc(iado)                       4d20s17
                 end do                                                 4d20s17
                end do                                                  4d20s17
               end if
              end if                                                    4d20s17
              i4o=i4o+nrowo                                              4d20s17
             end do                                                     4d25s17
             i10=1                                                       4d20s17
            end do                                                       4d20s17
           else if(isa.ne.isao.and.isblk(1,is).eq.isao.and.isblk(2,is). 10d30s17
     $        eq.isa.and.isblk(3,is).eq.isb.and.isblk(4,is).eq.isbo)then10d30s17
            igot=1                                                      10d30s17
            call ilimts(noc(isb),noc(isbo),mynprocg,mynowprog,ilo,iho,  10d30s17
     $         i1so,i1eo,i2so,i2eo)                                     10d30s17
            nrowo=noc(isa)*noc(isao)                                    10d30s17
            i10=i1so                                                    10d30s17
            i1n=noc(isb)                                                10d30s17
            i4o=ioooo(is)                                               10d30s17
            do i2=i2so,i2eo                                             10d30s17
             i2m=i2-1                                                   10d30s17
             if(i2.eq.i2eo)i1n=i1eo                                     10d30s17
             do i1=i10,i1n                                              10d30s17
              i1m=i1-1-idoub(isb)                                       10d30s17
              if(i2.le.idoub(isbo).and.i1.gt.idoub(isb))then            10d30s17
               do i4=0,idoub(isao)-1                                    10d30s17
                do i3=0,iacto(isa)-1                                    10d30s17
                 i3p=i3+idoub(isa)                                      10d30s17
                 iado=i4o+i4+noc(isao)*i3p                              10d30s17
                 iadh=ihessa(isa,isb)+i2m+idoub(isbo)*(i4+idoub(isao)*  4d20s17
     $                (i3+iacto(isa)*i1m))                              10d30s17
                 bc(iadh)=bc(iadh)+16d0*bc(iado)                        10d30s17
                end do                                                  10d30s17
               end do                                                   10d30s17
              end if                                                    10d30s17
              i4o=i4o+nrowo                                             10d30s17
             end do                                                     10d30s17
             i10=1                                                      10d30s17
            end do                                                      10d30s17
           end if                                                        4d20s17
           if(igot.eq.0.and.isblk(2,is).eq.isa.and.isblk(1,is).eq.isao
     $          .and.isblk(4,is).eq.isb.and.isblk(3,is).eq.isbo)then    4d20s17
            call ilimts(noc(isbo),noc(isb),mynprocg,mynowprog,ilo,iho,  4d20s17
     $         i1so,i1eo,i2so,i2eo)                                     4d20s17
            nrowo=noc(isa)*noc(isao)                                    4d20s17
            i10=i1so                                                     4d20s17
            i1n=noc(isbo)                                               4d20s17
            i4o=ioooo(is)                                                4d20s17
            do i2=i2so,i2eo                                                4d20s17
             i2m=i2-1-idoub(isb)                                        4d20s17
             if(i2.eq.i2eo)i1n=i1eo                                       4d20s17
             do i1=i10,i1n                                               4d20s17
              i1m=i1-1                                                  4d20s17
              if(i1.le.idoub(isbo).and.i2.gt.idoub(isb))then            4d20s17
               do i4=0,iacto(isa)-1                                     4d20s17
                i4p=i4+idoub(isa)                                       4d20s17
                do i3=0,idoub(isao)-1                                   4d20s17
                 iado=i4o+i3+noc(isao)*i4p                              4d20s17
                 iadh=ihessa(isa,isb)+i1m+idoub(isbo)*(i3+idoub(isao)*  4d20s17
     $                (i4+iacto(isa)*i2m))                              4d20s17
                 bc(iadh)=bc(iadh)+16d0*bc(iado)                        4d20s17
                end do                                                  4d20s17
               end do                                                   4d20s17
              end if                                                    4d20s17
              i4o=i4o+nrowo                                              4d20s17
             end do                                                     4d25s17
             i10=1                                                       4d20s17
            end do                                                       4d20s17
           end if                                                        4d20s17
c
c     dddd pseudo K part
c
           igot=0                                                       4d20s17
           if(isblk(1,is).eq.isa.and.isblk(4,is).eq.isao.and.           5d2s17
     $          isblk(3,is).eq.isb.and.isblk(2,is).eq.isbo)then         5d2s17
            igot=1                                                      4d20s17
            call ilimts(noc(isb),noc(isao),mynprocg,mynowprog,ilo,iho,  5d2s17
     $         i1so,i1eo,i2so,i2eo)                                     4d20s17
            if(isa.eq.isbo)then                                         5d2s17
             nrowo=(noc(isa)*(noc(isa)+1))/2                             4d20s17
            else                                                         4d20s17
             nrowo=noc(isa)*noc(isbo)                                   5d2s17
            end if                                                       4d20s17
            i10=i1so                                                     4d20s17
            i1n=noc(isb)                                                 4d20s17
            i4o=ioooo(is)                                                4d20s17
            do i2=i2so,i2eo                                                4d20s17
             i2m=i2-1                                                   4d20s17
             if(i2.eq.i2eo)i1n=i1eo                                       4d20s17
             do i1=i10,i1n                                               4d20s17
              i1m=i1-1-idoub(isb)                                       4d20s17
              if(i2.le.idoub(isao).and.i1.gt.idoub(isb))then            5d2s17
               if(isa.eq.isbo)then                                      5d2s17
                do i4=0,idoub(isbo)-1                                   5d2s17
                 do i3=0,iacto(isa)-1                                   4d20s17
                  i3p=i3+idoub(isa)                                     4d20s17
                  ix=max(i4,i3p)                                        4d20s17
                  in=min(i4,i3p)                                        4d20s17
                  iado=i4o+((ix*(ix+1))/2)+in                           4d20s17
                  iadh=ihessa(isa,isb)+i4+idoub(isbo)*(i2m+idoub(isao)* 5d2s17
     $                 (i3+iacto(isa)*i1m))                             4d20s17
                  bc(iadh)=bc(iadh)-4d0*bc(iado)                        5d2s17
                 end do                                                 4d20s17
                end do                                                  4d20s17
               else
                do i4=0,idoub(isbo)-1                                   5d2s17
                 do i3=0,iacto(isa)-1                                   4d20s17
                  i3p=i3+idoub(isa)                                     4d20s17
                  iado=i4o+i3p+noc(isa)*i4                              4d20s17
                  iadh=ihessa(isa,isb)+i4+idoub(isbo)*(i2m+idoub(isao)* 5d2s17
     $                 (i3+iacto(isa)*i1m))                             4d20s17
                  bc(iadh)=bc(iadh)-4d0*bc(iado)                        5d2s17
                 end do                                                 4d20s17
                end do                                                  4d20s17
               end if
              end if                                                    4d20s17
              i4o=i4o+nrowo                                              4d20s17
             end do                                                     4d25s17
             i10=1                                                       4d20s17
            end do                                                       4d20s17
           else if(isa.ne.isbo.and.isblk(1,is).eq.isbo.and.isblk(4,is). 10d30s17
     $           eq.isao.and.isblk(3,is).eq.isb.and.isblk(2,is).eq.isa) 10d30s17
     $           then                                                   10d30s17
            igot=1                                                      10d30s17
            call ilimts(noc(isb),noc(isao),mynprocg,mynowprog,ilo,iho,  10d30s17
     $         i1so,i1eo,i2so,i2eo)                                     10d30s17
            nrowo=noc(isa)*noc(isbo)                                    10d30s17
            i10=i1so                                                    10d30s17
            i1n=noc(isb)                                                10d30s17
            i4o=ioooo(is)                                               10d30s17
            do i2=i2so,i2eo                                             10d30s17
             i2m=i2-1                                                   10d30s17
             if(i2.eq.i2eo)i1n=i1eo                                     10d30s17
             do i1=i10,i1n                                              10d30s17
              i1m=i1-1-idoub(isb)                                       10d30s17
              if(i2.le.idoub(isao).and.i1.gt.idoub(isb))then            10d30s17
               do i4=0,idoub(isbo)-1                                    10d30s17
                do i3=0,iacto(isa)-1                                    10d30s17
                 i3p=i3+idoub(isa)                                      10d30s17
                 iado=i4o+i4+noc(isbo)*i3p                              10d30s17
                 iadh=ihessa(isa,isb)+i4+idoub(isbo)*(i2m+idoub(isao)*  10d30s17
     $                (i3+iacto(isa)*i1m))                              10d30s17
                 bc(iadh)=bc(iadh)-4d0*bc(iado)                         10d30s17
                end do                                                  10d30s17
               end do                                                   10d30s17
              end if                                                    10d30s17
              i4o=i4o+nrowo                                             10d30s17
             end do                                                     10d30s17
             i10=1                                                      10d30s17
            end do                                                      10d30s17
           end if                                                       10d30s17
           if(igot.eq.0.and.isblk(2,is).eq.isa.and.isblk(1,is).eq.isbo  5d2s17
     $          .and.isblk(4,is).eq.isb.and.isblk(3,is).eq.isao)then    5d2s17
            call ilimts(noc(isao),noc(isb),mynprocg,mynowprog,ilo,iho,  5d2s17
     $         i1so,i1eo,i2so,i2eo)                                     4d20s17
            nrowo=noc(isa)*noc(isbo)                                    5d2s17
            i10=i1so                                                     4d20s17
            i1n=noc(isao)                                               5d2s17
            i4o=ioooo(is)                                                4d20s17
            do i2=i2so,i2eo                                                4d20s17
             i2m=i2-1-idoub(isb)                                        4d20s17
             if(i2.eq.i2eo)i1n=i1eo                                       4d20s17
             do i1=i10,i1n                                               4d20s17
              i1m=i1-1                                                  4d20s17
              if(i1.le.idoub(isao).and.i2.gt.idoub(isb))then            5d2s17
               do i4=0,iacto(isa)-1                                     4d20s17
                i4p=i4+idoub(isa)                                       4d20s17
                do i3=0,idoub(isbo)-1                                   5d2s17
                 iado=i4o+i3+noc(isbo)*i4p                              5d2s17
                 iadh=ihessa(isa,isb)+i3+idoub(isbo)*(i1m+idoub(isao)*  5d2s17
     $                (i4+iacto(isa)*i2m))                              4d20s17
                 bc(iadh)=bc(iadh)-4d0*bc(iado)                         5d2s17
                end do                                                  4d20s17
               end do                                                   4d20s17
              end if                                                    4d20s17
              i4o=i4o+nrowo                                              4d20s17
             end do                                                     4d25s17
             i10=1                                                       4d20s17
            end do                                                       4d20s17
           end if                                                        4d20s17
c
c     dddd J part
c
           igot=0
           if(isblk(1,is).eq.isa.and.isblk(2,is).eq.isb.and.            4d20s17
     $          isblk(3,is).eq.isao.and.isblk(4,is).eq.isbo)then        4d20s17
            igot=1                                                       4d20s17
            call ilimts(noc(isao),noc(isbo),mynprocg,mynowprog,ilo,iho,  4d20s17
     $           i1so,i1eo,i2so,i2eo)                                    4d20s17
            i10=i1so                                                    4d20s17
            i1n=noc(isao)                                               4d20s17
            if(isa.eq.isb)then
             nrowo=(noc(isa)*(noc(isa)+1))/2                             4d20s17
            else                                                        4d20s17
             nrowo=noc(isa)*noc(isb)                                    4d20s17
            end if                                                      4d20s17
            i4o=ioooo(is)                                               4d20s17
            do i2=i2so,i2eo
             i2m=i2-1                                                   4d20s17
             if(i2.eq.i2eo)i1n=i1eo
             do i1=i10,i1n
              i1m=i1-1                                                  4d20s17
              if(i1.le.idoub(isao).and.i2.le.idoub(isbo))then           4d20s17
               if(isa.eq.isb)then                                       4d20s17
                do i4=0,iacto(isb)-1                                     4d20s17
                 i4p=i4+idoub(isb)                                      4d20s17
                 do i3=0,iacto(isa)-1                                    4d20s17
                  i3p=i3+idoub(isa)                                     4d20s17
                  ix=max(i4p,i3p)                                       4d20s17
                  in=min(i4p,i3p)                                       4d20s17
                  iado=i4o+((ix*(ix+1))/2)+in                           4d20s17
                  iadh=ihessa(isa,isb)+i2m+idoub(isbo)*(i1m+idoub(isao)*4d20s17
     $                 (i3+iacto(isa)*i4))                              4d20s17
                  bc(iadh)=bc(iadh)-4d0*bc(iado)                        4d20s17
                 end do                                                  4d20s17
                end do                                                   4d20s17
               else                                                     4d20s17
                do i4=0,iacto(isb)-1                                     4d20s17
                 i4p=i4+idoub(isb)                                      4d20s17
                 do i3=0,iacto(isa)-1                                    4d20s17
                  i3p=i3+idoub(isa)                                     4d20s17
                  iado=i4o+i3p+noc(isa)*i4p                             4d20s17
                  iadh=ihessa(isa,isb)+i2m+idoub(isbo)*(i1m+idoub(isao)*4d20s17
     $                 (i3+iacto(isa)*i4))                              4d20s17
c     my best guess
                  bc(iadh)=bc(iadh)-4d0*bc(iado)                        10d30s17
                 end do                                                  4d20s17
                end do                                                   4d20s17
               end if                                                   4d20s17
              end if                                                    4d20s17
              i4o=i4o+nrowo
             end do
             i10=1
            end do
           end if                                                       4d20s17
           if(igot.eq.0.and.isblk(2,is).eq.isa.and.isblk(1,is).eq.isb   4d20s17
     $          .and.isblk(3,is).eq.isao.and.isblk(4,is).eq.isbo)then   4d20s17
            igot=1
            call ilimts(noc(isao),noc(isbo),mynprocg,mynowprog,ilo,iho,  4d20s17
     $           i1so,i1eo,i2so,i2eo)                                    4d20s17
            i10=i1so                                                    4d20s17
            i1n=noc(isao)                                               4d20s17
            nrowo=noc(isa)*noc(isb)                                     4d20s17
            i4o=ioooo(is)                                               4d20s17
            do i2=i2so,i2eo
             i2m=i2-1                                                   4d20s17
             if(i2.eq.i2eo)i1n=i1eo
             do i1=i10,i1n
              i1m=i1-1                                                  4d20s17
              if(i1.le.idoub(isao).and.i2.le.idoub(isbo))then           4d20s17
                do i4=0,iacto(isa)-1                                     4d20s17
                 i4p=i4+idoub(isa)                                      4d20s17
                 do i3=0,iacto(isb)-1                                    4d20s17
                  i3p=i3+idoub(isb)                                     4d20s17
                  iado=i4o+i3p+noc(isb)*i4p                             4d20s17
                  iadh=ihessa(isa,isb)+i2m+idoub(isbo)*(i1m+idoub(isao)*4d20s17
     $                 (i4+iacto(isa)*i3))                              4d20s17
c     my best guess
                  bc(iadh)=bc(iadh)-4d0*bc(iado)                        10d30s17
                 end do                                                  4d20s17
                end do                                                   4d20s17
              end if                                                    4d20s17
              i4o=i4o+nrowo
             end do
             i10=1
            end do
           end if                                                       4d20s17
           if(igot.eq.0.and.isblk(2,is).eq.isa.and.isblk(1,is).eq.isb   4d20s17
     $          .and.isblk(4,is).eq.isao.and.isblk(3,is).eq.isbo)then   4d20s17
            igot=1
            call ilimts(noc(isbo),noc(isao),mynprocg,mynowprog,ilo,iho,  4d20s17
     $           i1so,i1eo,i2so,i2eo)                                    4d20s17
            i10=i1so                                                    4d20s17
            i1n=noc(isbo)                                               4d20s17
            nrowo=noc(isa)*noc(isb)                                     4d20s17
            i4o=ioooo(is)                                               4d20s17
            do i2=i2so,i2eo
             i2m=i2-1                                                   4d20s17
             if(i2.eq.i2eo)i1n=i1eo
             do i1=i10,i1n
              i1m=i1-1                                                  4d20s17
              if(i1.le.idoub(isbo).and.i2.le.idoub(isao))then           4d20s17
                do i4=0,iacto(isa)-1                                     4d20s17
                 i4p=i4+idoub(isa)                                      4d20s17
                 do i3=0,iacto(isb)-1                                    4d20s17
                  i3p=i3+idoub(isb)                                     4d20s17
                  iado=i4o+i3p+noc(isb)*i4p                             4d20s17
                  iadh=ihessa(isa,isb)+i1m+idoub(isbo)*(i2m+idoub(isao)*4d20s17
     $                 (i4+iacto(isa)*i3))                              4d20s17
                  bc(iadh)=bc(iadh)-4d0*bc(iado)                        10d30s17
                 end do                                                  4d20s17
                end do                                                   4d20s17
              end if                                                    4d20s17
              i4o=i4o+nrowo
             end do
             i10=1
            end do
           end if                                                       4d20s17
           if(igot.eq.0.and.isblk(1,is).eq.isa.and.isblk(2,is).eq.isb   4d20s17
     $          .and.isblk(4,is).eq.isao.and.isblk(3,is).eq.isbo)then   4d20s17
            igot=1
            call ilimts(noc(isbo),noc(isao),mynprocg,mynowprog,ilo,iho,  4d20s17
     $           i1so,i1eo,i2so,i2eo)                                    4d20s17
            i10=i1so                                                    4d20s17
            i1n=noc(isbo)                                               4d20s17
            nrowo=noc(isa)*noc(isb)                                     4d20s17
            i4o=ioooo(is)                                               4d20s17
            do i2=i2so,i2eo
             i2m=i2-1                                                   4d20s17
             if(i2.eq.i2eo)i1n=i1eo
             do i1=i10,i1n
              i1m=i1-1                                                  4d20s17
              if(i1.le.idoub(isbo).and.i2.le.idoub(isao))then           4d20s17
                do i4=0,iacto(isb)-1                                     4d20s17
                 i4p=i4+idoub(isb)                                      4d20s17
                 do i3=0,iacto(isa)-1                                    4d20s17
                  i3p=i3+idoub(isa)                                     4d20s17
                  iado=i4o+i3p+noc(isa)*i4p                             4d20s17
                  iadh=ihessa(isa,isb)+i1m+idoub(isbo)*(i2m+idoub(isao)*4d20s17
     $                 (i3+iacto(isa)*i4))                              4d20s17
c     my best guess
                  bc(iadh)=bc(iadh)-4d0*bc(iado)                        10d30s17
                 end do                                                  4d20s17
                end do                                                   4d20s17
              end if                                                    4d20s17
              i4o=i4o+nrowo
             end do
             i10=1
            end do
           end if                                                       4d20s17
           end if                                                       9d25s17
          end if                                                        7d20s17
         end do                                                         7d20s17
         if(iacto(isav)*iacto(isbv)*idoub(isao)*idoub(isbo).gt.0        11d30s17
     $     .and.lsym)then                                               11d30s17
c                                                                       9d15s17
c     isd: symmetry of 2-particle density                               9d15s17
c     isi: symmetry of integrals                                        9d15s17
c                                                                       9d15s17
c
c     aaaa j type part
c
          if(idoit(5).ne.0)then                                         9d25s17
           do isd=1,nsdlk                                                 9d15s17
            if(jdenpt(isd).gt.0)then
             if(isblk(1,isd).eq.isav.and.isblk(2,isd).eq.isbv)then      11d6s17
              if(isav.eq.isbv)then                                        9d15s17
               nrowd=(iacto(isav)*(iacto(isav)+1))/2                      9d15s17
              else                                                        9d15s17
               nrowd=iacto(isav)*iacto(isbv)                              9d15s17
              end if                                                      9d15s17
              igot=0                                                      9d15s17
              do isi=1,nsdlk                                              9d15s17
               if(isblk(1,isi).eq.isao.and.isblk(2,isi).eq.isbo.and.      9d15s17
     $          isblk(3,isi).eq.isblk(3,isd).and.                       9d15s17
     $            isblk(4,isi).eq.isblk(4,isd))then                     9d15s17
                lprint=isblk(1,isi).eq.2.and.isblk(2,isi).eq.2.and.
     $            isblk(3,isi)
     $            .eq.3.and.isblk(1,isd).eq.1.and.isblk(2,isd).eq.1
     $            .and.isblk(3,isd).eq.3
                igot=1                                                    9d15s17
                is3=isblk(3,isi)                                          9d15s17
                is4=isblk(4,isi)                                          9d15s17
                call ilimts(noc(is3),noc(is4),mynprocg,mynowprog,il,ih,   9d15s17
     $             i1s,i1e,i2s,i2e)                                     9d15s17
                i10=i1s                                                   9d15s17
                i1n=noc(is3)                                              9d15s17
                if(isao.eq.isbo)then                                      9d15s17
                 nrowi=(noc(isao)*(noc(isao)+1))/2                        9d15s17
                else                                                      9d15s17
                 nrowi=noc(isao)*noc(isbo)                                9d15s17
                end if                                                    9d15s17
                i4o=ioooo(isi)                                            9d15s17
                do i2=i2s,i2e                                             9d15s17
                 i2m=i2-1-idoub(is4)                                      9d15s17
                 if(i2.eq.i2e)i1n=i1e                                     9d15s17
                 do i1=i10,i1n                                            9d15s17
                  if(i1.gt.idoub(is3).and.i2.gt.idoub(is4))then           9d15s17
                   i1m=i1-1-idoub(is3)                                    9d15s17
                   in=min(i1m,i2m)                                        9d15s17
                   ix=max(i1m,i2m)                                        9d15s17
                   ieq=((ix*(ix+1))/2)+in                                 9d15s17
                   inot=i1m+iacto(is3)*i2m                                9d15s17
                   iswitch=iabs(is3-is4)                                  9d15s17
                   iswitch=iswitch/max(1,iswitch)                         9d15s17
                   iden=jdenpt(isd)+nrowd*((inot-ieq)*iswitch+ieq)        9d15s17
                   do i4=0,idoub(isbo)-1                                  9d15s17
                    do i3=0,idoub(isao)-1                                 9d15s17
                     in=min(i3,i4)                                        9d15s17
                     ix=max(i3,i4)                                        9d15s17
                     ieq=((ix*(ix+1))/2)+in                               9d15s17
                     inot=i3+noc(isao)*i4                                 9d15s17
                     iswitch=iabs(isao-isbo)                              9d15s17
                     iswitch=iswitch/max(1,iswitch)                       9d15s17
                     iad=(inot-ieq)*iswitch+ieq                           9d15s17
                     fact=bc(iad+i4o)*4d0                                 9d15s17
                     jhessa=ihessa(isa,isb)+i4+idoub(isbo)*i3             9d15s17
                     do i6=0,iacto(isbv)-1                                9d15s17
                      do i5=0,iacto(isav)-1                               9d15s17
                       in=min(i5,i6)                                      9d15s17
                       ix=max(i5,i6)                                      9d15s17
                       ieq=((ix*(ix+1))/2)+in                             9d15s17
                       inot=i5+iacto(isav)*i6                             9d15s17
                       iswitch=iabs(isav-isbv)                            9d15s17
                       iswitch=iswitch/max(1,iswitch)                     9d15s17
                       iadd=(inot-ieq)*iswitch+ieq+iden                   9d15s17
                       khessa=jhessa+nrowa*(i5+iacto(isav)*i6)            9d15s17
                       bc(khessa)=bc(khessa)+fact*bc(iadd)                9d15s17
                      end do                                              9d15s17
                     end do                                               9d15s17
                    end do                                                9d15s17
                   end do                                                 9d15s17
                  end if                                                  9d15s17
                  i4o=i4o+nrowi                                           9d15s17
                 end do                                                   9d15s17
                 i10=1                                                    9d15s17
                end do                                                    9d15s17
               end if                                                     9d15s17
              end do                                                      9d15s17
              if(igot.eq.0)then
               write(6,*)('fail no. 1')
               write(6,*)('failed to match integral symmetry type ')
               do isi=1,nsdlk
                write(6,12)isi,(isblk(j,isi),j=1,4)
               end do
               call dws_sync
               call dws_finalize
               stop
              end if
             else if(isblk(2,isd).eq.isav.and.isblk(1,isd).eq.isbv)then 11d6s17
              nrowd=iacto(isav)*iacto(isbv)                              9d15s17
              igot=0                                                      9d15s17
              do isi=1,nsdlk                                              9d15s17
               if(isblk(2,isi).eq.isao.and.isblk(1,isi).eq.isbo.and.    11d6s17
     $          isblk(3,isi).eq.isblk(3,isd).and.                       9d15s17
     $              isblk(4,isi).eq.isblk(4,isd))then                     9d15s17
                lprint=isblk(1,isd).eq.2.and.isblk(2,isd).eq.1.and.
     $              isblk(3,isd).eq.2.and.isblk(4,isd).eq.1
                igot=1                                                    9d15s17
                is3=isblk(3,isi)                                          9d15s17
                is4=isblk(4,isi)                                          9d15s17
                call ilimts(noc(is3),noc(is4),mynprocg,mynowprog,il,ih,   9d15s17
     $             i1s,i1e,i2s,i2e)                                     9d15s17
                i10=i1s                                                   9d15s17
                i1n=noc(is3)                                              9d15s17
                nrowi=noc(isao)*noc(isbo)                                9d15s17
                i4o=ioooo(isi)                                            9d15s17
                do i2=i2s,i2e                                             9d15s17
                 i2m=i2-1-idoub(is4)                                      9d15s17
                 if(i2.eq.i2e)i1n=i1e                                     9d15s17
                 do i1=i10,i1n                                            9d15s17
                  if(i1.gt.idoub(is3).and.i2.gt.idoub(is4))then           9d15s17
                   i1m=i1-1-idoub(is3)                                    9d15s17
                   inot=i1m+iacto(is3)*i2m                                9d15s17
                   iden=jdenpt(isd)+nrowd*inot                          11d6s17
                   do i4=0,idoub(isbo)-1                                  9d15s17
                    do i3=0,idoub(isao)-1                                 9d15s17
                     inot=i4+noc(isbo)*i3                               11d6s17
                     fact=bc(inot+i4o)*2d0                              11d6s17
                     jhessa=ihessa(isa,isb)+i4+idoub(isbo)*i3             9d15s17
                     do i6=0,iacto(isbv)-1                                9d15s17
                      do i5=0,iacto(isav)-1                               9d15s17
                       inot=i6+iacto(isbv)*i5                           11d6s17
                       iadd=inot+iden                                   11d6s17
                       khessa=jhessa+nrowa*(i5+iacto(isav)*i6)            9d15s17
                       bc(khessa)=bc(khessa)+fact*bc(iadd)                9d15s17
                      end do                                              9d15s17
                     end do                                               9d15s17
                    end do                                                9d15s17
                   end do                                                 9d15s17
                  end if                                                  9d15s17
                  i4o=i4o+nrowi                                           9d15s17
                 end do                                                   9d15s17
                 i10=1                                                    9d15s17
                end do                                                    9d15s17
               else if(isblk(1,isi).eq.isao.and.isblk(2,isi).eq.isbo.   11d6s17
     $          and.isblk(3,isi).eq.isblk(3,isd).and.                   11d6s17
     $              isblk(4,isi).eq.isblk(4,isd))then                     9d15s17
                igot=1                                                    9d15s17
                is3=isblk(3,isi)                                          9d15s17
                is4=isblk(4,isi)                                          9d15s17
                call ilimts(noc(is3),noc(is4),mynprocg,mynowprog,il,ih,   9d15s17
     $             i1s,i1e,i2s,i2e)                                     9d15s17
                i10=i1s                                                   9d15s17
                i1n=noc(is3)                                              9d15s17
                nrowi=noc(isao)*noc(isbo)                                9d15s17
                i4o=ioooo(isi)                                            9d15s17
                do i2=i2s,i2e                                             9d15s17
                 i2m=i2-1-idoub(is4)                                      9d15s17
                 if(i2.eq.i2e)i1n=i1e                                     9d15s17
                 do i1=i10,i1n                                            9d15s17
                  if(i1.gt.idoub(is3).and.i2.gt.idoub(is4))then           9d15s17
                   i1m=i1-1-idoub(is3)                                    9d15s17
                   inot=i1m+iacto(is3)*i2m                                9d15s17
                   iden=jdenpt(isd)+nrowd*inot                          11d6s17
                   do i4=0,idoub(isbo)-1                                  9d15s17
                    do i3=0,idoub(isao)-1                                 9d15s17
                     inot=i3+noc(isao)*i4                               11d6s17
                     fact=bc(inot+i4o)*2d0                              11d6s17
                     jhessa=ihessa(isa,isb)+i4+idoub(isbo)*i3             9d15s17
                     do i6=0,iacto(isbv)-1                                9d15s17
                      do i5=0,iacto(isav)-1                               9d15s17
                       inot=i6+iacto(isbv)*i5                           11d6s17
                       iadd=inot+iden                                   11d6s17
                       khessa=jhessa+nrowa*(i5+iacto(isav)*i6)            9d15s17
                       bc(khessa)=bc(khessa)+fact*bc(iadd)                9d15s17
                      end do                                              9d15s17
                     end do                                               9d15s17
                    end do                                                9d15s17
                   end do                                                 9d15s17
                  end if                                                  9d15s17
                  i4o=i4o+nrowi                                           9d15s17
                 end do                                                   9d15s17
                 i10=1                                                    9d15s17
                end do                                                    9d15s17
               end if                                                     9d15s17
              end do                                                      9d15s17
              if(igot.eq.0)then
               write(6,*)('aaaa hessa a2: '),
     $              ('failed to match integral symmetry type '),isao,
     $              isbo
               write(6,12)isd,(isblk(j,isd),j=1,4)
               do isi=1,nsdlk
                write(6,12)isi,(isblk(j,isi),j=1,4)
               end do
               call dws_sync
               call dws_finalize
               stop
              end if
             end if                                                       9d15s17
            end if                                                      11d6s17
           end do                                                        9d15s17
c
c     aaaa k type part
c
          do isd=1,nsdlk                                                 9d15s17
           if(jdenpt(isd).gt.0)then                                     10d27s17
           if(isblk(1,isd).eq.isblk(2,isd))then                         9d15s17
            if(isblk(1,isd).eq.isblk(3,isd))then                        9d15s17
             fmult=8d0                                                  9d15s17
            else                                                        9d15s17
c
c     this is my guess
c
             fmult=8d0                                                  11d6s17
            end if                                                      9d15s17
           else                                                         9d15s17
            fmult=2d0                                                   9d15s17
           end if                                                       9d15s17
           if(isblk(1,isd).eq.isav.and.isblk(3,isd).eq.isbv.and.        9d15s17
     $          jdenpt(isd).gt.0)then                                   9d15s17
            isd2=isblk(2,isd)
            isd4=isblk(4,isd)                                           9d15s17
            if(isav.eq.isd2)then                                        9d15s17
             nrowd=(iacto(isav)*(iacto(isav)+1))/2                      9d15s17
            else                                                        9d15s17
             nrowd=iacto(isav)*iacto(isd2)                              9d15s17
            end if                                                      9d15s17
            igot=0                                                      9d15s17
            do isi=1,nsdlk                                              9d15s17
             if(isblk(1,isi).eq.isao.and.isblk(3,isi).eq.isbo.and.      9d15s17
     $          isblk(2,isi).eq.isblk(2,isd).and.                       9d15s17
     $            isblk(4,isi).eq.isblk(4,isd))then                     9d15s17
              igot=1                                                    9d15s17
              call ilimts(noc(isbo),noc(isd4),mynprocg,mynowprog,il,ih,  9d15s17
     $             i1s,i1e,i2s,i2e)                                     9d15s17
              i10=i1s                                                   9d15s17
              i1n=noc(isbo)                                             9d15s17
              if(isao.eq.isd2)then                                       9d15s17
               nrowi=(noc(isao)*(noc(isao)+1))/2                        9d15s17
              else                                                      9d15s17
               nrowi=noc(isao)*noc(isd2)                                 9d15s17
              end if                                                    9d15s17
              i4o=ioooo(isi)                                            9d15s17
              do i2=i2s,i2e                                             9d15s17
               i2m=i2-1-idoub(isd4)                                      9d15s17
               if(i2.eq.i2e)i1n=i1e                                     9d15s17
               do i1=i10,i1n                                            9d15s17
                if(i1.le.idoub(isbo).and.i2.gt.idoub(isd4))then          9d15s17
                 i1m=i1-1                                               9d15s17
                 do i6=0,iacto(isbv)-1                                  9d15s17
                  in=min(i2m,i6)                                        9d15s17
                  ix=max(i2m,i6)                                        9d15s17
                  ieq=((ix*(ix+1))/2)+in                                9d15s17
                  inot=i6+iacto(isbv)*i2m                               9d15s17
                  iswitch=iabs(isd4-isbv)                               9d15s17
                  iswitch=iswitch/max(1,iswitch)                        9d15s17
                  iden=nrowd*((inot-ieq)*iswitch+ieq)+jdenpt(isd)       9d15s17
                  do i4=0,iacto(isd2)-1                                  9d15s17
                   i4p=i4+idoub(isd2)                                    9d15s17
                   do i5=0,iacto(isav)-1                                9d15s17
                    jhessa=ihessa(isa,isb)+nrowa*(i5+iacto(isav)*i6)+i1m9d15s17
                    in=min(i5,i4)                                       9d15s17
                    ix=max(i5,i4)                                       9d15s17
                    ieq=((ix*(ix+1))/2)+in                              9d15s17
                    inot=i5+iacto(isav)*i4                              9d15s17
                    iswitch=iabs(isd2-isav)                             9d15s17
                    iswitch=iswitch/max(1,iswitch)                      9d15s17
                    iadd=iden+(inot-ieq)*iswitch+ieq                    9d15s17
                    do i3=0,idoub(isao)-1                               9d15s17
                     in=min(i3,i4p)                                     9d15s17
                     ix=max(i3,i4p)                                     9d15s17
                     ieq=((ix*(ix+1))/2)+in                             9d15s17
                     inot=i3+noc(isao)*i4p                              9d15s17
                     iswitch=iabs(isao-isd2)                            9d15s17
                     iswitch=iswitch/max(1,iswitch)                     9d15s17
                     iad=(inot-ieq)*iswitch+ieq                         9d15s17
                     fact=bc(iad+i4o)*fmult                             9d15s17
                     khessa=jhessa+idoub(isbo)*i3                       9d15s17
                     bc(khessa)=bc(khessa)+fact*bc(iadd)                9d15s17
                    end do                                              9d15s17
                   end do                                               9d15s17
                  end do                                                9d15s17
                 end do                                                 9d15s17
                end if                                                  9d15s17
                i4o=i4o+nrowi                                           9d15s17
               end do                                                   9d15s17
               i10=1                                                    9d15s17
              end do                                                    9d15s17
             else if(isblk(2,isi).eq.isao.and.isblk(3,isi).eq.isbo.and. 11d7s17
     $          isblk(1,isi).eq.isblk(2,isd).and.                       11d7s17
     $            isblk(4,isi).eq.isblk(4,isd))then                     9d15s17
              igot=1                                                    9d15s17
              call ilimts(noc(isbo),noc(isd4),mynprocg,mynowprog,il,ih,  9d15s17
     $             i1s,i1e,i2s,i2e)                                     9d15s17
              i10=i1s                                                   9d15s17
              i1n=noc(isbo)                                             9d15s17
              if(isao.eq.isd2)then                                       9d15s17
               nrowi=(noc(isao)*(noc(isao)+1))/2                        9d15s17
              else                                                      9d15s17
               nrowi=noc(isao)*noc(isd2)                                 9d15s17
              end if                                                    9d15s17
              i4o=ioooo(isi)                                            9d15s17
              do i2=i2s,i2e                                             9d15s17
               i2m=i2-1-idoub(isd4)                                      9d15s17
               if(i2.eq.i2e)i1n=i1e                                     9d15s17
               do i1=i10,i1n                                            9d15s17
                if(i1.le.idoub(isbo).and.i2.gt.idoub(isd4))then          9d15s17
                 i1m=i1-1                                               9d15s17
                 do i6=0,iacto(isbv)-1                                  9d15s17
                  in=min(i2m,i6)                                        9d15s17
                  ix=max(i2m,i6)                                        9d15s17
                  ieq=((ix*(ix+1))/2)+in                                9d15s17
                  inot=i6+iacto(isbv)*i2m                               9d15s17
                  iswitch=iabs(isd4-isbv)                               9d15s17
                  iswitch=iswitch/max(1,iswitch)                        9d15s17
                  iden=nrowd*((inot-ieq)*iswitch+ieq)+jdenpt(isd)       9d15s17
                  do i4=0,iacto(isd2)-1                                  9d15s17
                   i4p=i4+idoub(isd2)                                    9d15s17
                   do i5=0,iacto(isav)-1                                9d15s17
                    jhessa=ihessa(isa,isb)+nrowa*(i5+iacto(isav)*i6)+i1m9d15s17
                    in=min(i5,i4)                                       9d15s17
                    ix=max(i5,i4)                                       9d15s17
                    ieq=((ix*(ix+1))/2)+in                              9d15s17
                    inot=i5+iacto(isav)*i4                              9d15s17
                    iswitch=iabs(isd2-isav)                             9d15s17
                    iswitch=iswitch/max(1,iswitch)                      9d15s17
                    iadd=iden+(inot-ieq)*iswitch+ieq                    9d15s17
                    do i3=0,idoub(isao)-1                               9d15s17
                     inot=i4p+noc(isd2)*i3                              11d7s17
                     fact=bc(inot+i4o)*fmult                            11d7s17
                     khessa=jhessa+idoub(isbo)*i3                       9d15s17
                     bc(khessa)=bc(khessa)+fact*bc(iadd)                9d15s17
                    end do                                              9d15s17
                   end do                                               9d15s17
                  end do                                                9d15s17
                 end do                                                 9d15s17
                end if                                                  9d15s17
                i4o=i4o+nrowi                                           9d15s17
               end do                                                   9d15s17
               i10=1                                                    9d15s17
              end do                                                    9d15s17
             end if                                                     9d15s17
             if(isblk(2,isi).eq.isao.and.isblk(4,isi).eq.isbo.and.      9d15s17
     $          isblk(1,isi).eq.isblk(2,isd).and.                       9d15s17
     $            isblk(3,isi).eq.isblk(4,isd).and.igot.eq.0)then       9d15s17
              igot=1                                                    9d15s17
              is1=isblk(1,isi)                                          9d15s17
              call ilimts(noc(isd4),noc(isbo),mynprocg,mynowprog,il,ih,  9d15s17
     $             i1s,i1e,i2s,i2e)                                     9d15s17
              i10=i1s                                                   9d15s17
              i1n=noc(isd4)                                              9d15s17
              if(isao.eq.isd2)then                                       9d15s17
               nrowi=(noc(isao)*(noc(isao)+1))/2                        9d15s17
              else                                                      9d15s17
               nrowi=noc(isao)*noc(isd2)                                 9d15s17
              end if                                                    9d15s17
              i4o=ioooo(isi)                                            9d15s17
              do i2=i2s,i2e                                             9d15s17
               i2m=i2-1                                                 9d15s17
               if(i2.eq.i2e)i1n=i1e                                     9d15s17
               do i1=i10,i1n                                            9d15s17
                if(i2.le.idoub(isbo).and.i1.gt.idoub(isd4))then          9d15s17
                 i1m=i1-1-idoub(isd4)                                    9d15s17
                 do i6=0,iacto(isbv)-1                                  9d15s17
                  in=min(i1m,i6)                                        9d15s17
                  ix=max(i1m,i6)                                        9d15s17
                  ieq=((ix*(ix+1))/2)+in                                9d15s17
                  inot=i6+iacto(isbv)*i1m                               9d15s17
                  iswitch=iabs(isbv-isd4)                               9d15s17
                  iswitch=iswitch/max(1,iswitch)                        9d15s17
                  iden=jdenpt(isd)+nrowd*((inot-ieq)*iswitch+ieq)       9d15s17
                  do i4=0,iacto(isd2)-1                                  9d15s17
                   i4p=i4+idoub(isd2)                                    9d15s17
                   do i5=0,iacto(isav)-1                                9d15s17
                    jhessa=ihessa(isa,isb)+i2m+nrowa*(i5+iacto(isav)*i6)9d15s17
                    in=min(i5,i4)                                       9d15s17
                    ix=max(i5,i4)                                       9d15s17
                    ieq=((ix*(ix+1))/2)+in                              9d15s17
                    inot=i5+iacto(isav)*i4                              9d15s17
                    iswitch=iabs(isav-isd2)                             9d15s17
                    iswitch=iswitch/max(1,iswitch)                      9d15s17
                    iadd=iden+(inot-ieq)*iswitch+ieq                    9d15s17
                    do i3=0,idoub(isao)-1                               9d15s17
                     in=min(i3,i4p)                                     9d15s17
                     ix=max(i3,i4p)                                     9d15s17
                     ieq=((ix*(ix+1))/2)+in                             9d15s17
                     inot=i4p+noc(isd2)*i3                               9d15s17
                     iswitch=iabs(isao-isd2)                            9d15s17
                     iswitch=iswitch/max(1,iswitch)                     9d15s17
                     iad=(inot-ieq)*iswitch+ieq                         9d15s17
                     fact=bc(iad+i4o)*fmult                             9d15s17
                     khessa=jhessa+idoub(isbo)*i3                       9d15s17
                     bc(khessa)=bc(khessa)+fact*bc(iadd)                9d15s17
                    end do                                              9d15s17
                   end do                                               9d15s17
                  end do                                                9d15s17
                 end do                                                 9d15s17
                end if                                                  9d15s17
                i4o=i4o+nrowi                                           9d15s17
               end do                                                   9d15s17
               i10=1                                                    9d15s17
              end do                                                    9d15s17
             else if(isblk(1,isi).eq.isao.and.isblk(4,isi).eq.isbo.and. 11d7s17
     $          isblk(2,isi).eq.isblk(2,isd).and.                       11d7s17
     $            isblk(3,isi).eq.isblk(4,isd).and.igot.eq.0)then       9d15s17
              igot=1                                                    9d15s17
              is1=isblk(1,isi)                                          9d15s17
              call ilimts(noc(isd4),noc(isbo),mynprocg,mynowprog,il,ih,  9d15s17
     $             i1s,i1e,i2s,i2e)                                     9d15s17
              i10=i1s                                                   9d15s17
              i1n=noc(isd4)                                              9d15s17
              if(isao.eq.isd2)then                                       9d15s17
               nrowi=(noc(isao)*(noc(isao)+1))/2                        9d15s17
              else                                                      9d15s17
               nrowi=noc(isao)*noc(isd2)                                 9d15s17
              end if                                                    9d15s17
              i4o=ioooo(isi)                                            9d15s17
              do i2=i2s,i2e                                             9d15s17
               i2m=i2-1                                                 9d15s17
               if(i2.eq.i2e)i1n=i1e                                     9d15s17
               do i1=i10,i1n                                            9d15s17
                if(i2.le.idoub(isbo).and.i1.gt.idoub(isd4))then          9d15s17
                 i1m=i1-1-idoub(isd4)                                    9d15s17
                 do i6=0,iacto(isbv)-1                                  9d15s17
                  in=min(i1m,i6)                                        9d15s17
                  ix=max(i1m,i6)                                        9d15s17
                  ieq=((ix*(ix+1))/2)+in                                9d15s17
                  inot=i6+iacto(isbv)*i1m                               9d15s17
                  iswitch=iabs(isbv-isd4)                               9d15s17
                  iswitch=iswitch/max(1,iswitch)                        9d15s17
                  iden=jdenpt(isd)+nrowd*((inot-ieq)*iswitch+ieq)       9d15s17
                  do i4=0,iacto(isd2)-1                                  9d15s17
                   i4p=i4+idoub(isd2)                                    9d15s17
                   do i5=0,iacto(isav)-1                                9d15s17
                    jhessa=ihessa(isa,isb)+i2m+nrowa*(i5+iacto(isav)*i6)9d15s17
                    in=min(i5,i4)                                       9d15s17
                    ix=max(i5,i4)                                       9d15s17
                    ieq=((ix*(ix+1))/2)+in                              9d15s17
                    inot=i5+iacto(isav)*i4                              9d15s17
                    iswitch=iabs(isav-isd2)                             9d15s17
                    iswitch=iswitch/max(1,iswitch)                      9d15s17
                    iadd=iden+(inot-ieq)*iswitch+ieq                    9d15s17
                    do i3=0,idoub(isao)-1                               9d15s17
                     inot=i3+noc(isao)*i4p                              11d7s17
                     fact=bc(inot+i4o)*fmult                            11d7s17
                     khessa=jhessa+idoub(isbo)*i3                       9d15s17
                     bc(khessa)=bc(khessa)+fact*bc(iadd)                9d15s17
                    end do                                              9d15s17
                   end do                                               9d15s17
                  end do                                                9d15s17
                 end do                                                 9d15s17
                end if                                                  9d15s17
                i4o=i4o+nrowi                                           9d15s17
               end do                                                   9d15s17
               i10=1                                                    9d15s17
              end do                                                    9d15s17
             end if                                                     9d15s17
            end do                                                      9d15s17
            if(igot.eq.0)then
               write(6,*)('fail no. 3')
             write(6,*)('for k:'),isao,isbo,isbv,isav
             write(6,*)('2 part density: '),(isblk(j,isd),j=1,4)
             write(6,*)('failed to match integral symmetry type ')
             do isi=1,nsdlk
              write(6,12)isi,(isblk(j,isi),j=1,4)
             end do
             call dws_sync
             call dws_finalize
             stop
            end if
           else if(isblk(2,isd).eq.isav.and.isblk(4,isd).eq.isbv.and.   9d15s17
     $          jdenpt(isd).gt.0)then                                   9d15s17
            isd1=isblk(1,isd)
            isd3=isblk(3,isd)                                           9d15s17
            if(isav.eq.isd1)then                                        9d15s17
             nrowd=(iacto(isav)*(iacto(isav)+1))/2                      9d15s17
            else                                                        9d15s17
             nrowd=iacto(isav)*iacto(isd1)                              9d15s17
            end if                                                      9d15s17
            igot=0                                                      9d15s17
            do isi=1,nsdlk                                              9d15s17
             if(isblk(1,isi).eq.isao.and.isblk(3,isi).eq.isbo.and.      9d15s17
     $          isblk(2,isi).eq.isblk(1,isd).and.                       9d15s17
     $            isblk(4,isi).eq.isblk(3,isd))then                     9d15s17
              igot=1                                                    9d15s17
              call ilimts(noc(isbo),noc(isd3),mynprocg,mynowprog,il,ih,  9d15s17
     $             i1s,i1e,i2s,i2e)                                     9d15s17
              i10=i1s                                                   9d15s17
              i1n=noc(isbo)                                             9d15s17
              if(isao.eq.isd1)then                                       9d15s17
               nrowi=(noc(isao)*(noc(isao)+1))/2                        9d15s17
              else                                                      9d15s17
               nrowi=noc(isao)*noc(isd1)                                 9d15s17
              end if                                                    9d15s17
              i4o=ioooo(isi)                                            9d15s17
              do i2=i2s,i2e                                             9d15s17
               i2m=i2-1-idoub(isd3)                                      9d15s17
               if(i2.eq.i2e)i1n=i1e                                     9d15s17
               do i1=i10,i1n                                            9d15s17
                if(i1.le.idoub(isbo).and.i2.gt.idoub(isd3))then          9d15s17
                 i1m=i1-1                                               9d15s17
                 do i6=0,iacto(isbv)-1                                  9d15s17
                  in=min(i2m,i6)                                        9d15s17
                  ix=max(i2m,i6)                                        9d15s17
                  ieq=((ix*(ix+1))/2)+in                                9d15s17
                  inot=i2m+iacto(isd3)*i6                               11d6s17
                  iswitch=iabs(isd3-isbv)                               9d15s17
                  iswitch=iswitch/max(1,iswitch)                        9d15s17
                  iden=nrowd*((inot-ieq)*iswitch+ieq)+jdenpt(isd)       9d15s17
                  do i4=0,iacto(isd1)-1                                  9d15s17
                   i4p=i4+idoub(isd1)                                    9d15s17
                   do i5=0,iacto(isav)-1                                9d15s17
                    jhessa=ihessa(isa,isb)+nrowa*(i5+iacto(isav)*i6)+i1m9d15s17
                    in=min(i5,i4)                                       9d15s17
                    ix=max(i5,i4)                                       9d15s17
                    ieq=((ix*(ix+1))/2)+in                              9d15s17
                    inot=i4+iacto(isd1)*i5                              9d15s17
                    iswitch=iabs(isd1-isav)                             9d15s17
                    iswitch=iswitch/max(1,iswitch)                      9d15s17
                    iadd=iden+(inot-ieq)*iswitch+ieq                    9d15s17
                    do i3=0,idoub(isao)-1                               9d15s17
                     in=min(i3,i4p)                                     9d15s17
                     ix=max(i3,i4p)                                     9d15s17
                     ieq=((ix*(ix+1))/2)+in                             9d15s17
                     inot=i3+noc(isao)*i4p                              9d15s17
                     iswitch=iabs(isao-isd1)                            9d15s17
                     iswitch=iswitch/max(1,iswitch)                     9d15s17
                     iad=(inot-ieq)*iswitch+ieq                         9d15s17
                     fact=bc(iad+i4o)*fmult                             9d15s17
                     khessa=jhessa+idoub(isbo)*i3                       9d15s17
                     bc(khessa)=bc(khessa)+fact*bc(iadd)                9d15s17
                    end do                                              9d15s17
                   end do                                               9d15s17
                  end do                                                9d15s17
                 end do                                                 9d15s17
                end if                                                  9d15s17
                i4o=i4o+nrowi                                           9d15s17
               end do                                                   9d15s17
               i10=1                                                    9d15s17
              end do                                                    9d15s17
             end if                                                     9d15s17
             if(isblk(2,isi).eq.isao.and.isblk(4,isi).eq.isbo.and.      9d15s17
     $          isblk(1,isi).eq.isblk(1,isd).and.                       9d15s17
     $            isblk(3,isi).eq.isblk(3,isd).and.igot.eq.0)then       9d15s17
              igot=1                                                    9d15s17
              call ilimts(noc(isd3),noc(isbo),mynprocg,mynowprog,il,ih,  9d15s17
     $             i1s,i1e,i2s,i2e)                                     9d15s17
              i10=i1s                                                   9d15s17
              i1n=noc(isd3)                                              9d15s17
              if(isao.eq.isd1)then                                       9d15s17
               nrowi=(noc(isao)*(noc(isao)+1))/2                        9d15s17
              else                                                      9d15s17
               nrowi=noc(isao)*noc(isd1)                                 9d15s17
              end if                                                    9d15s17
              i4o=ioooo(isi)                                            9d15s17
              do i2=i2s,i2e                                             9d15s17
               i2m=i2-1                                                 9d15s17
               if(i2.eq.i2e)i1n=i1e                                     9d15s17
               do i1=i10,i1n                                            9d15s17
                if(i2.le.idoub(isbo).and.i1.gt.idoub(isd3))then          9d15s17
                 i1m=i1-1-idoub(isd3)                                    9d15s17
                 do i6=0,iacto(isbv)-1                                  9d15s17
                  in=min(i1m,i6)                                        9d15s17
                  ix=max(i1m,i6)                                        9d15s17
                  ieq=((ix*(ix+1))/2)+in                                9d15s17
                  inot=i1m+iacto(isd3)*i6                               9d15s17
                  iswitch=iabs(isbv-isd3)                               9d15s17
                  iswitch=iswitch/max(1,iswitch)                        9d15s17
                  iden=jdenpt(isd)+nrowd*((inot-ieq)*iswitch+ieq)       9d15s17
                  do i4=0,iacto(isd1)-1                                  9d15s17
                   i4p=i4+idoub(isd1)                                    9d15s17
                   do i5=0,iacto(isav)-1                                9d15s17
                    jhessa=ihessa(isa,isb)+i2m+nrowa*(i5+iacto(isav)*i6)9d15s17
                    in=min(i5,i4)                                       9d15s17
                    ix=max(i5,i4)                                       9d15s17
                    ieq=((ix*(ix+1))/2)+in                              9d15s17
                    inot=i4+iacto(isd1)*i5                              9d15s17
                    iswitch=iabs(isav-isd1)                             9d15s17
                    iswitch=iswitch/max(1,iswitch)                      9d15s17
                    iadd=iden+(inot-ieq)*iswitch+ieq                    9d15s17
                    do i3=0,idoub(isao)-1                               9d15s17
                     in=min(i3,i4p)                                     9d15s17
                     ix=max(i3,i4p)                                     9d15s17
                     ieq=((ix*(ix+1))/2)+in                             9d15s17
                     inot=i4p+noc(isd1)*i3                               9d15s17
                     iswitch=iabs(isao-isd1)                            9d15s17
                     iswitch=iswitch/max(1,iswitch)                     9d15s17
                     iad=(inot-ieq)*iswitch+ieq                         9d15s17
                     fact=bc(iad+i4o)*fmult                             9d15s17
                     khessa=jhessa+idoub(isbo)*i3                       9d15s17
                     bc(khessa)=bc(khessa)+fact*bc(iadd)                9d15s17
                    end do                                              9d15s17
                   end do                                               9d15s17
                  end do                                                9d15s17
                 end do                                                 9d15s17
                end if                                                  9d15s17
                i4o=i4o+nrowi                                           9d15s17
               end do                                                   9d15s17
               i10=1                                                    9d15s17
              end do                                                    9d15s17
             end if                                                     9d15s17
            end do                                                      9d15s17
            if(igot.eq.0)then
               write(6,*)('fail no. 4')
             write(6,*)('for k:')
             write(6,*)('failed to match integral symmetry type ')
             do isi=1,nsdlk
              write(6,12)isi,(isblk(j,isi),j=1,4)
             end do
             call dws_sync
             call dws_finalize
             stop
            end if
           else if(isblk(1,isd).eq.isav.and.isblk(4,isd).eq.isbv.and.   11d6s17
     $          jdenpt(isd).gt.0)then                                   11d6s17
            isd2=isblk(2,isd)                                           11d6s17
            isd3=isblk(3,isd)                                           9d15s17
            if(isav.eq.isd2)then                                        11d6s17
             nrowd=(iacto(isav)*(iacto(isav)+1))/2                      9d15s17
            else                                                        9d15s17
             nrowd=iacto(isav)*iacto(isd2)                              11d6s17
            end if                                                      9d15s17
            igot=0                                                      9d15s17
            do isi=1,nsdlk                                              9d15s17
             if(isblk(1,isi).eq.isao.and.isblk(3,isi).eq.isbo.and.      9d15s17
     $          isblk(2,isi).eq.isd2.and.                               11d6s17
     $            isblk(4,isi).eq.isd3)then                             11d6s17
              igot=1                                                    9d15s17
              call ilimts(noc(isbo),noc(isd3),mynprocg,mynowprog,il,ih,  9d15s17
     $             i1s,i1e,i2s,i2e)                                     9d15s17
              i10=i1s                                                   9d15s17
              i1n=noc(isbo)                                             9d15s17
              if(isao.eq.isd2)then                                      11d6s17
               nrowi=(noc(isao)*(noc(isao)+1))/2                        9d15s17
              else                                                      9d15s17
               nrowi=noc(isao)*noc(isd2)                                11d6s17
              end if                                                    9d15s17
              i4o=ioooo(isi)                                            9d15s17
              do i2=i2s,i2e                                             9d15s17
               i2m=i2-1-idoub(isd3)                                      9d15s17
               if(i2.eq.i2e)i1n=i1e                                     9d15s17
               do i1=i10,i1n                                            9d15s17
                if(i1.le.idoub(isbo).and.i2.gt.idoub(isd3))then          9d15s17
                 i1m=i1-1                                               9d15s17
                 do i6=0,iacto(isbv)-1                                  9d15s17
                  in=min(i2m,i6)                                        9d15s17
                  ix=max(i2m,i6)                                        9d15s17
                  ieq=((ix*(ix+1))/2)+in                                9d15s17
                  inot=i2m+iacto(isd3)*i6                               11d6s17
                  iswitch=iabs(isd3-isbv)                               9d15s17
                  iswitch=iswitch/max(1,iswitch)                        9d15s17
                  iden=nrowd*((inot-ieq)*iswitch+ieq)+jdenpt(isd)       9d15s17
                  do i4=0,iacto(isd2)-1                                 11d6s17
                   i4p=i4+idoub(isd2)                                   11d6s17
                   do i5=0,iacto(isav)-1                                9d15s17
                    jhessa=ihessa(isa,isb)+nrowa*(i5+iacto(isav)*i6)+i1m9d15s17
                    in=min(i5,i4)                                       9d15s17
                    ix=max(i5,i4)                                       9d15s17
                    ieq=((ix*(ix+1))/2)+in                              9d15s17
                    inot=i5+iacto(isav)*i4                              11d6s17
                    iswitch=iabs(isd2-isav)                             11d6s17
                    iswitch=iswitch/max(1,iswitch)                      9d15s17
                    iadd=iden+(inot-ieq)*iswitch+ieq                    9d15s17
                    do i3=0,idoub(isao)-1                               9d15s17
                     in=min(i3,i4p)                                     9d15s17
                     ix=max(i3,i4p)                                     9d15s17
                     ieq=((ix*(ix+1))/2)+in                             9d15s17
                     inot=i3+noc(isao)*i4p                              9d15s17
                     iswitch=iabs(isao-isd2)                            11d6s17
                     iswitch=iswitch/max(1,iswitch)                     9d15s17
                     iad=(inot-ieq)*iswitch+ieq                         9d15s17
                     fact=bc(iad+i4o)*fmult                             9d15s17
                     khessa=jhessa+idoub(isbo)*i3                       9d15s17
                     bc(khessa)=bc(khessa)+fact*bc(iadd)                9d15s17
                    end do                                              9d15s17
                   end do                                               9d15s17
                  end do                                                9d15s17
                 end do                                                 9d15s17
                end if                                                  9d15s17
                i4o=i4o+nrowi                                           9d15s17
               end do                                                   9d15s17
               i10=1                                                    9d15s17
              end do                                                    9d15s17
             end if                                                     9d15s17
             if(isblk(2,isi).eq.isao.and.isblk(3,isi).eq.isbo.and.      12d25s17
     $          isblk(1,isi).eq.isd2.and.                               12d25s17
     $            isblk(4,isi).eq.isd3)then                             11d6s17
              igot=1                                                    9d15s17
              call ilimts(noc(isbo),noc(isd3),mynprocg,mynowprog,il,ih,  9d15s17
     $             i1s,i1e,i2s,i2e)                                     9d15s17
              i10=i1s                                                   9d15s17
              i1n=noc(isbo)                                             9d15s17
              nrowi=noc(isao)*noc(isd2)                                 12d25s17
              i4o=ioooo(isi)                                            9d15s17
              do i2=i2s,i2e                                             9d15s17
               i2m=i2-1-idoub(isd3)                                      9d15s17
               if(i2.eq.i2e)i1n=i1e                                     9d15s17
               do i1=i10,i1n                                            9d15s17
                if(i1.le.idoub(isbo).and.i2.gt.idoub(isd3))then          9d15s17
                 i1m=i1-1                                               9d15s17
                 do i6=0,iacto(isbv)-1                                  9d15s17
                  in=min(i2m,i6)                                        9d15s17
                  ix=max(i2m,i6)                                        9d15s17
                  ieq=((ix*(ix+1))/2)+in                                9d15s17
                  inot=i2m+iacto(isd3)*i6                               11d6s17
                  iswitch=iabs(isd3-isbv)                               9d15s17
                  iswitch=iswitch/max(1,iswitch)                        9d15s17
                  iden=nrowd*((inot-ieq)*iswitch+ieq)+jdenpt(isd)       9d15s17
                  do i4=0,iacto(isd2)-1                                 11d6s17
                   i4p=i4+idoub(isd2)                                   11d6s17
                   do i5=0,iacto(isav)-1                                9d15s17
                    jhessa=ihessa(isa,isb)+nrowa*(i5+iacto(isav)*i6)+i1m9d15s17
                    in=min(i5,i4)                                       9d15s17
                    ix=max(i5,i4)                                       9d15s17
                    ieq=((ix*(ix+1))/2)+in                              9d15s17
                    inot=i5+iacto(isav)*i4                              11d6s17
                    iswitch=iabs(isd2-isav)                             11d6s17
                    iswitch=iswitch/max(1,iswitch)                      9d15s17
                    iadd=iden+(inot-ieq)*iswitch+ieq                    9d15s17
                    do i3=0,idoub(isao)-1                               9d15s17
                     inot=i3+noc(isao)*i4p                              9d15s17
                     fact=bc(i4p+noc(isd2)*i3+i4o)*fmult                12d25s17
                     khessa=jhessa+idoub(isbo)*i3                       9d15s17
                     bc(khessa)=bc(khessa)+fact*bc(iadd)                9d15s17
                    end do                                              9d15s17
                   end do                                               9d15s17
                  end do                                                9d15s17
                 end do                                                 9d15s17
                end if                                                  9d15s17
                i4o=i4o+nrowi                                           9d15s17
               end do                                                   9d15s17
               i10=1                                                    9d15s17
              end do                                                    9d15s17
             end if                                                     9d15s17
             if(isblk(2,isi).eq.isao.and.isblk(4,isi).eq.isbo.and.      9d15s17
     $          isblk(1,isi).eq.isd2.and.                               11d6s17
     $            isblk(3,isi).eq.isd3.and.igot.eq.0)then               11d6s17
              igot=1                                                    9d15s17
              call ilimts(noc(isd3),noc(isbo),mynprocg,mynowprog,il,ih,  9d15s17
     $             i1s,i1e,i2s,i2e)                                     9d15s17
              i10=i1s                                                   9d15s17
              i1n=noc(isd3)                                              9d15s17
              if(isao.eq.isd2)then                                      11d6s17
               nrowi=(noc(isao)*(noc(isao)+1))/2                        9d15s17
              else                                                      9d15s17
               nrowi=noc(isao)*noc(isd2)                                11d6s17
              end if                                                    9d15s17
              i4o=ioooo(isi)                                            9d15s17
              do i2=i2s,i2e                                             9d15s17
               i2m=i2-1                                                 9d15s17
               if(i2.eq.i2e)i1n=i1e                                     9d15s17
               do i1=i10,i1n                                            9d15s17
                if(i2.le.idoub(isbo).and.i1.gt.idoub(isd3))then          9d15s17
                 i1m=i1-1-idoub(isd3)                                    9d15s17
                 do i6=0,iacto(isbv)-1                                  9d15s17
                  in=min(i1m,i6)                                        9d15s17
                  ix=max(i1m,i6)                                        9d15s17
                  ieq=((ix*(ix+1))/2)+in                                9d15s17
                  inot=i1m+iacto(isd3)*i6                               9d15s17
                  iswitch=iabs(isbv-isd3)                               9d15s17
                  iswitch=iswitch/max(1,iswitch)                        9d15s17
                  iden=jdenpt(isd)+nrowd*((inot-ieq)*iswitch+ieq)       9d15s17
                  do i4=0,iacto(isd2)-1                                 11d6s17
                   i4p=i4+idoub(isd2)                                   11d6s17
                   do i5=0,iacto(isav)-1                                9d15s17
                    jhessa=ihessa(isa,isb)+i2m+nrowa*(i5+iacto(isav)*i6)9d15s17
                    in=min(i5,i4)                                       9d15s17
                    ix=max(i5,i4)                                       9d15s17
                    ieq=((ix*(ix+1))/2)+in                              9d15s17
                    inot=i5+iacto(isav)*i4                              11d6s17
                    iswitch=iabs(isav-isd2)                             11d6s17
                    iswitch=iswitch/max(1,iswitch)                      9d15s17
                    iadd=iden+(inot-ieq)*iswitch+ieq                    9d15s17
                    do i3=0,idoub(isao)-1                               9d15s17
                     in=min(i3,i4p)                                     9d15s17
                     ix=max(i3,i4p)                                     9d15s17
                     ieq=((ix*(ix+1))/2)+in                             9d15s17
                     inot=i4p+noc(isd2)*i3                              11d6s17
                     iswitch=iabs(isao-isd2)                            11d6s17
                     iswitch=iswitch/max(1,iswitch)                     9d15s17
                     iad=(inot-ieq)*iswitch+ieq                         9d15s17
                     fact=bc(iad+i4o)*fmult                             9d15s17
                     khessa=jhessa+idoub(isbo)*i3                       9d15s17
                     bc(khessa)=bc(khessa)+fact*bc(iadd)                9d15s17
                    end do                                              9d15s17
                   end do                                               9d15s17
                  end do                                                9d15s17
                 end do                                                 9d15s17
                end if                                                  9d15s17
                i4o=i4o+nrowi                                           9d15s17
               end do                                                   9d15s17
               i10=1                                                    9d15s17
              end do                                                    9d15s17
             end if                                                     9d15s17
             if(isblk(1,isi).eq.isao.and.isblk(4,isi).eq.isbo.and.      12d25s17
     $          isblk(2,isi).eq.isd2.and.                               12d25s17
     $            isblk(3,isi).eq.isd3.and.igot.eq.0)then               11d6s17
              igot=1                                                    9d15s17
              call ilimts(noc(isd3),noc(isbo),mynprocg,mynowprog,il,ih,  9d15s17
     $             i1s,i1e,i2s,i2e)                                     9d15s17
              i10=i1s                                                   9d15s17
              i1n=noc(isd3)                                              9d15s17
              nrowi=noc(isao)*noc(isd2)                                 12d25s17
              i4o=ioooo(isi)                                            9d15s17
              do i2=i2s,i2e                                             9d15s17
               i2m=i2-1                                                 9d15s17
               if(i2.eq.i2e)i1n=i1e                                     9d15s17
               do i1=i10,i1n                                            9d15s17
                if(i2.le.idoub(isbo).and.i1.gt.idoub(isd3))then          9d15s17
                 i1m=i1-1-idoub(isd3)                                    9d15s17
                 do i6=0,iacto(isbv)-1                                  9d15s17
                  in=min(i1m,i6)                                        9d15s17
                  ix=max(i1m,i6)                                        9d15s17
                  ieq=((ix*(ix+1))/2)+in                                9d15s17
                  inot=i1m+iacto(isd3)*i6                               9d15s17
                  iswitch=iabs(isbv-isd3)                               9d15s17
                  iswitch=iswitch/max(1,iswitch)                        9d15s17
                  iden=jdenpt(isd)+nrowd*((inot-ieq)*iswitch+ieq)       9d15s17
                  do i4=0,iacto(isd2)-1                                 11d6s17
                   i4p=i4+idoub(isd2)                                   11d6s17
                   do i5=0,iacto(isav)-1                                9d15s17
                    jhessa=ihessa(isa,isb)+i2m+nrowa*(i5+iacto(isav)*i6)9d15s17
                    in=min(i5,i4)                                       9d15s17
                    ix=max(i5,i4)                                       9d15s17
                    ieq=((ix*(ix+1))/2)+in                              9d15s17
                    inot=i5+iacto(isav)*i4                              11d6s17
                    iswitch=iabs(isav-isd2)                             11d6s17
                    iswitch=iswitch/max(1,iswitch)                      9d15s17
                    iadd=iden+(inot-ieq)*iswitch+ieq                    9d15s17
                    do i3=0,idoub(isao)-1                               9d15s17
                     inot=i4p+noc(isd2)*i3                              11d6s17
                     fact=bc(i3+noc(isao)*i4p+i4o)*fmult                12d25s17
                     khessa=jhessa+idoub(isbo)*i3                       9d15s17
                     bc(khessa)=bc(khessa)+fact*bc(iadd)                9d15s17
                    end do                                              9d15s17
                   end do                                               9d15s17
                  end do                                                9d15s17
                 end do                                                 9d15s17
                end if                                                  9d15s17
                i4o=i4o+nrowi                                           9d15s17
               end do                                                   9d15s17
               i10=1                                                    9d15s17
              end do                                                    9d15s17
             end if                                                     9d15s17
            end do                                                      9d15s17
            if(igot.eq.0)then
               write(6,*)('fail no. 5')
             write(6,*)('for k:')
             write(6,*)('failed to match integral symmetry type ')
             do isi=1,nsdlk
              write(6,12)isi,(isblk(j,isi),j=1,4)
             end do
             call dws_sync
             call dws_finalize
             stop
            end if
           else if(isblk(2,isd).eq.isav.and.isblk(3,isd).eq.isbv.and.   11d6s17
     $          jdenpt(isd).gt.0)then                                   9d15s17
            isd1=isblk(1,isd)
            isd4=isblk(4,isd)                                           11d6s17
            if(isav.eq.isd1)then                                        9d15s17
             nrowd=(iacto(isav)*(iacto(isav)+1))/2                      9d15s17
            else                                                        9d15s17
             nrowd=iacto(isav)*iacto(isd1)                              9d15s17
            end if                                                      9d15s17
            igot=0                                                      9d15s17
            do isi=1,nsdlk                                              9d15s17
             if(isblk(1,isi).eq.isao.and.isblk(3,isi).eq.isbo.and.      9d15s17
     $          isblk(2,isi).eq.isd1.and.                               11d6s17
     $            isblk(4,isi).eq.isd4)then                             11d6s17
              igot=1                                                    9d15s17
              call ilimts(noc(isbo),noc(isd4),mynprocg,mynowprog,il,ih, 11d6s17
     $             i1s,i1e,i2s,i2e)                                     9d15s17
              i10=i1s                                                   9d15s17
              i1n=noc(isbo)                                             9d15s17
              if(isao.eq.isd1)then                                       9d15s17
               nrowi=(noc(isao)*(noc(isao)+1))/2                        9d15s17
              else                                                      9d15s17
               nrowi=noc(isao)*noc(isd1)                                 9d15s17
              end if                                                    9d15s17
              i4o=ioooo(isi)                                            9d15s17
              do i2=i2s,i2e                                             9d15s17
               i2m=i2-1-idoub(isd4)                                     11d6s17
               if(i2.eq.i2e)i1n=i1e                                     9d15s17
               do i1=i10,i1n                                            9d15s17
                if(i1.le.idoub(isbo).and.i2.gt.idoub(isd4))then         11d6s17
                 i1m=i1-1                                               9d15s17
                 do i6=0,iacto(isbv)-1                                  9d15s17
                  in=min(i2m,i6)                                        9d15s17
                  ix=max(i2m,i6)                                        9d15s17
                  ieq=((ix*(ix+1))/2)+in                                9d15s17
                  inot=i6+iacto(isbv)*i2m                               11d6s17
                  iswitch=iabs(isd4-isbv)                               11d6s17
                  iswitch=iswitch/max(1,iswitch)                        9d15s17
                  iden=nrowd*((inot-ieq)*iswitch+ieq)+jdenpt(isd)       9d15s17
                  do i4=0,iacto(isd1)-1                                  9d15s17
                   i4p=i4+idoub(isd1)                                    9d15s17
                   do i5=0,iacto(isav)-1                                9d15s17
                    jhessa=ihessa(isa,isb)+nrowa*(i5+iacto(isav)*i6)+i1m9d15s17
                    in=min(i5,i4)                                       9d15s17
                    ix=max(i5,i4)                                       9d15s17
                    ieq=((ix*(ix+1))/2)+in                              9d15s17
                    inot=i4+iacto(isd1)*i5                              9d15s17
                    iswitch=iabs(isd1-isav)                             9d15s17
                    iswitch=iswitch/max(1,iswitch)                      9d15s17
                    iadd=iden+(inot-ieq)*iswitch+ieq                    9d15s17
                    do i3=0,idoub(isao)-1                               9d15s17
                     in=min(i3,i4p)                                     9d15s17
                     ix=max(i3,i4p)                                     9d15s17
                     ieq=((ix*(ix+1))/2)+in                             9d15s17
                     inot=i3+noc(isao)*i4p                              9d15s17
                     iswitch=iabs(isao-isd1)                            9d15s17
                     iswitch=iswitch/max(1,iswitch)                     9d15s17
                     iad=(inot-ieq)*iswitch+ieq                         9d15s17
                     fact=bc(iad+i4o)*fmult                             9d15s17
                     khessa=jhessa+idoub(isbo)*i3                       9d15s17
                     bc(khessa)=bc(khessa)+fact*bc(iadd)                9d15s17
                    end do                                              9d15s17
                   end do                                               9d15s17
                  end do                                                9d15s17
                 end do                                                 9d15s17
                end if                                                  9d15s17
                i4o=i4o+nrowi                                           9d15s17
               end do                                                   9d15s17
               i10=1                                                    9d15s17
              end do                                                    9d15s17
             end if                                                     9d15s17
             if(isblk(2,isi).eq.isao.and.isblk(4,isi).eq.isbo.and.      9d15s17
     $          isblk(1,isi).eq.isd1.and.                               11d6s17
     $            isblk(3,isi).eq.isd4.and.igot.eq.0)then               11d6s17
              igot=1                                                    9d15s17
              call ilimts(noc(isd4),noc(isbo),mynprocg,mynowprog,il,ih, 11d6s17
     $             i1s,i1e,i2s,i2e)                                     9d15s17
              i10=i1s                                                   9d15s17
              i1n=noc(isd4)                                             11d6s17
              if(isao.eq.isd1)then                                       9d15s17
               nrowi=(noc(isao)*(noc(isao)+1))/2                        9d15s17
              else                                                      9d15s17
               nrowi=noc(isao)*noc(isd1)                                 9d15s17
              end if                                                    9d15s17
              i4o=ioooo(isi)                                            9d15s17
              do i2=i2s,i2e                                             9d15s17
               i2m=i2-1                                                 9d15s17
               if(i2.eq.i2e)i1n=i1e                                     9d15s17
               do i1=i10,i1n                                            9d15s17
                if(i2.le.idoub(isbo).and.i1.gt.idoub(isd4))then         11d6s17
                 i1m=i1-1-idoub(isd4)                                   11d6s17
                 do i6=0,iacto(isbv)-1                                  9d15s17
                  in=min(i1m,i6)                                        9d15s17
                  ix=max(i1m,i6)                                        9d15s17
                  ieq=((ix*(ix+1))/2)+in                                9d15s17
                  inot=i6+iacto(isbv)*i1m                               11d6s17
                  iswitch=iabs(isbv-isd4)                               11d6s17
                  iswitch=iswitch/max(1,iswitch)                        9d15s17
                  iden=jdenpt(isd)+nrowd*((inot-ieq)*iswitch+ieq)       9d15s17
                  do i4=0,iacto(isd1)-1                                  9d15s17
                   i4p=i4+idoub(isd1)                                    9d15s17
                   do i5=0,iacto(isav)-1                                9d15s17
                    jhessa=ihessa(isa,isb)+i2m+nrowa*(i5+iacto(isav)*i6)9d15s17
                    in=min(i5,i4)                                       9d15s17
                    ix=max(i5,i4)                                       9d15s17
                    ieq=((ix*(ix+1))/2)+in                              9d15s17
                    inot=i4+iacto(isd1)*i5                              9d15s17
                    iswitch=iabs(isav-isd1)                             9d15s17
                    iswitch=iswitch/max(1,iswitch)                      9d15s17
                    iadd=iden+(inot-ieq)*iswitch+ieq                    9d15s17
                    do i3=0,idoub(isao)-1                               9d15s17
                     in=min(i3,i4p)                                     9d15s17
                     ix=max(i3,i4p)                                     9d15s17
                     ieq=((ix*(ix+1))/2)+in                             9d15s17
                     inot=i4p+noc(isd1)*i3                               9d15s17
                     iswitch=iabs(isao-isd1)                            9d15s17
                     iswitch=iswitch/max(1,iswitch)                     9d15s17
                     iad=(inot-ieq)*iswitch+ieq                         9d15s17
                     fact=bc(iad+i4o)*fmult                             9d15s17
                     khessa=jhessa+idoub(isbo)*i3                       9d15s17
                     bc(khessa)=bc(khessa)+fact*bc(iadd)                9d15s17
                    end do                                              9d15s17
                   end do                                               9d15s17
                  end do                                                9d15s17
                 end do                                                 9d15s17
                end if                                                  9d15s17
                i4o=i4o+nrowi                                           9d15s17
               end do                                                   9d15s17
               i10=1                                                    9d15s17
              end do                                                    9d15s17
             end if                                                     9d15s17
             if(isblk(2,isi).eq.isao.and.isblk(3,isi).eq.isbo.and.      11d6s17
     $          isblk(1,isi).eq.isd1.and.                               11d6s17
     $            isblk(4,isi).eq.isd4.and.igot.eq.0)then               11d6s17
              igot=1                                                    9d15s17
              call ilimts(noc(isbo),noc(isd4),mynprocg,mynowprog,il,ih, 11d6s17
     $             i1s,i1e,i2s,i2e)                                     9d15s17
              i10=i1s                                                   9d15s17
              i1n=noc(isbo)                                             11d6s17
              if(isao.eq.isd1)then                                       9d15s17
               nrowi=(noc(isao)*(noc(isao)+1))/2                        9d15s17
              else                                                      9d15s17
               nrowi=noc(isao)*noc(isd1)                                 9d15s17
              end if                                                    9d15s17
              i4o=ioooo(isi)                                            9d15s17
              do i2=i2s,i2e                                             9d15s17
               i2m=i2-1-idoub(isd4)                                     11d6s17
               if(i2.eq.i2e)i1n=i1e                                     9d15s17
               do i1=i10,i1n                                            9d15s17
                if(i1.le.idoub(isbo).and.i2.gt.idoub(isd4))then         11d6s17
                 i1m=i1-1                                               11d6s17
                 do i6=0,iacto(isbv)-1                                  9d15s17
                  in=min(i2m,i6)                                        11d6s17
                  ix=max(i2m,i6)                                        11d6s17
                  ieq=((ix*(ix+1))/2)+in                                9d15s17
                  inot=i6+iacto(isbv)*i2m                               11d6s17
                  iswitch=iabs(isbv-isd4)                               11d6s17
                  iswitch=iswitch/max(1,iswitch)                        9d15s17
                  iden=jdenpt(isd)+nrowd*((inot-ieq)*iswitch+ieq)       9d15s17
                  do i4=0,iacto(isd1)-1                                  9d15s17
                   i4p=i4+idoub(isd1)                                    9d15s17
                   do i5=0,iacto(isav)-1                                9d15s17
                    jhessa=ihessa(isa,isb)+i1m+nrowa*(i5+iacto(isav)*i6)11d6s17
                    in=min(i5,i4)                                       9d15s17
                    ix=max(i5,i4)                                       9d15s17
                    ieq=((ix*(ix+1))/2)+in                              9d15s17
                    inot=i4+iacto(isd1)*i5                              9d15s17
                    iswitch=iabs(isav-isd1)                             9d15s17
                    iswitch=iswitch/max(1,iswitch)                      9d15s17
                    iadd=iden+(inot-ieq)*iswitch+ieq                    9d15s17
                    do i3=0,idoub(isao)-1                               9d15s17
                     in=min(i3,i4p)                                     9d15s17
                     ix=max(i3,i4p)                                     9d15s17
                     ieq=((ix*(ix+1))/2)+in                             9d15s17
                     inot=i4p+noc(isd1)*i3                               9d15s17
                     iswitch=iabs(isao-isd1)                            9d15s17
                     iswitch=iswitch/max(1,iswitch)                     9d15s17
                     iad=(inot-ieq)*iswitch+ieq                         9d15s17
                     fact=bc(iad+i4o)*fmult                             9d15s17
                     khessa=jhessa+idoub(isbo)*i3                       9d15s17
                     bc(khessa)=bc(khessa)+fact*bc(iadd)                9d15s17
                    end do                                              9d15s17
                   end do                                               9d15s17
                  end do                                                9d15s17
                 end do                                                 9d15s17
                end if                                                  9d15s17
                i4o=i4o+nrowi                                           9d15s17
               end do                                                   9d15s17
               i10=1                                                    9d15s17
              end do                                                    9d15s17
             end if                                                     9d15s17
            end do                                                      9d15s17
            if(igot.eq.0)then
             write(6,*)('for k4:')
             write(6,*)('failed to match integral symmetry type ')
             write(6,12)isd,(isblk(j,isd),j=1,4)
             do isi=1,nsdlk
              write(6,12)isi,(isblk(j,isi),j=1,4)
             end do
             call dws_sync
             call dws_finalize
             stop
            end if
           end if                                                       9d15s17
           end if                                                       10d27s17
          end do                                                        9d15s17
          end if                                                        9d25s17
         end if                                                         9d15s17
c
c     doub-active contributions to one particle density terms
c
         if(idoub(isao)*idoub(isbo).ne.0.and.lsym)then                  11d30s17
          igot=0                                                         7d20s17
          igot2=0                                                       7d20s17
          igot3=0                                                       7d21s17
          igot4=0
          igot5=0                                                       7d21s17
          igot6=0                                                       7d21s17
          do is=1,nsdlk                                                  7d20s17
c
c     ddaa part k
c
           if(idoit(4).ne.0)then                                        9d25s17
           if(isblk(1,is).eq.isa.and.isblk(2,is).eq.isao.and.
     $          isblk(3,is).eq.isbo)then                                7d20s17
            igot=1                                                      7d20s17
            is4=isblk(4,is)                                             7d20s17
            call ilimts(noc(isbo),noc(is4),mynprocg,mynowprog,ilo,iho,  7d20s17
     $           i1so,i1eo,i2so,i2eo)                                   7d20s17
            i10=i1so                                                    7d20s17
            i1n=noc(isbo)                                               7d20s17
            i4o=ioooo(is)                                               7d20s17
            if(isa.eq.isao)then                                         7d20s17
             nrowo=(noc(isa)*(noc(isa)+1))/2                            7d20s17
            else                                                        7d20s17
             nrowo=noc(isa)*noc(isao)                                   7d20s17
            end if                                                      7d20s17
            nnn=idoub(isbo)*idoub(isao)*iacto(isa)                      7d20s17
            do i2=i2so,i2eo                                             7d20s17
             if(i2.eq.i2eo)i1n=i1eo                                     7d20s17
             i2m=i2-1-idoub(is4)                                        7d20s17
             do i1=i10,i1n                                              7d20s17
              if(i2.gt.idoub(is4).and.i1.le.idoub(isbo))then            7d20s17
               i1m=i1-1                                                 7d20s17
               if(isa.eq.isao)then                                      7d20s17
                do i4=0,idoub(isao)-1                                    7d20s17
                 do i3=0,iacto(isa)-1                                   7d20s17
                  i3p=i3+idoub(isa)                                     7d20s17
                  ix=max(i3p,i4)                                        8d8s17
                  in=min(i3p,i4)                                        8d8s17
                  i4ou=i4o+((ix*(ix+1))/2)+in                           8d8s17
                  term=8d0*bc(i4ou)                                     7d20s17
                  iadd=iden1(isb)+iacto(isb)*i2m                        7d20s17
                  iadh=ihessa(isa,isb)+i1m+idoub(isbo)*(i4+idoub(isao)  7d20s17
     $                 *i3)                                             7d20s17
                  do i5=0,iacto(isb)-1                                  7d20s17
                   jadh=iadh+i5*nnn                                     7d20s17
                   bc(jadh)=bc(jadh)-term*bc(iadd+i5)                   7d20s17
                  end do                                                7d20s17
                 end do                                                 7d20s17
                end do                                                   7d20s17
               else                                                     7d20s17
                do i4=0,idoub(isao)-1                                    7d20s17
                 do i3=0,iacto(isa)-1                                   7d20s17
                  i3p=i3+idoub(isa)                                     7d20s17
                  i4ou=i4o+i3p+noc(isa)*i4                              7d20s17
                  term=8d0*bc(i4ou)                                     7d20s17
                  iadd=iden1(isb)+iacto(isb)*i2m                        7d20s17
                  iadh=ihessa(isa,isb)+i1m+idoub(isbo)*(i4+idoub(isao)  7d20s17
     $                 *i3)                                             7d20s17
                  do i5=0,iacto(isb)-1                                  7d20s17
                   jadh=iadh+i5*nnn                                     7d20s17
                   bc(jadh)=bc(jadh)-term*bc(iadd+i5)                   7d20s17
                  end do                                                7d20s17
                 end do                                                 7d20s17
                end do                                                   7d20s17
               end if                                                   7d20s17
              end if                                                    7d20s17
              i4o=i4o+nrowo                                             7d20s17
             end do                                                     7d20s17
             i10=1                                                      7d20s17
            end do                                                      7d20s17
           end if                                                       7d20s17
           if(isblk(1,is).eq.isao.and.isblk(2,is).eq.isa.and.           7d20s17
     $          isblk(3,is).eq.isbo.and.igot.eq.0)then                  7d20s17
            igot=1                                                      7d20s17
            is4=isblk(4,is)                                             7d20s17
            call ilimts(noc(isbo),noc(is4),mynprocg,mynowprog,ilo,iho,  7d20s17
     $           i1so,i1eo,i2so,i2eo)                                   7d20s17
            i10=i1so                                                    7d20s17
            i1n=noc(isbo)                                               7d20s17
            i4o=ioooo(is)                                               7d20s17
            nrowo=noc(isa)*noc(isao)                                    7d20s17
            nnn=idoub(isbo)*idoub(isao)*iacto(isa)                      7d20s17
            do i2=i2so,i2eo                                             7d20s17
             if(i2.eq.i2eo)i1n=i1eo                                     7d20s17
             i2m=i2-1-idoub(is4)                                        7d20s17
             do i1=i10,i1n                                              7d20s17
              if(i2.gt.idoub(is4).and.i1.le.idoub(isbo))then            7d20s17
              i1m=i1-1                                                  7d20s17
               do i4=0,iacto(isa)-1                                     7d20s17
                i4p=i4+idoub(isa)                                       7d20s17
                do i3=0,idoub(isao)-1                                   7d20s17
                 i4ou=i4o+i3+noc(isao)*i4p                              7d21s17
                 term=8d0*bc(i4ou)                                      7d20s17
                 iadd=iden1(isb)+iacto(isb)*i2m                         7d20s17
                 iadh=ihessa(isa,isb)+i1m+idoub(isbo)*(i3+idoub(isao)   11d3s17
     $                *i4)                                              7d20s17
                 do i5=0,iacto(isb)-1                                   7d20s17
                  jadh=iadh+i5*nnn                                      7d20s17
                  bc(jadh)=bc(jadh)-term*bc(iadd+i5)                    7d20s17
                 end do                                                 7d20s17
                end do                                                  7d20s17
               end do                                                    7d20s17
              end if                                                    7d20s17
              i4o=i4o+nrowo                                             7d20s17
             end do                                                     7d20s17
             i10=1                                                      7d20s17
            end do                                                      7d20s17
           end if                                                       7d20s17
           if(isblk(1,is).eq.isa.and.isblk(2,is).eq.isao.and.
     $          isblk(4,is).eq.isbo.and.igot.eq.0)then                  7d20s17
            igot=1                                                      7d20s17
            is3=isblk(3,is)                                             7d20s17
            call ilimts(noc(is3),noc(isbo),mynprocg,mynowprog,ilo,iho,  7d20s17
     $           i1so,i1eo,i2so,i2eo)                                   7d20s17
            i10=i1so                                                    7d20s17
            i1n=noc(is3)                                                7d20s17
            i4o=ioooo(is)                                               7d20s17
            nrowo=noc(isa)*noc(isao)                                    7d20s17
            nnn=idoub(isbo)*idoub(isao)*iacto(isa)                      7d20s17
            do i2=i2so,i2eo                                             7d20s17
             if(i2.eq.i2eo)i1n=i1eo                                     7d20s17
             i2m=i2-1                                                   7d20s17
             do i1=i10,i1n                                              7d20s17
              if(i1.gt.idoub(is3).and.i2.le.idoub(isbo))then            7d20s17
               i1m=i1-1-idoub(is3)                                      7d20s17
               do i4=0,idoub(isao)-1                                    7d20s17
                do i3=0,iacto(isa)-1                                    7d20s17
                 i3p=i3+idoub(isa)                                      7d20s17
                 i4ou=i4o+i3p+noc(isa)*i4                               7d20s17
                 term=8d0*bc(i4ou)                                      7d20s17
                 iadd=iden1(isb)+iacto(isb)*i1m                         7d20s17
                 iadh=ihessa(isa,isb)+i2m+idoub(isbo)*(i4+idoub(isao)   7d20s17
     $                *i3)                                              7d20s17
                 do i5=0,iacto(isb)-1                                   7d20s17
                  jadh=iadh+i5*nnn                                      7d20s17
                  bc(jadh)=bc(jadh)-term*bc(iadd+i5)                    7d20s17
                 end do                                                 7d20s17
                end do                                                  7d20s17
               end do                                                    7d20s17
              end if                                                    7d20s17
              i4o=i4o+nrowo                                             7d20s17
             end do                                                     7d20s17
             i10=1                                                      7d20s17
            end do                                                      7d20s17
           end if                                                       7d20s17
           if(isblk(1,is).eq.isao.and.isblk(2,is).eq.isa.and.           7d20s17
     $          isblk(4,is).eq.isbo.and.igot.eq.0)then                  7d20s17
            igot=1                                                      7d20s17
            is3=isblk(3,is)                                             7d20s17
            call ilimts(noc(is3),noc(isbo),mynprocg,mynowprog,ilo,iho,  7d20s17
     $           i1so,i1eo,i2so,i2eo)                                   7d20s17
            i10=i1so                                                    7d20s17
            i1n=noc(is3)                                                7d20s17
            i4o=ioooo(is)                                               7d20s17
            nrowo=noc(isa)*noc(isao)                                    7d20s17
            nnn=idoub(isbo)*idoub(isao)*iacto(isa)                      7d20s17
            do i2=i2so,i2eo                                             7d20s17
             if(i2.eq.i2eo)i1n=i1eo                                     7d20s17
             i2m=i2-1                                                   7d20s17
             do i1=i10,i1n                                              7d20s17
              if(i1.gt.idoub(is3).and.i2.le.idoub(isbo))then            7d20s17
               i1m=i1-1-idoub(is3)                                      7d20s17
               do i4=0,iacto(isa)-1                                     7d20s17
                i4p=i4+idoub(isa)                                       7d20s17
                do i3=0,idoub(isao)-1                                   7d20s17
                 i4ou=i4o+i3+noc(isao)*i4p                              7d20s17
                 term=8d0*bc(i4ou)                                      7d20s17
                 iadd=iden1(isb)+iacto(isb)*i1m                         7d20s17
                 iadh=ihessa(isa,isb)+i2m+idoub(isbo)*(i3+idoub(isao)   7d20s17
     $                *i4)                                              7d20s17
                 do i5=0,iacto(isb)-1                                   7d20s17
                  jadh=iadh+i5*nnn                                      7d20s17
                  bc(jadh)=bc(jadh)-term*bc(iadd+i5)                    7d20s17
                 end do                                                 7d20s17
                end do                                                  7d20s17
               end do                                                    7d20s17
              end if                                                    7d20s17
              i4o=i4o+nrowo                                             7d20s17
             end do                                                     7d20s17
             i10=1                                                      7d20s17
            end do                                                      7d20s17
           end if                                                       7d20s17
c
c     ddaa part l
c
           if(isblk(1,is).eq.isbo.and.isblk(2,is).eq.isa.and.           7d20s17
     $          isblk(4,is).eq.isao)then                                7d20s17
            is3=isblk(3,is)                                             7d20s17
            igot2=1                                                     7d20s17
            call ilimts(noc(is3),noc(isao),mynprocg,mynowprog,ilo,iho,  7d20s17
     $           i1so,i1eo,i2so,i2eo)                                   7d20s17
            i10=i1so                                                    7d20s17
            i1n=noc(is3)                                                7d20s17
            i4o=ioooo(is)                                               7d20s17
            if(isbo.eq.isa)then                                         7d20s17
             nrowo=((noc(isa)*(noc(isa)+1))/2)                          7d20s17
            else                                                        7d20s17
             nrowo=noc(isa)*noc(isbo)                                   7d20s17
            end if                                                      7d20s17
            nnn=idoub(isbo)*idoub(isao)*iacto(isa)                      7d20s17
            do i2=i2so,i2eo                                             7d20s17
             if(i2.eq.i2eo)i1n=i1eo                                     7d20s17
             i2m=i2-1
             do i1=i10,i1n                                              7d20s17
              if(i1.gt.idoub(is3).and.i2.le.idoub(isao))then            7d20s17
               i1m=i1-1-idoub(is3)                                       7d20s17
               if(isbo.eq.isa)then                                      7d20s17
                do i4=0,iacto(isa)-1                                     7d20s17
                 i4p=i4+idoub(isa)                                      7d20s17
                 do i3=0,idoub(isbo)-1                                  7d20s17
                  i4ou=i4o+((i4p*(i4p+1))/2)+i3                         7d20s17
                  term=2d0*bc(i4ou)                                     7d20s17
                  iadd=iden1(isb)+iacto(isb)*i1m                        7d20s17
                  iadh=ihessa(isa,isb)+i3+idoub(isbo)*(i2m+idoub(isao)  7d20s17
     $                 *i4)                                             7d20s17
                  do i5=0,iacto(isb)-1                                  7d20s17
                   jadh=iadh+nnn*i5                                     7d20s17
                   bc(jadh)=bc(jadh)+bc(iadd+i5)*term                   7d20s17
                  end do                                                7d20s17
                 end do                                                 7d20s17
                end do                                                  7d20s17
               else                                                     7d20s17
                do i4=0,iacto(isa)-1                                     7d20s17
                 i4p=i4+idoub(isa)                                      7d20s17
                 do i3=0,idoub(isbo)-1                                  7d20s17
                  i4ou=i4o+i3+noc(isbo)*i4p                             7d20s17
                  term=2d0*bc(i4ou)                                     7d20s17
                  iadd=iden1(isb)+iacto(isb)*i1m                        7d20s17
                  iadh=ihessa(isa,isb)+i3+idoub(isbo)*(i2m+idoub(isao)  7d20s17
     $                 *i4)                                             7d20s17
                  do i5=0,iacto(isb)-1                                  7d20s17
                   jadh=iadh+nnn*i5                                     7d20s17
                   bc(jadh)=bc(jadh)+bc(iadd+i5)*term                   7d20s17
                  end do                                                7d20s17
                 end do                                                 7d20s17
                end do                                                  7d20s17
               end if                                                   7d20s17
              end if                                                    7d20s17
              i4o=i4o+nrowo                                             7d20s17
             end do                                                     7d20s17
             i10=1                                                      7d20s17
            end do                                                      7d20s17
           end if                                                       7d20s17
           if(isblk(2,is).eq.isbo.and.isblk(1,is).eq.isa.and.           7d20s17
     $          isblk(4,is).eq.isao.and.igot2.eq.0)then                 7d20s17
            is3=isblk(3,is)                                             7d20s17
            igot2=1                                                     7d20s17
            call ilimts(noc(is3),noc(isao),mynprocg,mynowprog,ilo,iho,  7d20s17
     $           i1so,i1eo,i2so,i2eo)                                   7d20s17
            i10=i1so                                                    7d20s17
            i1n=noc(is3)                                                7d20s17
            i4o=ioooo(is)                                               7d20s17
            nrowo=noc(isa)*noc(isbo)                                    7d20s17
            nnn=idoub(isbo)*idoub(isao)*iacto(isa)                      7d20s17
            do i2=i2so,i2eo                                             7d20s17
             if(i2.eq.i2eo)i1n=i1eo                                     7d20s17
             i2m=i2-1
             do i1=i10,i1n                                              7d20s17
              if(i1.gt.idoub(is3).and.i2.le.idoub(isao))then            7d20s17
               i1m=i1-1-idoub(is3)                                       7d20s17
               do i4=0,idoub(isbo)-1                                    7d20s17
                do i3=0,iacto(isa)-1                                    7d20s17
                 i3p=i3+idoub(isa)                                      7d20s17
                 i4ou=i4o+i3p+noc(isa)*i4                               7d20s17
                 term=2d0*bc(i4ou)                                      7d20s17
                 iadd=iden1(isb)+iacto(isb)*i1m                         7d20s17
                 iadh=ihessa(isa,isb)+i4+idoub(isbo)*(i2m+idoub(isao)   7d20s17
     $                *i3)                                              7d20s17
                 do i5=0,iacto(isb)-1                                   7d20s17
                  jadh=iadh+nnn*i5                                      7d20s17
                  bc(jadh)=bc(jadh)+bc(iadd+i5)*term                    7d20s17
                 end do                                                 7d20s17
                end do                                                  7d20s17
               end do                                                   7d20s17
              end if                                                     7d20s17
              i4o=i4o+nrowo                                             7d20s17
             end do                                                     7d20s17
             i10=1                                                      7d20s17
            end do                                                      7d20s17
           end if                                                       7d20s17
           if(isblk(2,is).eq.isbo.and.isblk(1,is).eq.isa.and.           7d20s17
     $          isblk(3,is).eq.isao.and.igot2.eq.0)then                 7d20s17
            is4=isblk(4,is)                                             7d20s17
            igot2=1                                                     7d20s17
            call ilimts(noc(isao),noc(is4),mynprocg,mynowprog,ilo,iho,  7d20s17
     $           i1so,i1eo,i2so,i2eo)                                   7d20s17
            i10=i1so                                                    7d20s17
            i1n=noc(isao)                                               7d20s17
            i4o=ioooo(is)                                               7d20s17
            nrowo=noc(isa)*noc(isbo)                                    7d20s17
            nnn=idoub(isbo)*idoub(isao)*iacto(isa)                      7d20s17
            do i2=i2so,i2eo                                             7d20s17
             if(i2.eq.i2eo)i1n=i1eo                                     7d20s17
             i2m=i2-1-idoub(is4)                                        7d20s17
             do i1=i10,i1n                                              7d20s17
              if(i2.gt.idoub(is4).and.i1.le.idoub(isao))then            7d20s17
               i1m=i1-1                                                 7d20s17
               do i4=0,idoub(isbo)-1                                    7d20s17
                do i3=0,iacto(isa)-1                                    7d20s17
                 i3p=i3+idoub(isa)                                      7d20s17
                 i4ou=i4o+i3p+noc(isa)*i4                               7d20s17
                 term=2d0*bc(i4ou)                                      7d20s17
                 iadd=iden1(isb)+iacto(isb)*i2m                         7d20s17
                 iadh=ihessa(isa,isb)+i4+idoub(isbo)*(i1m+idoub(isao)   7d20s17
     $                *i3)                                              7d20s17
                 do i5=0,iacto(isb)-1                                   7d20s17
                  jadh=iadh+nnn*i5                                      7d20s17
                  bc(jadh)=bc(jadh)+bc(iadd+i5)*term                    7d20s17
                 end do                                                 7d20s17
                end do                                                  7d20s17
               end do                                                   7d20s17
              end if                                                     7d20s17
              i4o=i4o+nrowo                                             7d20s17
             end do                                                     7d20s17
             i10=1                                                      7d20s17
            end do                                                      7d20s17
           end if                                                       7d20s17
           if(isblk(1,is).eq.isbo.and.isblk(2,is).eq.isa.and.           7d20s17
     $          isblk(3,is).eq.isao.and.igot2.eq.0)then                 7d20s17
            is4=isblk(4,is)                                             7d20s17
            igot2=1                                                     7d20s17
            call ilimts(noc(isao),noc(is4),mynprocg,mynowprog,ilo,iho,  7d20s17
     $           i1so,i1eo,i2so,i2eo)                                   7d20s17
            i10=i1so                                                    7d20s17
            i1n=noc(isao)                                               7d20s17
            i4o=ioooo(is)                                               7d20s17
            nrowo=noc(isa)*noc(isbo)                                    7d20s17
            nnn=idoub(isbo)*idoub(isao)*iacto(isa)                      7d20s17
            do i2=i2so,i2eo                                             7d20s17
             if(i2.eq.i2eo)i1n=i1eo                                     7d20s17
             i2m=i2-1-idoub(is4)                                        7d20s17
             do i1=i10,i1n                                              7d20s17
              if(i2.gt.idoub(is4).and.i1.le.idoub(isao))then            7d20s17
               i1m=i1-1                                                 7d20s17
               do i4=0,iacto(isa)-1                                     7d20s17
                i4p=i4+idoub(isa)                                       7d20s17
                do i3=0,idoub(isbo)-1                                   7d20s17
                 i4ou=i4o+i3+noc(isbo)*i4p                              7d20s17
                 term=2d0*bc(i4ou)                                      7d20s17
                 iadd=iden1(isb)+iacto(isb)*i2m                         7d20s17
                 iadh=ihessa(isa,isb)+i3+idoub(isbo)*(i1m+idoub(isao)   7d20s17
     $                *i4)                                              7d20s17
                 do i5=0,iacto(isb)-1                                   7d20s17
                  jadh=iadh+nnn*i5                                      7d20s17
                  bc(jadh)=bc(jadh)+bc(iadd+i5)*term                    7d20s17
                 end do                                                 7d20s17
                end do                                                  7d20s17
               end do                                                   7d20s17
              end if                                                     7d20s17
              i4o=i4o+nrowo                                             7d20s17
             end do                                                     7d20s17
             i10=1                                                      7d20s17
            end do                                                      7d20s17
           end if                                                       7d20s17
c
c     ddaa part m
c
           if(isblk(1,is).eq.isbo.and.isblk(2,is).eq.isao.and.          7d21s17
     $          isblk(4,is).eq.isa)then                                 7d21s17
            is3=isblk(3,is)                                             7d20s17
            igot3=1                                                     7d20s17
            call ilimts(noc(is3),noc(isa),mynprocg,mynowprog,ilo,iho,   7d21s17
     $           i1so,i1eo,i2so,i2eo)                                   7d20s17
            i10=i1so                                                    7d20s17
            i1n=noc(is3)                                                7d20s17
            i4o=ioooo(is)                                               7d20s17
            if(isbo.eq.isao)then                                        7d21s17
             nrowo=(noc(isbo)*(noc(isbo)+1))/2                          7d21s17
            else                                                        7d21s17
             nrowo=noc(isbo)*noc(isao)                                  7d21s17
            end if                                                      7d21s17
            nnn=idoub(isbo)*idoub(isao)*iacto(isa)                      7d21s17
            do i2=i2so,i2eo                                             7d21s17
             if(i2.eq.i2eo)i1n=i1eo                                     7d21s17
             i2m=i2-1-idoub(isa)                                        7d21s17
             do i1=i10,i1n                                              7d21s17
              if(i1.gt.idoub(is3).and.i2.gt.idoub(isa))then             7d21s17
               i1m=i1-1-idoub(is3)                                      7d21s17
               iadd=iden1(isb)+iacto(isb)*i1m                           7d21s17
               if(isbo.eq.isao)then                                     7d21s17
                do i4=0,idoub(isao)-1                                   7d21s17
                 do i3=0,i4-1                                           7d21s17
                  i4ou=i4o+((i4*(i4+1))/2)+i3                           7d21s17
                  term=bc(i4ou)*2d0                                     7d21s17
                  iadh34=ihessa(isa,isb)+i3+idoub(isbo)                 7d21s17
     $                 *(i4+idoub(isao)*i2m)                            7d21s17
                  iadh43=ihessa(isa,isb)+i4+idoub(isbo)                 7d21s17
     $                 *(i3+idoub(isao)*i2m)                            7d21s17
                  do i5=0,iacto(isb)-1                                  7d21s17
                   part=term*bc(iadd+i5)                                7d21s17
                   jadh34=iadh34+i5*nnn                                 7d21s17
                   bc(jadh34)=bc(jadh34)+part                           7d21s17
                   jadh43=iadh43+i5*nnn                                 7d21s17
                   bc(jadh43)=bc(jadh43)+part                           7d21s17
                  end do                                                7d21s17
                 end do                                                 7d21s17
                 i4ou=i4o+((i4*(i4+1))/2)+i4                            7d21s17
                 term=bc(i4ou)*2d0                                      7d21s17
                 iadh=ihessa(isa,isb)+i4+idoub(isbo)                    7d21s17
     $                 *(i4+idoub(isao)*i2m)                            7d21s17
                 do i5=0,iacto(isb)-1                                   7d21s17
                  jadh=iadh+i5*nnn                                      7d21s17
                  bc(jadh)=bc(jadh)+term*bc(iadd+i5)                    7d21s17
                 end do                                                 7d21s17
                end do                                                  7d21s17
               else                                                     7d21s17
                do i4=0,idoub(isao)-1                                   7d21s17
                 do i3=0,idoub(isbo)-1                                  7d21s17
                  i4ou=i4o+i3+noc(isbo)*i4                              7d21s17
                  term=bc(i4ou)*2d0                                     7d21s17
                  iadh=ihessa(isa,isb)+i3+idoub(isbo)*                  7d21s17
     $                 (i4+idoub(isao)*i2m)                             7d21s17
                  do i5=0,iacto(isb)-1                                  7d21s17
                   jadh=iadh+i5*nnn                                     7d21s17
                   bc(jadh)=bc(jadh)+term*bc(iadd+i5)                   7d21s17
                  end do                                                7d21s17
                 end do                                                 7d21s17
                end do                                                  7d21s17
               end if                                                   7d21s17
              end if                                                    7d21s17
              i4o=i4o+nrowo                                             7d21s17
             end do                                                     7d21s17
             i10=1                                                      7d21s17
            end do                                                      7d21s17
           end if                                                       7d21s17
           if(isblk(2,is).eq.isbo.and.isblk(1,is).eq.isao.and.          7d21s17
     $          isblk(4,is).eq.isa.and.igot3.eq.0)then                  7d21s17
            is3=isblk(3,is)                                             7d20s17
            igot3=1                                                     7d20s17
            call ilimts(noc(is3),noc(isa),mynprocg,mynowprog,ilo,iho,   7d21s17
     $           i1so,i1eo,i2so,i2eo)                                   7d20s17
            i10=i1so                                                    7d20s17
            i1n=noc(is3)                                                7d20s17
            i4o=ioooo(is)                                               7d20s17
            nrowo=noc(isbo)*noc(isao)                                   7d21s17
            nnn=idoub(isbo)*idoub(isao)*iacto(isa)                      7d21s17
            do i2=i2so,i2eo                                             7d21s17
             if(i2.eq.i2eo)i1n=i1eo                                     7d21s17
             i2m=i2-1-idoub(isa)                                        7d21s17
             do i1=i10,i1n                                              7d21s17
              if(i1.gt.idoub(is3).and.i2.gt.idoub(isa))then             7d21s17
               i1m=i1-1-idoub(is3)                                      7d21s17
               iadd=iden1(isb)+iacto(isb)*i1m                           7d21s17
               do i4=0,idoub(isbo)-1                                    7d21s17
                do i3=0,idoub(isao)-1                                   7d21s17
                 i4ou=i4o+i3+noc(isao)*i4                               7d21s17
                 term=bc(i4ou)*2d0                                      7d21s17
                 iadh=ihessa(isa,isb)+i4+idoub(isbo)*                   7d21s17
     $                (i3+idoub(isao)*i2m)                              7d21s17
                 do i5=0,iacto(isb)-1                                   7d21s17
                  jadh=iadh+i5*nnn                                      7d21s17
                  bc(jadh)=bc(jadh)+term*bc(iadd+i5)                    7d21s17
                 end do                                                 7d21s17
                end do                                                  7d21s17
               end do                                                   7d21s17
              end if                                                    7d21s17
              i4o=i4o+nrowo                                             7d21s17
             end do                                                     7d21s17
             i10=1                                                      7d21s17
            end do                                                      7d21s17
           end if                                                       7d21s17
           if(isblk(2,is).eq.isbo.and.isblk(1,is).eq.isao.and.          7d21s17
     $          isblk(3,is).eq.isa.and.igot3.eq.0)then                  7d21s17
            is4=isblk(4,is)                                             7d21s17
            igot3=1                                                     7d20s17
            call ilimts(noc(isa),noc(is4),mynprocg,mynowprog,ilo,iho,   7d21s17
     $           i1so,i1eo,i2so,i2eo)                                   7d20s17
            i10=i1so                                                    7d20s17
            i1n=noc(isa)                                                7d21s17
            i4o=ioooo(is)                                               7d20s17
            nrowo=noc(isbo)*noc(isao)                                   7d21s17
            nnn=idoub(isbo)*idoub(isao)*iacto(isa)                      7d21s17
            do i2=i2so,i2eo                                             7d21s17
             if(i2.eq.i2eo)i1n=i1eo                                     7d21s17
             i2m=i2-1-idoub(is4)                                        7d21s17
             do i1=i10,i1n                                              7d21s17
              if(i1.gt.idoub(isa).and.i2.gt.idoub(is4))then             7d21s17
               i1m=i1-1-idoub(isa)                                      7d21s17
               iadd=iden1(isb)+iacto(isb)*i2m                           7d21s17
               do i4=0,idoub(isbo)-1                                    7d21s17
                do i3=0,idoub(isao)-1                                   7d21s17
                 i4ou=i4o+i3+noc(isao)*i4                               7d21s17
                 term=bc(i4ou)*2d0                                      7d21s17
                 iadh=ihessa(isa,isb)+i4+idoub(isbo)*                   7d21s17
     $                (i3+idoub(isao)*i1m)                              7d21s17
                 do i5=0,iacto(isb)-1                                   7d21s17
                  jadh=iadh+i5*nnn                                      7d21s17
                  bc(jadh)=bc(jadh)+term*bc(iadd+i5)                    7d21s17
                 end do                                                 7d21s17
                end do                                                  7d21s17
               end do                                                   7d21s17
              end if                                                    7d21s17
              i4o=i4o+nrowo                                             7d21s17
             end do                                                     7d21s17
             i10=1                                                      7d21s17
            end do                                                      7d21s17
           end if                                                       7d21s17
           if(isblk(1,is).eq.isbo.and.isblk(2,is).eq.isao.and.          7d21s17
     $          isblk(3,is).eq.isa.and.igot3.eq.0)then                  7d21s17
            is4=isblk(4,is)                                             7d21s17
            igot3=1                                                     7d20s17
            call ilimts(noc(isa),noc(is4),mynprocg,mynowprog,ilo,iho,   7d21s17
     $           i1so,i1eo,i2so,i2eo)                                   7d20s17
            i10=i1so                                                    7d20s17
            i1n=noc(isa)                                                7d21s17
            i4o=ioooo(is)                                               7d20s17
            nrowo=noc(isbo)*noc(isao)                                   7d21s17
            nnn=idoub(isbo)*idoub(isao)*iacto(isa)                      7d21s17
            do i2=i2so,i2eo                                             7d21s17
             if(i2.eq.i2eo)i1n=i1eo                                     7d21s17
             i2m=i2-1-idoub(is4)                                        7d21s17
             do i1=i10,i1n                                              7d21s17
              if(i1.gt.idoub(isa).and.i2.gt.idoub(is4))then             7d21s17
               i1m=i1-1-idoub(isa)                                      7d21s17
               iadd=iden1(isb)+iacto(isb)*i2m                           7d21s17
               do i4=0,idoub(isao)-1                                    7d21s17
                do i3=0,idoub(isbo)-1                                   7d21s17
                 i4ou=i4o+i3+noc(isbo)*i4                               7d21s17
                 term=bc(i4ou)*2d0                                      7d21s17
                 iadh=ihessa(isa,isb)+i3+idoub(isbo)*                   7d21s17
     $                (i4+idoub(isao)*i1m)                              7d21s17
                 do i5=0,iacto(isb)-1                                   7d21s17
                  jadh=iadh+i5*nnn                                      7d21s17
                  bc(jadh)=bc(jadh)+term*bc(iadd+i5)                    7d21s17
                 end do                                                 7d21s17
                end do                                                  7d21s17
               end do                                                   7d21s17
              end if                                                    7d21s17
              i4o=i4o+nrowo                                             7d21s17
             end do                                                     7d21s17
             i10=1                                                      7d21s17
            end do                                                      7d21s17
           end if                                                       7d21s17
c
c     ddaa part n
c
           if(isblk(1,is).eq.isb.and.isblk(2,is).eq.isbo.and.           7d21s17
     $          isblk(3,is).eq.isao)then                                7d21s17
            igot4=1                                                     7d21s17
            call ilimts(noc(isao),noc(isa),mynprocg,mynowprog,ilo,iho,  7d21s17
     $           i1so,i1eo,i2so,i2eo)                                   7d20s17
            i10=i1so                                                    7d20s17
            i1n=noc(isao)                                               7d21s17
            i4o=ioooo(is)                                               7d20s17
            if(isb.eq.isbo)then                                         7d21s17
             nrowo=(noc(isb)*(noc(isb)+1))/2                            7d21s17
            else                                                        7d21s17
             nrowo=noc(isbo)*noc(isb)                                   7d21s17
            end if                                                      7d21s17
            nnn=idoub(isbo)*idoub(isao)                                 7d21s17
            do i2=i2so,i2eo                                             7d21s17
             if(i2.eq.i2eo)i1n=i1eo                                     7d21s17
             i2m=i2-1-idoub(isa)                                        7d21s17
             iadd=iden1(isa)+iacto(isa)*i2m                             7d21s17
             do i1=i10,i1n                                              7d21s17
              if(i1.le.idoub(isao).and.i2.gt.idoub(isa))then            7d21s17
               i1m=i1-1                                                 7d21s17
               if(isb.eq.isbo)then                                      7d21s17
                do i4=0,idoub(isbo)-1                                   7d21s17
                 do i3=0,iacto(isb)-1                                   7d21s17
                  i3p=i3+idoub(isb)                                     7d21s17
                  i4ou=i4o+((i3p*(i3p+1))/2)+i4                         7d21s17
                  term=8d0*bc(i4ou)                                     7d21s17
                  iadh=ihessa(isa,isb)+i4+idoub(isbo)                   7d21s17
     $                 *(i1m+idoub(isao)*iacto(isa)*i3)                 7d21s17
                  do i5=0,iacto(isa)-1                                  7d21s17
                   jadh=iadh+nnn*i5                                     7d21s17
                   bc(jadh)=bc(jadh)-term*bc(iadd+i5)                   7d21s17
                  end do                                                7d21s17
                 end do                                                 7d21s17
                end do                                                  7d21s17
               else                                                     7d21s17
                do i4=0,idoub(isbo)-1                                   7d21s17
                 do i3=0,iacto(isb)-1                                   7d21s17
                  i3p=i3+idoub(isb)                                     7d21s17
                  i4ou=i4o+i3p+noc(isb)*i4                              7d21s17
                  term=8d0*bc(i4ou)                                     7d21s17
                  iadh=ihessa(isa,isb)+i4+idoub(isbo)                   7d21s17
     $                 *(i1m+idoub(isao)*iacto(isa)*i3)                 7d21s17
                  do i5=0,iacto(isa)-1                                  7d21s17
                   jadh=iadh+nnn*i5                                     7d21s17
                   bc(jadh)=bc(jadh)-term*bc(iadd+i5)                   7d21s17
                  end do                                                7d21s17
                 end do                                                 7d21s17
                end do                                                  7d21s17
               end if                                                   7d21s17
              end if                                                    7d21s17
              i4o=i4o+nrowo                                             7d21s17
             end do                                                     7d21s17
             i10=1                                                      7d21s17
            end do                                                      7d21s17
           end if                                                       7d21s17
           if(isblk(2,is).eq.isb.and.isblk(1,is).eq.isbo.and.           7d21s17
     $          isblk(3,is).eq.isao.and.igot4.eq.0)then                 7d21s17
            igot4=1                                                     7d21s17
            call ilimts(noc(isao),noc(isa),mynprocg,mynowprog,ilo,iho,  7d21s17
     $           i1so,i1eo,i2so,i2eo)                                   7d20s17
            i10=i1so                                                    7d20s17
            i1n=noc(isao)                                               7d21s17
            i4o=ioooo(is)                                               7d20s17
            nrowo=noc(isbo)*noc(isb)                                    7d21s17
            nnn=idoub(isbo)*idoub(isao)                                 7d21s17
            do i2=i2so,i2eo                                             7d21s17
             if(i2.eq.i2eo)i1n=i1eo                                     7d21s17
             i2m=i2-1-idoub(isa)                                        7d21s17
             iadd=iden1(isa)+iacto(isa)*i2m                             7d21s17
             do i1=i10,i1n                                              7d21s17
              if(i1.le.idoub(isao).and.i2.gt.idoub(isa))then            7d21s17
               i1m=i1-1                                                 7d21s17
               do i4=0,iacto(isb)-1                                     7d21s17
                i4p=i4+idoub(isb)                                       7d21s17
                do i3=0,idoub(isbo)-1                                   7d21s17
                 i4ou=i4o+i3+noc(isbo)*i4p                              7d21s17
                 term=8d0*bc(i4ou)                                      7d21s17
                 iadh=ihessa(isa,isb)+i3+idoub(isbo)                    7d21s17
     $                *(i1m+idoub(isao)*iacto(isa)*i4)                  7d21s17
                 do i5=0,iacto(isa)-1                                   7d21s17
                  jadh=iadh+nnn*i5                                      7d21s17
                  bc(jadh)=bc(jadh)-term*bc(iadd+i5)                    7d21s17
                 end do                                                 7d21s17
                end do                                                  7d21s17
               end do                                                   7d21s17
              end if                                                    7d21s17
              i4o=i4o+nrowo                                             7d21s17
             end do                                                     7d21s17
             i10=1                                                      7d21s17
            end do                                                      7d21s17
           end if                                                       7d21s17
           if(isblk(2,is).eq.isb.and.isblk(1,is).eq.isbo.and.           7d21s17
     $          isblk(4,is).eq.isao.and.igot4.eq.0)then                 7d21s17
            igot4=1                                                     7d21s17
            call ilimts(noc(isa),noc(isao),mynprocg,mynowprog,ilo,iho,  7d21s17
     $           i1so,i1eo,i2so,i2eo)                                   7d20s17
            i10=i1so                                                    7d20s17
            i1n=noc(isa)                                                7d21s17
            i4o=ioooo(is)                                               7d20s17
            nrowo=noc(isbo)*noc(isb)                                    7d21s17
            nnn=idoub(isbo)*idoub(isao)                                 7d21s17
            do i2=i2so,i2eo                                             7d21s17
             if(i2.eq.i2eo)i1n=i1eo                                     7d21s17
             i2m=i2-1
             do i1=i10,i1n                                              7d21s17
              if(i2.le.idoub(isao).and.i1.gt.idoub(isa))then            7d21s17
               i1m=i1-1-idoub(isa)                                      7d21s17
               iadd=iden1(isa)+iacto(isa)*i1m                           7d21s17
               do i4=0,iacto(isb)-1                                     7d21s17
                i4p=i4+idoub(isb)                                       7d21s17
                do i3=0,idoub(isbo)-1                                   7d21s17
                 i4ou=i4o+i3+noc(isbo)*i4p                              7d21s17
                 term=8d0*bc(i4ou)                                      7d21s17
                 iadh=ihessa(isa,isb)+i3+idoub(isbo)                    7d21s17
     $                *(i2m+idoub(isao)*iacto(isa)*i4)                  7d21s17
                 do i5=0,iacto(isa)-1                                   7d21s17
                  jadh=iadh+nnn*i5                                      7d21s17
                  bc(jadh)=bc(jadh)-term*bc(iadd+i5)                    7d21s17
                 end do                                                 7d21s17
                end do                                                  7d21s17
               end do                                                   7d21s17
              end if                                                    7d21s17
              i4o=i4o+nrowo                                             7d21s17
             end do                                                     7d21s17
             i10=1                                                      7d21s17
            end do                                                      7d21s17
           end if                                                       7d21s17
           if(isblk(1,is).eq.isb.and.isblk(2,is).eq.isbo.and.           7d21s17
     $          isblk(4,is).eq.isao.and.igot4.eq.0)then                 7d21s17
            igot4=1                                                     7d21s17
            call ilimts(noc(isa),noc(isao),mynprocg,mynowprog,ilo,iho,  7d21s17
     $           i1so,i1eo,i2so,i2eo)                                   7d20s17
            i10=i1so                                                    7d20s17
            i1n=noc(isa)                                                7d21s17
            i4o=ioooo(is)                                               7d20s17
            nrowo=noc(isbo)*noc(isb)                                    7d21s17
            nnn=idoub(isbo)*idoub(isao)                                 7d21s17
            do i2=i2so,i2eo                                             7d21s17
             if(i2.eq.i2eo)i1n=i1eo                                     7d21s17
             i2m=i2-1
             do i1=i10,i1n                                              7d21s17
              if(i2.le.idoub(isao).and.i1.gt.idoub(isa))then            7d21s17
               i1m=i1-1-idoub(isa)                                      7d21s17
               iadd=iden1(isa)+iacto(isa)*i1m                           7d21s17
               do i4=0,idoub(isbo)-1                                    7d21s17
                do i3=0,iacto(isb)-1                                    7d21s17
                 i3p=i3+idoub(isb)                                       7d21s17
                 i4ou=i4o+i3p+noc(isb)*i4                               7d21s17
                 term=8d0*bc(i4ou)                                      7d21s17
                 iadh=ihessa(isa,isb)+i4+idoub(isbo)                    7d21s17
     $                *(i2m+idoub(isao)*iacto(isa)*i3)                  7d21s17
                 do i5=0,iacto(isa)-1                                   7d21s17
                  jadh=iadh+nnn*i5                                      7d21s17
                  bc(jadh)=bc(jadh)-term*bc(iadd+i5)                    7d21s17
                 end do                                                 7d21s17
                end do                                                  7d21s17
               end do                                                   7d21s17
              end if                                                    7d21s17
              i4o=i4o+nrowo                                             7d21s17
             end do                                                     7d21s17
             i10=1                                                      7d21s17
            end do                                                      7d21s17
           end if                                                       7d21s17
c
c     ddaa part o
c
           if(isblk(1,is).eq.isao.and.isblk(2,is).eq.isb.and.           7d21s17
     $          isblk(4,is).eq.isbo)then                                7d21s17
            igot5=1                                                     7d21s17
            call ilimts(noc(isa),noc(isbo),mynprocg,mynowprog,ilo,iho,  11d3s17
     $           i1so,i1eo,i2so,i2eo)                                   7d20s17
            i10=i1so                                                    7d20s17
            i1n=noc(isa)                                                7d21s17
            i4o=ioooo(is)                                               7d20s17
            if(isb.eq.isao)then                                         7d21s17
             nrowo=(noc(isb)*(noc(isb)+1))/2                            7d21s17
            else                                                        7d21s17
             nrowo=noc(isao)*noc(isb)                                   7d21s17
            end if                                                      7d21s17
            nnn=idoub(isbo)*idoub(isao)                                 7d21s17
            do i2=i2so,i2eo                                             7d21s17
             if(i2.eq.i2eo)i1n=i1eo                                     7d21s17
             i2m=i2-1                                                   7d21s17
             do i1=i10,i1n                                              7d21s17
              if(i1.gt.idoub(isa).and.i2.le.idoub(isbo))then            7d21s17
               i1m=i1-1-idoub(isa)
               iadd=iden1(isa)+iacto(isa)*i1m                           7d21s17
               if(isb.eq.isao)then                                      7d21s17
                do i4=0,iacto(isb)-1                                    7d21s17
                 i4p=i4+idoub(isb)                                      7d21s17
                 do i3=0,idoub(isao)-1                                  7d21s17
                  i4ou=i4o+((i4p*(i4p+1))/2)+i3                         7d21s17
                  term=2d0*bc(i4ou)                                     7d21s17
                  iadh=ihessa(isa,isb)+i2m+idoub(isbo)*(i3+idoub(isao)  7d21s17
     $                 *iacto(isa)*i4)                                  7d21s17
                  do i5=0,iacto(isa)-1                                  7d21s17
                   jadh=iadh+i5*nnn                                     7d21s17
                   bc(jadh)=bc(jadh)+term*bc(iadd+i5)                   7d21s17
                  end do                                                7d21s17
                 end do                                                 7d21s17
                end do                                                  7d21s17
               else
                do i4=0,iacto(isb)-1                                    7d21s17
                 i4p=i4+idoub(isb)                                      7d21s17
                 do i3=0,idoub(isao)-1                                  7d21s17
                  i4ou=i4o+i3+noc(isao)*i4p                             7d21s17
                  term=2d0*bc(i4ou)                                     7d21s17
                  iadh=ihessa(isa,isb)+i2m+idoub(isbo)*(i3+idoub(isao)  7d21s17
     $                 *iacto(isa)*i4)                                  7d21s17
                  do i5=0,iacto(isa)-1                                  7d21s17
                   jadh=iadh+i5*nnn                                     7d21s17
                   bc(jadh)=bc(jadh)+term*bc(iadd+i5)                   7d21s17
                  end do                                                7d21s17
                 end do                                                 7d21s17
                end do                                                  7d21s17
               end if
              end if
              i4o=i4o+nrowo
             end do
             i10=1
            end do
           end if
           if(isblk(2,is).eq.isao.and.isblk(1,is).eq.isb.and.           7d21s17
     $          isblk(4,is).eq.isbo.and.igot5.eq.0)then                 7d21s17
            igot5=1                                                     7d21s17
            call ilimts(noc(isa),noc(isbo),mynprocg,mynowprog,ilo,iho,  11d3s17
     $           i1so,i1eo,i2so,i2eo)                                   7d20s17
            i10=i1so                                                    7d20s17
            i1n=noc(isa)                                                7d21s17
            i4o=ioooo(is)                                               7d20s17
            nrowo=noc(isao)*noc(isb)                                    7d21s17
            nnn=idoub(isbo)*idoub(isao)                                 7d21s17
            do i2=i2so,i2eo                                             7d21s17
             if(i2.eq.i2eo)i1n=i1eo                                     7d21s17
             i2m=i2-1                                                   7d21s17
             do i1=i10,i1n                                              7d21s17
              if(i1.gt.idoub(isa).and.i2.le.idoub(isbo))then            7d21s17
               i1m=i1-1-idoub(isa)
               iadd=iden1(isa)+iacto(isa)*i1m                           7d21s17
               do i4=0,idoub(isao)-1                                     7d21s17
                do i3=0,iacto(isb)-1                                    7d21s17
                 i3p=i3+idoub(isb)                                       7d21s17
                 i4ou=i4o+i3p+noc(isb)*i4                               7d21s17
                 term=2d0*bc(i4ou)                                      7d21s17
                 iadh=ihessa(isa,isb)+i2m+idoub(isbo)*(i4+idoub(isao)   7d21s17
     $                *iacto(isa)*i3)                                   7d21s17
                 do i5=0,iacto(isa)-1                                   7d21s17
                  jadh=iadh+i5*nnn                                      7d21s17
                  bc(jadh)=bc(jadh)+term*bc(iadd+i5)                    7d21s17
                 end do                                                 7d21s17
                end do                                                  7d21s17
               end do                                                   7d21s17
              end if
              i4o=i4o+nrowo
             end do
             i10=1
            end do
           end if
           if(isblk(2,is).eq.isao.and.isblk(1,is).eq.isb.and.           7d21s17
     $          isblk(3,is).eq.isbo.and.igot5.eq.0)then                 7d21s17
            igot5=1                                                     7d21s17
            call ilimts(noc(isbo),noc(isa),mynprocg,mynowprog,ilo,iho,  7d21s17
     $           i1so,i1eo,i2so,i2eo)                                   7d20s17
            i10=i1so                                                    7d20s17
            i1n=noc(isbo)                                                7d21s17
            i4o=ioooo(is)                                               7d20s17
            nrowo=noc(isao)*noc(isb)                                     7d21s17
            nnn=idoub(isbo)*idoub(isao)                                 7d21s17
            do i2=i2so,i2eo                                             7d21s17
             if(i2.eq.i2eo)i1n=i1eo                                     7d21s17
             i2m=i2-1-idoub(isa)                                        7d21s17
             iadd=iden1(isa)+iacto(isa)*i2m                             7d21s17
             do i1=i10,i1n                                              7d21s17
              if(i2.gt.idoub(isa).and.i1.le.idoub(isbo))then            7d21s17
               i1m=i1-1                                                 7d21s17
               do i4=0,idoub(isao)-1                                     7d21s17
                do i3=0,iacto(isb)-1                                    7d21s17
                 i3p=i3+idoub(isb)                                       7d21s17
                 i4ou=i4o+i3p+noc(isb)*i4                               7d21s17
                 term=2d0*bc(i4ou)                                      7d21s17
                 iadh=ihessa(isa,isb)+i1m+idoub(isbo)*(i4+idoub(isao)   7d21s17
     $                *iacto(isa)*i3)                                   7d21s17
                 do i5=0,iacto(isa)-1                                   7d21s17
                  jadh=iadh+i5*nnn                                      7d21s17
                  bc(jadh)=bc(jadh)+term*bc(iadd+i5)                    7d21s17
                 end do                                                 7d21s17
                end do                                                  7d21s17
               end do                                                   7d21s17
              end if
              i4o=i4o+nrowo
             end do
             i10=1
            end do
           end if
           if(isblk(1,is).eq.isao.and.isblk(2,is).eq.isb.and.           7d21s17
     $          isblk(3,is).eq.isbo.and.igot5.eq.0)then                 7d21s17
            igot5=1                                                     7d21s17
            call ilimts(noc(isbo),noc(isa),mynprocg,mynowprog,ilo,iho,  7d21s17
     $           i1so,i1eo,i2so,i2eo)                                   7d20s17
            i10=i1so                                                    7d20s17
            i1n=noc(isbo)                                                7d21s17
            i4o=ioooo(is)                                               7d20s17
            nrowo=noc(isao)*noc(isb)                                     7d21s17
            nnn=idoub(isbo)*idoub(isao)                                 7d21s17
            do i2=i2so,i2eo                                             7d21s17
             if(i2.eq.i2eo)i1n=i1eo                                     7d21s17
             i2m=i2-1-idoub(isa)                                        7d21s17
             iadd=iden1(isa)+iacto(isa)*i2m                             7d21s17
             do i1=i10,i1n                                              7d21s17
              if(i2.gt.idoub(isa).and.i1.le.idoub(isbo))then            7d21s17
               i1m=i1-1                                                 7d21s17
               do i4=0,iacto(isb)-1                                     7d21s17
                i4p=i4+idoub(isb)                                       7d21s17
                do i3=0,idoub(isao)-1                                    7d21s17
                 i4ou=i4o+i3+noc(isao)*i4p                              7d21s17
                 term=2d0*bc(i4ou)                                      7d21s17
                 iadh=ihessa(isa,isb)+i1m+idoub(isbo)*(i3+idoub(isao)   7d21s17
     $                *iacto(isa)*i4)                                   7d21s17
                 do i5=0,iacto(isa)-1                                   7d21s17
                  jadh=iadh+i5*nnn                                      7d21s17
                  bc(jadh)=bc(jadh)+term*bc(iadd+i5)                    7d21s17
                 end do                                                 7d21s17
                end do                                                  7d21s17
               end do                                                   7d21s17
              end if
              i4o=i4o+nrowo
             end do
             i10=1
            end do
           end if
c
c     ddaa part p
c
           if(isblk(1,is).eq.isbo.and.isblk(2,is).eq.isao.and.           7d21s17
     $          isblk(4,is).eq.isb)then                                 7d21s17
            igot6=1                                                     7d21s17
            call ilimts(noc(isa),noc(isb),mynprocg,mynowprog,ilo,iho,   7d21s17
     $           i1so,i1eo,i2so,i2eo)                                   7d20s17
            i10=i1so                                                    7d20s17
            i1n=noc(isa)                                                7d21s17
            i4o=ioooo(is)                                               7d20s17
            if(isbo.eq.isao)then                                         7d21s17
             nrowo=(noc(isbo)*(noc(isbo)+1))/2                          7d21s17
            else                                                        7d21s17
             nrowo=noc(isbo)*noc(isao)                                  7d21s17
            end if                                                      7d21s17
            nnn=idoub(isao)*idoub(isbo)                                 7d21s17
            do i2=i2so,i2eo                                             7d21s17
             if(i2.eq.i2eo)i1n=i1eo                                     7d21s17
             i2m=i2-1-idoub(isb)                                        7d21s17
             do i1=i10,i1n                                              7d21s17
              if(i1.gt.idoub(isa).and.i2.gt.idoub(isb))then             7d21s17
               i1m=i1-1-idoub(isa)                                      7d21s17
               iadd=iden1(isa)+iacto(isa)*i1m                           7d21s17
               if(isbo.eq.isao)then                                     7d21s17
                do i4=0,idoub(isao)-1                                   7d21s17
                 do i3=0,i4-1                                           7d21s17
                  i4ou=i4o+((i4*(i4+1))/2)+i3                           7d21s17
                  term=bc(i4ou)*2d0                                     7d21s17
                  iadh34=ihessa(isa,isb)+i3+idoub(isbo)*(i4+idoub(isao) 7d21s17
     $                 *iacto(isa)*i2m)                                 7d21s17
                  iadh43=ihessa(isa,isb)+i4+idoub(isbo)*(i3+idoub(isao) 7d21s17
     $                 *iacto(isa)*i2m)                                 7d21s17
                  do i5=0,iacto(isa)-1                                  7d21s17
                   jadh34=iadh34+nnn*i5                                 7d21s17
                   jadh43=iadh43+nnn*i5                                 7d21s17
                   part=term*bc(iadd+i5)                                7d21s17
                   bc(jadh34)=bc(jadh34)+part                           7d21s17
                   bc(jadh43)=bc(jadh43)+part                           7d21s17
                  end do                                                7d21s17
                 end do                                                 7d21s17
                 i4ou=i4o+((i4*(i4+1))/2)+i4                            7d21s17
                 term=bc(i4ou)*2d0                                      7d21s17
                 iadh=ihessa(isa,isb)+i4+idoub(isbo)*(i4+idoub(isao)    7d21s17
     $                 *iacto(isa)*i2m)                                 7d21s17
                 do i5=0,iacto(isa)-1                                   7d21s17
                  jadh=iadh+nnn*i5                                      7d21s17
                  bc(jadh)=bc(jadh)+term*bc(iadd+i5)                    7d21s17
                 end do                                                 7d21s17
                end do                                                  7d21s17
               else                                                     7d21s17
                do i4=0,idoub(isao)-1                                   7d21s17
                 do i3=0,idoub(isbo)-1                                  7d21s17
                  i4ou=i4o+i3+noc(isbo)*i4                              7d21s17
                  term=bc(i4ou)*2d0                                     7d21s17
                  iadh=ihessa(isa,isb)+i3+idoub(isbo)*(i4+idoub(isao)   7d21s17
     $                 *iacto(isa)*i2m)                                 7d21s17
                  do i5=0,iacto(isa)-1                                  7d21s17
                   jadh=iadh+nnn*i5                                     7d21s17
                   bc(jadh)=bc(jadh)+term*bc(iadd+i5)                   7d21s17
                  end do                                                7d21s17
                 end do                                                 7d21s17
                end do                                                  7d21s17
               end if                                                   7d21s17
              end if                                                    7d21s17
              i4o=i4o+nrowo                                             7d21s17
             end do                                                     7d21s17
             i10=1                                                      7d21s17
            end do                                                      7d21s17
           else if(isblk(1,is).eq.isbo.and.isblk(2,is).eq.isao.and.     11d3s17
     $          isblk(3,is).eq.isb)then                                 11d3s17
            call ilimts(noc(isb),noc(isa),mynprocg,mynowprog,ilo,iho,   11d3s17
     $           i1so,i1eo,i2so,i2eo)                                   11d3s17
            i10=i1so                                                    11d3s17
            i1n=noc(isb)                                                11d3s17
            i4o=ioooo(is)                                               11d3s17
            nrowo=noc(isbo)*noc(isao)                                   11d3s17
            nnn=idoub(isao)*idoub(isbo)                                 11d3s17
            do i2=i2so,i2eo                                             11d3s17
             if(i2.eq.i2eo)i1n=i1eo                                     11d3s17
             i2m=i2-1-idoub(isa)                                        11d3s17
             do i1=i10,i1n                                              11d3s17
              if(i1.gt.idoub(isb).and.i2.gt.idoub(isa))then             11d3s17
               i1m=i1-1-idoub(isb)                                      11d3s17
               iadd=iden1(isa)+iacto(isa)*i2m                           11d3s17
               do i4=0,idoub(isao)-1                                    11d3s17
                do i3=0,idoub(isbo)-1                                   11d3s17
                 i4ou=i4o+i3+noc(isbo)*i4                               11d3s17
                 term=bc(i4ou)*2d0                                      11d3s17
                 iadh=ihessa(isa,isb)+i3+idoub(isbo)*(i4+idoub(isao)    11d3s17
     $                *iacto(isa)*i1m)                                  11d3s17
                 do i5=0,iacto(isa)-1                                   11d3s17
                  jadh=iadh+nnn*i5                                      11d3s17
                  bc(jadh)=bc(jadh)+term*bc(iadd+i5)                    11d3s17
                 end do                                                 11d3s17
                end do                                                  11d3s17
               end do                                                   11d3s17
              end if                                                    7d21s17
              i4o=i4o+nrowo                                             7d21s17
             end do                                                     7d21s17
             i10=1                                                      7d21s17
            end do                                                      7d21s17
           else if(isblk(2,is).eq.isbo.and.isblk(1,is).eq.isao.and.     11d3s17
     $          isblk(3,is).eq.isb)then                                 11d3s17
            call ilimts(noc(isb),noc(isa),mynprocg,mynowprog,ilo,iho,   11d3s17
     $           i1so,i1eo,i2so,i2eo)                                   11d3s17
            i10=i1so                                                    11d3s17
            i1n=noc(isb)                                                11d3s17
            i4o=ioooo(is)                                               11d3s17
            nrowo=noc(isbo)*noc(isao)                                   11d3s17
            nnn=idoub(isao)*idoub(isbo)                                 11d3s17
            do i2=i2so,i2eo                                             11d3s17
             if(i2.eq.i2eo)i1n=i1eo                                     11d3s17
             i2m=i2-1-idoub(isa)                                        11d3s17
             do i1=i10,i1n                                              11d3s17
              if(i1.gt.idoub(isb).and.i2.gt.idoub(isa))then             11d3s17
               i1m=i1-1-idoub(isb)                                      11d3s17
               iadd=iden1(isa)+iacto(isa)*i2m                           11d3s17
               do i4=0,idoub(isao)-1                                    11d3s17
                do i3=0,idoub(isbo)-1                                   11d3s17
                 i4ou=i4o+i4+noc(isao)*i3                               11d3s17
                 term=bc(i4ou)*2d0                                      11d3s17
                 iadh=ihessa(isa,isb)+i3+idoub(isbo)*(i4+idoub(isao)    11d3s17
     $                *iacto(isa)*i1m)                                  11d3s17
                 do i5=0,iacto(isa)-1                                   11d3s17
                  jadh=iadh+nnn*i5                                      11d3s17
                  bc(jadh)=bc(jadh)+term*bc(iadd+i5)                    11d3s17
                 end do                                                 11d3s17
                end do                                                  11d3s17
               end do                                                   11d3s17
              end if                                                    7d21s17
              i4o=i4o+nrowo                                             7d21s17
             end do                                                     7d21s17
             i10=1                                                      7d21s17
            end do                                                      7d21s17
           else if(isblk(2,is).eq.isbo.and.isblk(1,is).eq.isao.and.     11d3s17
     $          isblk(4,is).eq.isb)then                                 11d3s17
            call ilimts(noc(isa),noc(isb),mynprocg,mynowprog,ilo,iho,   11d3s17
     $           i1so,i1eo,i2so,i2eo)                                   11d3s17
            i10=i1so                                                    11d3s17
            i1n=noc(isa)                                                11d3s17
            i4o=ioooo(is)                                               11d3s17
            nrowo=noc(isbo)*noc(isao)                                   11d3s17
            nnn=idoub(isao)*idoub(isbo)                                 11d3s17
            do i2=i2so,i2eo                                             11d3s17
             if(i2.eq.i2eo)i1n=i1eo                                     11d3s17
             i2m=i2-1-idoub(isb)                                        11d3s17
             do i1=i10,i1n                                              11d3s17
              if(i1.gt.idoub(isa).and.i2.gt.idoub(isb))then             11d3s17
               i1m=i1-1-idoub(isa)                                      11d3s17
               iadd=iden1(isa)+iacto(isa)*i1m                           11d3s17
               do i4=0,idoub(isao)-1                                    11d3s17
                do i3=0,idoub(isbo)-1                                   11d3s17
                 i4ou=i4o+i4+noc(isao)*i3                               11d3s17
                 term=bc(i4ou)*2d0                                      11d3s17
                 iadh=ihessa(isa,isb)+i3+idoub(isbo)*(i4+idoub(isao)    11d3s17
     $                *iacto(isa)*i2m)                                  11d3s17
                 do i5=0,iacto(isa)-1                                   11d3s17
                  jadh=iadh+nnn*i5                                      11d3s17
                  bc(jadh)=bc(jadh)+term*bc(iadd+i5)                    11d3s17
                 end do                                                 11d3s17
                end do                                                  11d3s17
               end do                                                   11d3s17
              end if                                                    7d21s17
              i4o=i4o+nrowo                                             7d21s17
             end do                                                     7d21s17
             i10=1                                                      7d21s17
            end do                                                      7d21s17
           end if                                                       7d21s17
           end if                                                       9d25s17
          end do                                                         7d20s17
   12       format(i3,' integral type: ',4i1)
         end if                                                         7d20s17
c
c     doub-active contributions to two particle density terms
c
c
c     for rotations between active and virtual orbitals
c
         if(idoit(5).ne.0.and.lsym)then                                 11d30s17
         do is=1,nsdlk                                                  7d20s17
          if(isblk(3,is).eq.isbv.and.isblk(4,is).eq.isav)then           4d11s17
           igotj=1                                                      4d11s17
           js1=isblk(1,is)                                              4d11s17
           js2=isblk(2,is)                                              4d11s17
           if(js1.eq.js2)then                                           4d11s17
            nrowj=(noc(js1)*(noc(js1)+1))/2                             4d11s17
           else                                                         4d11s17
            nrowj=noc(js1)*noc(js2)                                     4d11s17
           end if                                                       4d11s17
           igotden=0                                                    4d11s17
           do is2=1,nsdlk                                               4d11s17
            if(isblk(1,is2).eq.isblk(1,is).and.
     $           isblk(2,is2).eq.isblk(2,is).and.isblk(3,is2).eq.isao   4d11s17
     $           .and.isblk(4,is2).eq.isbo)then                         4d11s17
             igotden=1                                                  4d11s17
             if(isbo.eq.isao)then
              nrowd=(iacto(js1)*(iacto(js1)+1))/2                       4d11s17
             else                                                       4d11s17
              nrowd=iacto(js1)*iacto(js2)                               4d11s17
             end if                                                     4d11s17
             i10=i1sh
             i1n=nvirtc(isbv)                                           4d12s17
             ih=ihessc(isa,isb)                                          4d11s17
             ij=jmats(is)                                               4d11s17
             do i2=i2sh,i2eh
              if(i2.eq.i2eh)i1n=i1eh
              do i1=i10,i1n
               if(isbv.eq.isav)then
                do ir=0,iacto(isbo)-1                                   4d11s17
                 irp=ir+idoub(isbo)                                     4d11s17
                 jadh=ih+idoub(isao)+noc(isao)*irp                      4d11s17
                 do iq=0,iacto(isao)-1                                  4d11s17
                  ix=max(ir,iq)                                         4d11s17
                  in=min(ir,iq)                                         4d11s17
                  icol=jdenpt(is2)+nrowd*(((ix*(ix+1))/2)+in)           4d11s17
                  iadh=jadh+iq                                          4d11s17
                  do ip=0,iacto(js2)-1
                   ipp=ip+idoub(js2)
                   do io=0,iacto(js2)-1
                    ix=max(io,ip)
                    in=min(io,ip)
                    iadd=icol+((ix*(ix+1))/2)+in                          4d11s17
                    ix=ix+idoub(js2)                                      4d11s17
                    in=in+idoub(js2)                                      4d11s17
                    iadj=ij+((ix*(ix+1))/2)+in                            4d11s17
                    bc(iadh)=bc(iadh)+bc(iadj)*bc(iadd)*4d0               4d11s17
                   end do                                               4d11s17
                  end do                                                4d11s17
                 end do                                                 4d11s17
                end do                                                  4d11s17
               else
                do ir=0,iacto(isbo)-1                                   4d11s17
                 irp=ir+idoub(isbo)                                     4d11s17
                 jadh=ih+idoub(isao)+noc(isao)*irp                      4d11s17
                 do iq=0,iacto(isao)-1                                  4d11s17
                  icol=jdenpt(is2)+nrowd*(iq+iacto(isao)*ir)            4d11s17
                  iadh=jadh+iq                                          4d11s17
                  do ip=0,iacto(js2)-1                                    4d11s17
                   ipp=ip+idoub(js2)
                   iadj=ij+idoub(js1)+noc(js1)*ipp                      4d12s17
                   do io=0,iacto(js1)-1                                   4d11s17
c
c     get 2 rather than 4 because density is 4 times too large, for
c     we only sum over one of abcd bacd abdc badc, but we are not summing
c     over js1,js2 and js2,js1
c
                    bc(iadh)=bc(iadh)+bc(iadj+io)*bc(icol+io)*2d0       4d12s17
                   end do                                                 4d11s17
                   icol=icol+iacto(js1)                                   4d11s17
                  end do                                                  4d11s17
                 end do                                                 4d11s17
                end do                                                  4d11s17
               end if
               ih=ih+nrowh
               ij=ij+nrowj
              end do
              i10=1
             end do
            end if                                                      4d11s17
            if(isblk(1,is2).eq.isblk(1,is).and.
     $           isblk(2,is2).eq.isblk(2,is).and.isblk(4,is2).eq.isao   4d11s17
     $           .and.isblk(3,is2).eq.isbo.and.igotden.eq.0)then        4d11s17
             igotden=1                                                  4d11s17
             nrowd=iacto(js1)*iacto(js2)                                4d11s17
             i10=i1sh
             i1n=nvirtc(isbv)                                           4d12s17
             ih=ihessc(isa,isb)                                          4d11s17
             ij=jmats(is)                                               4d11s17
             do i2=i2sh,i2eh
              if(i2.eq.i2eh)i1n=i1eh
              do i1=i10,i1n
               do ir=0,iacto(isbo)-1                                    4d11s17
                irp=ir+idoub(isbo)                                      4d11s17
                jadh=ih+idoub(isao)+noc(isao)*irp                       4d11s17
                do iq=0,iacto(isao)-1                                   4d11s17
                 iadh=jadh+iq
                 icol=jdenpt(is2)+nrowd*(ir+iacto(isbo)*iq)             4d11s17
                 do ip=0,iacto(js2)-1                                     4d11s17
                  ipp=ip+idoub(js2)
                  iadj=ij+idoub(js1)+noc(js1)*ipp                           4d11s17
                  do io=0,iacto(js1)-1                                    4d11s17
c
c     get 2 rather than 4 because density is 4 times too large, for
c     we only sum over one of abcd bacd abdc badc, but we are not summing
c     over js1,js2 and js2,js1
c
c
                   bc(iadh)=bc(iadh)+bc(iadj+io)*bc(icol+io)*2d0        4d12s17
                  end do                                                  4d11s17
                  icol=icol+iacto(js1)                                    4d11s17
                 end do                                                 4d11s17
                end do                                                  4d11s17
               end do                                                   4d11s17
               ih=ih+nrowh
               ij=ij+nrowj
              end do
              i10=1
             end do
            end if                                                      4d11s17
           end do                                                       4d11s17
           if(igotden.eq.0)then
            write(6,*)('missed 2prt density for '),js1,js2,isao,isbo
            do is2=1,nsdlk
             write(6,*)(isblk(j,is2),j=1,4)
            end do
            call dws_sync
            call dws_finalize
            stop
           end if
          else if(isblk(4,is).eq.isbv.and.isblk(3,is).eq.isav)then      12d26s17
c
c     here we store in ihessc(isb,isa) instead to match up with         12d26s17
c     J indicies.                                                       12d26s17
c                                                                       12d26s17
           call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),         12d26s17
     $          mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)               12d26s17
           igotj=1                                                      4d11s17
           js1=isblk(1,is)                                              4d11s17
           js2=isblk(2,is)                                              4d11s17
           nrowj=noc(js1)*noc(js2)                                      12d26s17
           nrowhx=noc(isao)*noc(isbo)                                   12d26s17
           igotden=0                                                    4d11s17
           do is2=1,nsdlk                                               4d11s17
            if(isblk(1,is2).eq.isblk(1,is).and.
     $           isblk(2,is2).eq.isblk(2,is).and.isblk(3,is2).eq.isao   4d11s17
     $           .and.isblk(4,is2).eq.isbo)then                         4d11s17
             igotden=1                                                  4d11s17
             if(isbo.eq.isao)then
              nrowd=(iacto(js1)*(iacto(js1)+1))/2                       4d11s17
             else                                                       4d11s17
              nrowd=iacto(js1)*iacto(js2)                               4d11s17
             end if                                                     4d11s17
             i10=i1s                                                    12d26s17
             i1n=nvirtc(isav)                                           12d26s17
             ih=ihessc(isb,isa)                                         12d26s17
             ij=jmats(is)                                               4d11s17
             do i2=i2s,i2e                                              12d26s17
              if(i2.eq.i2e)i1n=i1e                                      12d26s17
              do i1=i10,i1n
               if(isbv.eq.isav)then
                do ir=0,iacto(isbo)-1                                   4d11s17
                 irp=ir+idoub(isbo)                                     4d11s17
                 jadh=ih+irp+noc(isbo)*idoub(isao)                      12d26s17
                 do iq=0,iacto(isao)-1                                  4d11s17
                  ix=max(ir,iq)                                         4d11s17
                  in=min(ir,iq)                                         4d11s17
                  icol=jdenpt(is2)+nrowd*(((ix*(ix+1))/2)+in)           4d11s17
                  iadh=jadh+iq*noc(isbo)                                12d26s17
                  do ip=0,iacto(js2)-1
                   ipp=ip+idoub(js2)
                   do io=0,iacto(js2)-1
                    ix=max(io,ip)
                    in=min(io,ip)
                    iadd=icol+((ix*(ix+1))/2)+in                          4d11s17
                    ix=ix+idoub(js2)                                      4d11s17
                    in=in+idoub(js2)                                      4d11s17
                    iadj=ij+((ix*(ix+1))/2)+in                            4d11s17
                    bc(iadh)=bc(iadh)+bc(iadj)*bc(iadd)*4d0               4d11s17
                   end do                                               4d11s17
                  end do                                                4d11s17
                 end do                                                 4d11s17
                end do                                                  4d11s17
               else
                do ir=0,iacto(isbo)-1                                   4d11s17
                 irp=ir+idoub(isbo)                                     4d11s17
                 jadh=ih+irp+noc(isbo)*idoub(isao)                      12d26s17
                 do iq=0,iacto(isao)-1                                  4d11s17
                  icol=jdenpt(is2)+nrowd*(iq+iacto(isao)*ir)            4d11s17
                  iadh=jadh+iq*noc(isbo)                                12d26s17
                  do ip=0,iacto(js2)-1                                    4d11s17
                   ipp=ip+idoub(js2)
                   iadj=ij+idoub(js1)+noc(js1)*ipp                      4d12s17
                   do io=0,iacto(js1)-1                                   4d11s17
c
c     get 2 rather than 4 because density is 4 times too large, for
c     we only sum over one of abcd bacd abdc badc, but we are not summing
c     over js1,js2 and js2,js1
c
                    bc(iadh)=bc(iadh)+bc(iadj+io)*bc(icol+io)*2d0       4d12s17
                   end do                                                 4d11s17
                   icol=icol+iacto(js1)                                   4d11s17
                  end do                                                  4d11s17
                 end do                                                 4d11s17
                end do                                                  4d11s17
               end if                                                   12d26s17
               ih=ih+nrowhx                                             12d26s17
               ij=ij+nrowj
              end do
              i10=1
             end do
            end if                                                      4d11s17
            if(isblk(1,is2).eq.isblk(1,is).and.
     $           isblk(2,is2).eq.isblk(2,is).and.isblk(4,is2).eq.isao   4d11s17
     $           .and.isblk(3,is2).eq.isbo.and.igotden.eq.0)then        4d11s17
             igotden=1                                                  4d11s17
             nrowd=iacto(js1)*iacto(js2)                                4d11s17
             i10=i1s                                                    12d26s17
             i1n=nvirtc(isav)                                           12d26s17
             ih=ihessc(isb,isa)                                         12d26s17
             ij=jmats(is)                                               4d11s17
             do i2=i2s,i2e                                              12d26s17
              if(i2.eq.i2e)i1n=i1e                                      12d26s17
              do i1=i10,i1n
               do ir=0,iacto(isbo)-1                                    4d11s17
                irp=ir+idoub(isbo)                                      4d11s17
                jadh=ih+irp+noc(isbo)*idoub(isao)                       12d26s17
                do iq=0,iacto(isao)-1                                   4d11s17
                 iadh=jadh+iq*noc(isbo)                                 12d26s17
                 icol=jdenpt(is2)+nrowd*(ir+iacto(isbo)*iq)             4d11s17
                 do ip=0,iacto(js2)-1                                     4d11s17
                  ipp=ip+idoub(js2)
                  iadj=ij+idoub(js1)+noc(js1)*ipp                           4d11s17
                  do io=0,iacto(js1)-1                                    4d11s17
c
c     get 2 rather than 4 because density is 4 times too large, for
c     we only sum over one of abcd bacd abdc badc, but we are not summing
c     over js1,js2 and js2,js1
c
c
                   bc(iadh)=bc(iadh)+bc(iadj+io)*bc(icol+io)*2d0        4d12s17
                  end do                                                  4d11s17
                  icol=icol+iacto(js1)                                    4d11s17
                 end do                                                 4d11s17
                end do                                                  4d11s17
               end do                                                   4d11s17
               ih=ih+nrowhx                                             12d26s17
               ij=ij+nrowj
              end do
              i10=1
             end do
            end if                                                      4d11s17
           end do                                                       4d11s17
           if(igotden.eq.0)then
            write(6,*)('missed 2prt density for '),js1,js2,isao,isbo
            do is2=1,nsdlk
             write(6,*)is2,(isblk(j,is2),j=1,4)
            end do
            call dws_sync
            call dws_finalize
            stop
           end if
          end if                                                        4d11s17
         end do                                                         4d11s17
         if(igotj.eq.0.and.min(nvirt(isbv),nvirt(isav)).gt.0)then       7d11s22
          write(6,*)('missed Jmat for '),isao,isbo,isbv,isav
          do is=1,nsdlk
           write(6,*)is,(isblk(j,is),j=1,4)
          end do
          call dws_sync
          call dws_finalize
          stop
         end if                                                         4d11s17
         end if                                                         9d25s17
         if(lsym)then                                                   12d4s17
          call dws_gsumf(bc(ihessa(isa,isb)),nrowa*ncola)               12d4s17
         end if                                                         12d4s17
         if(idoit(3).ne.0.and.lsym)then                                 11d30s17
         do is=1,nsdlk                                                  4d11s17
          if(isblk(1,is).eq.isao.and.isblk(2,is).eq.isbo.and.
     $         isblk(3,is).eq.isbv.and.isblk(4,is).eq.isav)then         3d29s16
           ii=jmats(is)
           if(isao.eq.isbo)then
            nrowj=(noc(isbo)*(noc(isbo)+1))/2
           else
            nrowj=noc(isbo)*noc(isao)
           end if
           kk=ihessc(isa,isb)
           kk0=kk
           i10=i1sh
           i1n=nvirtc(isbv)                                             3d29s16
           do i2=i2sh,i2eh
            if(i2.eq.i2eh)i1n=i1eh
            do i1=i10,i1n
             if(isa.eq.isb)then
              do i4=0,idoub(isbo)-1                                     4d11s17
               kku=kk                                                   4d11s17
               do i3=0,idoub(isao)-1                                    4d11s17
                ix=max(i3,i4)
                in=min(i3,i4)
                iadj=ii+((ix*(ix+1))/2)+in
                bc(kku)=-4d0*bc(iadj)                                   4d11s17
                kku=kku+1                                               4d11s17
               end do
               kk=kk+noc(isao)                                          4d11s17
              end do
              kk=kk+iacto(isbo)*noc(isao)                               7d19s17
              ii=ii+nrowj
             else
              do i4=0,idoub(isbo)-1                                     4d11s17
               kku=kk                                                   4d11s17
               iiu=ii                                                   4d11s17
               do i3=0,idoub(isao)-1                                    4d11s17
                bc(kku)=-4d0*bc(iiu)                                    4d11s17
                kku=kku+1                                               4d11s17
                iiu=iiu+1                                               4d11s17
               end do
               kk=kk+noc(isao)                                          4d11s17
               ii=ii+noc(isao)                                          4d11s17
              end do
              kk=kk+iacto(isbo)*noc(isao)                               4d11s17
              ii=ii+iacto(isbo)*noc(isao)                               4d11s17
             end if
            end do
            i10=1
           end do
           go to 1
          end if
          if(isblk(1,is).eq.isbo.and.isblk(2,is).eq.isao.and.
     $         isblk(3,is).eq.isbv.and.isblk(4,is).eq.isav)then         3d29s16
           ii=jmats(is)
           kk=ihessc(isa,isb)
           i10=i1sh
           i1n=nvirtc(isbv)                                             3d29s16
           do i2=i2sh,i2eh
            if(i2.eq.i2eh)i1n=i1eh
            do i1=i10,i1n
             do i4=0,idoub(isbo)-1                                      4d11s17
              kku=kk                                                    4d11s17
              do i3=0,idoub(isao)-1                                     4d11s17
               iadj=ii+i4+noc(isbo)*i3
               bc(kku)=-4d0*bc(iadj)
               kku=kku+1                                                4d11s17
              end do
              kk=kk+noc(isao)                                           4d11s17
             end do
             kk=kk+iacto(isbo)*noc(isao)                                4d11s17
             ii=ii+nrowh
            end do
            i10=1
           end do
           go to 1
          end if
          if(isblk(1,is).eq.isbo.and.isblk(2,is).eq.isao.and.
     $         isblk(4,is).eq.isbv.and.isblk(3,is).eq.isav)then         12d26s17
           call ilimts(nvirtc(isav),nvirtc(isbv),mynprocg,mynowprog,    12d26s17
     $         il,ih,i1s,i1e,i2s,i2e)                                   12d26s17
           ii=jmats(is)
           nrow=noc(isao)*noc(isbo)                                     12d26s17
           kk=ihessc(isb,isa)                                           12d26s17
           i10=i1s                                                      12d26s17
           i1n=nvirtc(isav)                                             12d26s17
           do i2=i2s,i2e                                                12d26s17
            if(i2.eq.i2e)i1n=i1e                                        12d26s17
            do i1=i10,i1n
             do i4=0,idoub(isbo)-1                                      4d11s17
              do i3=0,idoub(isao)-1                                     4d11s17
               iadj=ii+i4+noc(isbo)*i3
               kku=kk+i3+noc(isao)*i4                                   12d26s17
               bc(kku)=-4d0*bc(iadj)
              end do
             end do
             kk=kk+nrow                                                 12d26s17
             ii=ii+nrow                                                 12d26s17
            end do
            i10=1
           end do
           go to 1
          end if
          if(isblk(2,is).eq.isbo.and.isblk(1,is).eq.isao.and.           12d26s17
     $         isblk(4,is).eq.isbv.and.isblk(3,is).eq.isav)then         12d26s17
           call ilimts(nvirtc(isav),nvirtc(isbv),mynprocg,mynowprog,    12d26s17
     $         il,ih,i1s,i1e,i2s,i2e)                                   12d26s17
           ii=jmats(is)
           nrow=noc(isao)*noc(isbo)                                     12d26s17
           kk=ihessc(isb,isa)                                           12d26s17
           i10=i1s                                                      12d26s17
           i1n=nvirtc(isav)                                             12d26s17
           do i2=i2s,i2e                                                12d26s17
            if(i2.eq.i2e)i1n=i1e                                        12d26s17
            do i1=i10,i1n
             do i4=0,idoub(isbo)-1                                      4d11s17
              do i3=0,idoub(isao)-1                                     4d11s17
               iadj=ii+i3+noc(isao)*i4                                  12d26s17
               kku=kk+i3+noc(isao)*i4                                   12d26s17
               bc(kku)=-4d0*bc(iadj)
              end do
             end do
             kk=kk+nrow                                                 12d26s17
             ii=ii+nrow                                                 12d26s17
            end do
            i10=1
           end do
           go to 1
          end if
         end do
         if(min(nvirt(isbv),nvirt(isav)).gt.0)then                      7d11s22
          write(6,*)('in buildhess')
          write(6,*)('could not find jmat for '),isao,isbo,isbv,isav       3d29s16
          do is=1,nsdlk
           write(6,*)is,(isblk(j,is),j=1,4)
          end do
          call dws_sync
          call dws_finalize
          stop
         end if                                                         7d11s22
    1    continue
         end if                                                         9d25s17
         if(idoit(5).ne.0.and.lsym)then                                 11d30s17
         do is=1,nsdlkk
          if(isblkk(3,is).eq.isbv.and.isblkk(4,is).eq.isav)then         4d11s17
           js1=isblkk(1,is)                                             4d11s17
           js2=isblkk(2,is)                                             4d11s17
           igotden=0                                                    4d11s17
           nrowk=noc(js1)*noc(js2)                                      4d11s17
           do is2=1,nsdlk                                               4d11s17
            if(isblk(1,is2).eq.isblk(2,is2))then                        4d12s17
             nrowd=(iacto(js1)*(iacto(js1)+1))/2                        10d20s17
            else                                                        4d11s17
             nrowd=iacto(js1)*iacto(isao)                               10d20s17
            end if                                                      4d11s17
            if(isblk(1,is2).eq.js1.and.isblk(2,is2).eq.isao             10d20s17
     $           .and.isblk(3,is2).eq.js2                               10d20s17
     $           .and.isblk(4,is2).eq.isbo)then                         4d11s17
             igotden=1                                                  4d11s17
             i10=i1sh                                                   4d11s17
             i1n=nvirtc(isbv)                                           4d11s17
             ik=kmats(is)                                               4d11s17
             ih=ihessc(isa,isb)                                          4d11s17
             do i2=i2sh,i2eh                                            4d11s17
              if(i2.eq.i2eh)i1n=i1eh                                    4d11s17
              do i1=i10,i1n                                             4d11s17
               if(js1.eq.isao)then                                      10d19s17
                do ir=0,iacto(isbo)-1                                   4d11s17
                 irp=ir+idoub(isbo)                                     4d11s17
                 jadh=ih+idoub(isao)+noc(isao)*irp                      4d11s17
                 do iq=0,iacto(isao)-1                                  4d11s17
                  iadh=jadh+iq                                          4d11s17
                  do io=0,iacto(js1)-1                                    4d11s17
                   iop=io+idoub(js1)                                    4d11s17
                   ix=max(io,iq)                                        4d11s17
                   in=min(io,iq)                                        4d11s17
                   noq=(((ix*(ix+1))/2)+in)
                   jadd=jdenpt(is2)+(((ix*(ix+1))/2)+in)                4d13s17
                   iadk=ik+iop+noc(js1)*idoub(js2)                       10d20s17
                   do ip=0,iacto(js2)-1                                  10d20s17
                    ix=max(ir,ip)                                         4d11s17
                    in=min(ir,ip)                                         4d11s17
                    nrp=(((ix*(ix+1))/2)+in)
                    iadd=jadd+nrowd*(((ix*(ix+1))/2)+in)                4d13s17
                    bc(iadh)=bc(iadh)+8d0*bc(iadk+ip*noc(js1))*bc(iadd) 10d20s17
                   end do                                               4d11s17
                  end do                                                4d11s17
                 end do                                                 4d11s17
                end do                                                  4d11s17
               else                                                     4d11s17
                do ir=0,iacto(isbo)-1                                   4d11s17
                 irp=ir+idoub(isbo)                                     4d11s17
                 jadh=ih+idoub(isao)+noc(isao)*irp                      4d11s17
                 do iq=0,iacto(isao)-1                                  4d11s17
                  iadh=jadh+iq                                          4d11s17
                  do io=0,iacto(js1)-1                                  10d20s17
                   iop=io+idoub(js1)                                    10d20s17
                   jadd=jdenpt(is2)+io+iacto(js1)*iq+nrowd*iacto(js2)*ir10d20s17
                   iadk=ik+iop+idoub(js2)*noc(js1)                      10d20s17
                   do ip=0,iacto(js2)-1                                 10d20s17
                    iadd=jadd+nrowd*ip                                  4d13s17
c
c     get 2 rather than 8 because density is 4 times too large, for
c     we only sum over one of abcd bacd abdc badc, but we are summing
c     over js1,js2 and js2,js1
c
                    bc(iadh)=bc(iadh)+2d0*bc(iadk+ip*noc(js1))*bc(iadd) 10d20s17
                   end do                                               4d11s17
                  end do                                                4d11s17
                 end do                                                 4d11s17
                end do                                                  4d11s17
               end if                                                   10d19s17
               ik=ik+nrowk                                              4d12s17
               ih=ih+nrowh                                              4d12s17
              end do                                                    4d11s17
              i10=1                                                     4d11s17
             end do                                                     4d11s17
            end if                                                      4d11s17
            if(isblk(1,is2).eq.js1.and.isblk(2,is2).eq.isao             10d20s17
     $           .and.isblk(4,is2).eq.js2                               10d20s17
     $           .and.isblk(3,is2).eq.isbo.and.igotden.eq.0)then        4d12s17
             igotden=1                                                  4d11s17
             i10=i1sh                                                   4d11s17
             i1n=nvirtc(isbv)                                           4d11s17
             ik=kmats(is)                                               4d11s17
             ih=ihessc(isa,isb)                                          4d11s17
             npdim=nrowd*iacto(isbo)                                    4d12s17
             do i2=i2sh,i2eh                                            4d11s17
              if(i2.eq.i2eh)i1n=i1eh                                    4d11s17
              do i1=i10,i1n                                             4d11s17
               do ir=0,iacto(isbo)-1                                    4d11s17
                irp=ir+idoub(isbo)                                      4d11s17
                jadh=ih+idoub(isao)+noc(isao)*irp                       4d11s17
                do iq=0,iacto(isao)-1                                   4d11s17
                 iadh=jadh+iq                                           4d11s17
                 do io=0,iacto(js1)-1                                   10d20s17
                  iop=io+idoub(js1)                                     10d20s17
                  iadk=ik+iop+idoub(js2)*noc(js1)                       10d20s17
                  jadd=jdenpt(is2)+io+iacto(js1)*iq+nrowd*ir            10d20s17
                  do ip=0,iacto(js2)-1                                  10d20s17
                   iadd=jadd+npdim*ip                                   4d13s17
c
c     get 2 rather than 8 because density is 4 times too large, for
c     we only sum over one of abcd bacd abdc badc, but we are summing
c     over js1,js2 and js2,js1
c
                   bc(iadh)=bc(iadh)+2d0*bc(iadk+ip*noc(js1))*bc(iadd)  10d20s17
                  end do                                                4d11s17
                 end do                                                 4d11s17
                end do                                                  4d11s17
               end do                                                   4d11s17
               ik=ik+nrowk                                              4d12s17
               ih=ih+nrowh                                              4d12s17
              end do                                                    4d11s17
              i10=1                                                     4d11s17
             end do                                                     4d11s17
            end if                                                      4d11s17
            if(isblk(2,is2).eq.js1.and.isblk(1,is2).eq.isao             10d20s17
     $           .and.isblk(4,is2).eq.js2                               10d20s17
     $           .and.isblk(3,is2).eq.isbo.and.igotden.eq.0)then        4d12s17
             igotden=1                                                  4d11s17
             i10=i1sh                                                   4d11s17
             i1n=nvirtc(isbv)                                           4d11s17
             ik=kmats(is)                                               4d11s17
             ih=ihessc(isa,isb)                                          4d11s17
             npdim=nrowd*iacto(isbo)                                    4d12s17
             do i2=i2sh,i2eh                                            4d11s17
              if(i2.eq.i2eh)i1n=i1eh                                    4d11s17
              do i1=i10,i1n                                             4d11s17
               do ir=0,iacto(isbo)-1                                    4d11s17
                irp=ir+idoub(isbo)                                      4d11s17
                jadh=ih+idoub(isao)+noc(isao)*irp                       4d11s17
                do iq=0,iacto(isao)-1                                   4d11s17
                 iadh=jadh+iq                                           4d11s17
                 do io=0,iacto(js1)-1                                   10d20s17
                  iop=io+idoub(js1)                                     10d20s17
                  iadk=ik+iop+idoub(js2)*noc(js1)                       10d20s17
                  jadd=jdenpt(is2)+iq+iacto(isao)*io+nrowd*ir           4d13s17
                  do ip=0,iacto(js2)-1                                  10d20s17
                   iadd=jadd+npdim*ip                                   4d13s17
c
c     get 2 rather than 8 because density is 4 times too large, for
c     we only sum over one of abcd bacd abdc badc, but we are summing
c     over js1,js2 and js2,js1
c
                   bc(iadh)=bc(iadh)+2d0*bc(iadk+ip*noc(js1))*bc(iadd)  10d20s17
                  end do                                                4d11s17
                 end do                                                 4d11s17
                end do                                                  4d11s17
               end do                                                   4d11s17
               ik=ik+nrowk                                              4d12s17
               ih=ih+nrowh                                              4d12s17
              end do                                                    4d11s17
              i10=1                                                     4d11s17
             end do                                                     4d11s17
            end if                                                      4d11s17
            if(isblk(2,is2).eq.js1.and.isblk(1,is2).eq.isao             10d20s17
     $           .and.isblk(3,is2).eq.js2                               10d20s17
     $           .and.isblk(4,is2).eq.isbo.and.igotden.eq.0)then        4d12s17
             igotden=1                                                  4d11s17
             i10=i1sh                                                   4d11s17
             i1n=nvirtc(isbv)                                           4d11s17
             ik=kmats(is)                                               4d11s17
             ih=ihessc(isa,isb)                                          4d11s17
             npdim=nrowd*iacto(isbo)                                    4d12s17
             do i2=i2sh,i2eh                                            4d11s17
              if(i2.eq.i2eh)i1n=i1eh                                    4d11s17
              do i1=i10,i1n                                             4d11s17
               do ir=0,iacto(isbo)-1                                    4d11s17
                irp=ir+idoub(isbo)                                      4d11s17
                jadh=ih+idoub(isao)+noc(isao)*irp                       4d11s17
                do iq=0,iacto(isao)-1                                   4d11s17
                 iadh=jadh+iq                                           4d11s17
                 do io=0,iacto(js1)-1                                   4d11s17
                  iop=io+idoub(js1)                                     4d11s17
                  iadk=ik+iop+idoub(js2)*noc(js1)                       10d20s17
                  jadd=jdenpt(is2)+iq+iacto(isao)*io+nrowd*iacto(js2)   10d20s17
     $                 *ir                                              10d19s17
                  do ip=0,iacto(js2)-1                                  10d20s17
                   iadd=jadd+nrowd*ip                                   4d13s17
c
c     get 2 rather than 8 because density is 4 times too large, for
c     we only sum over one of abcd bacd abdc badc, but we are summing
c     over js1,js2 and js2,js1
c
                   bc(iadh)=bc(iadh)+2d0*bc(iadk+ip*noc(js1))*bc(iadd)  10d20s17
                  end do                                                4d11s17
                 end do                                                 4d11s17
                end do                                                  4d11s17
               end do                                                   4d11s17
               ik=ik+nrowk                                              4d12s17
               ih=ih+nrowh                                              4d12s17
              end do                                                    4d11s17
              i10=1                                                     4d11s17
             end do                                                     4d11s17
            end if                                                      4d11s17
           end do                                                       4d11s17
           if(igotden.eq.0)then
            write(6,*)('missed density ')
            call dws_sync
            call dws_finalize
            stop
           end if
          end if                                                        4d11s17
         end do                                                         4d12s17
         end if                                                         9d25s17
         if(idoit(3).ne.0.and.lsym)then                                 11d30s17
         do is=1,nsdlkk                                                 4d12s17
          if(isblkk(1,is).eq.isao.and.isblkk(2,is).eq.isbo.and.
     $         isblkk(3,is).eq.isbv.and.isblkk(4,is).eq.isav)then       3d29s16
           kk=ihessc(isa,isb)
           ii=kmats(is)
           i10=i1sh
           i1n=nvirtc(isbv)                                             3d29s16
           do i2=i2sh,i2eh
            if(i2.eq.i2eh)i1n=i1eh
            do i1=i10,i1n
             do i4=0,idoub(isbo)-1                                      4d11s17
              do i3=0,idoub(isao)-1                                     4d11s17
               orig=bc(kk+i3)
               bc(kk+i3)=bc(kk+i3)+16d0*bc(ii+i3)
              end do
              kk=kk+noc(isao)                                           4d11s17
              ii=ii+noc(isao)                                           4d11s17
             end do
             kk=kk+iacto(isbo)*noc(isao)                                4d11s17
             ii=ii+iacto(isbo)*noc(isao)                                4d11s17
            end do
            i10=1
           end do
           go to 2
          end if
          if(isblkk(2,is).eq.isao.and.isblkk(1,is).eq.isbo.and.         12d26s17
     $         isblkk(4,is).eq.isbv.and.isblkk(3,is).eq.isav)then       12d26s17
           call ilimts(nvirtc(isav),nvirtc(isbv),mynprocg,mynowprog,    12d26s17
     $         il,ih,i1s,i1e,i2s,i2e)                                   12d26s17
           kk=ihessc(isb,isa)                                           12d26s17
           ii=kmats(is)
           i10=i1s                                                      12d26s17
           i1n=nvirtc(isav)                                             12d26s17
           do i2=i2s,i2e                                                12d26s17
            if(i2.eq.i2e)i1n=i1e                                        12d26s17
            do i1=i10,i1n
             do i4=0,idoub(isao)-1                                      12d26s17
              do i3=0,idoub(isbo)-1                                     12d26s17
               bc(kk+i3)=bc(kk+i3)+16d0*bc(ii+i3)
              end do
              kk=kk+noc(isbo)                                           12d26s17
              ii=ii+noc(isbo)                                           12d26s17
             end do
             kk=kk+iacto(isao)*noc(isbo)                                12d26s17
             ii=ii+iacto(isao)*noc(isbo)                                12d26s17
            end do
            i10=1
           end do
           go to 2
          end if
         end do
         if(min(nvirt(isbv),nvirt(isav)).gt.0)then                      7d11s22
          write(6,*)('in buildhess')
          write(6,*)('could not find kmat for '),isao,isbo,isbv,isav
          do is=1,nsdlkk
           write(6,*)(isblkk(j,is),j=1,4)
          end do
          call dws_sync
          call dws_finalize
          stop
         end if                                                         7d11s22
    2    continue
         do is=1,nsdlkk
          if(isblkk(1,is).eq.isbo.and.isblkk(2,is).eq.isao.and.
     $         isblkk(3,is).eq.isbv.and.isblkk(4,is).eq.isav)then       3d29s16
           kk=ihessc(isa,isb)
           ii=kmats(is)
           i10=i1sh
           i1n=nvirtc(isbv)                                             3d29s16
           do i2=i2sh,i2eh
            if(i2.eq.i2eh)i1n=i1eh
            do i1=i10,i1n
             do i4=0,idoub(isbo)-1
              do i3=0,idoub(isao)-1
               iadk=ii+i4+noc(isbo)*i3
               bc(kk+i3)=bc(kk+i3)-4d0*bc(iadk)                         4d11s17
              end do
              kk=kk+noc(isao)                                           4d11s17
             end do
             kk=kk+iacto(isbo)*noc(isao)                                4d11s17
             ii=ii+nrowh
            end do
            i10=1
           end do
           go to 3
          end if
          if(isblkk(2,is).eq.isbo.and.isblkk(1,is).eq.isao.and.         12d26s17
     $         isblkk(4,is).eq.isbv.and.isblkk(3,is).eq.isav)then       12d26s17
           kk=ihessc(isb,isa)                                           12d26s17
           ii=kmats(is)
           i10=i1s                                                      12d26s17
           i1n=nvirtc(isav)                                             12d26s17
           nrow=noc(isao)*noc(isbo)                                     12d26s17
           do i2=i2s,i2e                                                12d26s17
            if(i2.eq.i2e)i1n=i1e                                        12d26s17
            do i1=i10,i1n
             do i4=0,idoub(isao)-1                                      12d26s17
              do i3=0,idoub(isbo)-1                                     12d26s17
               iadk=ii+i4+noc(isao)*i3                                  12d26s17
               bc(kk+i3)=bc(kk+i3)-4d0*bc(iadk)                         4d11s17
              end do
              kk=kk+noc(isbo)                                           12d26s17
             end do
             kk=kk+iacto(isao)*noc(isbo)                                12d26s17
             ii=ii+nrow                                                 12d26s17
            end do
            i10=1
           end do
           go to 3
          end if
         end do
         if(min(nvirt(isbv),nvirt(isav)).gt.0)then                      7d11s22
          write(6,*)('in buildhess')
          write(6,*)('could not find kmatT for '),isao,isbo,isbv,isav      3d29s16
          do is=1,nsdlkk
           write(6,*)(isblkk(j,is),j=1,4)
          end do
          call dws_sync
          call dws_finalize
          stop
         end if                                                         7d11s22
    3    continue
         end if                                                         9d25s17
         if(isav.eq.isbv.and.nvirt(isbv).gt.0)then                      7d11s22
          i10=i1sh                                                      4d13s17
          i1n=nvirtc(isbv)                                              4d13s17
          ih=ihessc(isa,isb)                                             4d13s17
          nbas=nbasdws(isbv)*nvirtc(isbv)/nvirt(isbv)                   9d25s17
          if(idoit(2).ne.0)then                                         9d25s17
          do i2=i2sh,i2eh                                               4d13s17
           i2m=i2-1+noc(isbv)                                           4d13s17
           if(i2.eq.i2eh)i1n=i1eh                                       4d13s17
           do i1=i10,i1n                                                4d13s17
            iadh0=ipt(isbv)+i1+noc(isbv)+nbas*i2m                       4d25s17
            h0val=h0mo(iadh0)*2d0
            do ir=0,iacto(isbo)-1                                         4d13s17
             irp=ir+idoub(isbo)                                         4d13s17
             iadh=ih+idoub(isbo)+noc(isbo)*irp                          4d13s17
             iadd=iden1(isbo)+iacto(isbo)*ir                             4d13s17
             do iq=0,iacto(isbo)-1                                        4d13s17
              bc(iadh+iq)=bc(iadh+iq)+h0val*bc(iadd+iq)                 4d13s17
             end do                                                       4d13s17
            end do                                                        4d13s17
            ih=ih+nrowh                                                 4d13s17
           end do                                                       4d13s17
           i10=1                                                        4d13s17
          end do                                                        4d13s17
          end if                                                        9d25s17
          naa=(iacto(isbo)*(iacto(isbo)+1))/2                           4d14s17
          nvv=(nvirtc(isbv)*(nvirtc(isbv)+1))/2                         10d2s17
          iaa=ibcoff                                                    4d14s17
          ivv=iaa+naa                                                   10d2s17
          iaa2=ivv+nvv                                                  10d2s17
          ibcoff=iaa2+naa                                               10d2s17
          call enough('buildhesscas.  5',bc,ibc)
          do i=0,naa+nvv-1                                              10d2s17
           bc(iaa+i)=0d0                                                4d14s17
          end do                                                        4d14s17
          jaa2=iaa2                                                     10d2s17
          do i=0,iacto(isbo)-1                                          10d2s17
           do j=0,i-1                                                   10d2s17
            bc(jaa2+j)=1d0                                              10d2s17
           end do                                                       10d2s17
           jaa2=jaa2+i                                                  10d2s17
           bc(jaa2)=2d0                                                 10d2s17
           jaa2=jaa2+1                                                  10d2s17
          end do                                                        10d2s17
          if(idoit(4).ne.0)then                                         10d2s17
c
c     ddaa vn, vn terms
c     part e
c
           do is=1,nsdlk                                                10d2s17
            if(isblk(3,is).eq.isbv.and.isblk(4,is).eq.isbv)then         10d2s17
             is12=isblk(1,is)                                           10d2s17
             call ilimts(nvirtc(isbv),nvirtc(isbv),mynprocg,mynowprog,  10d2s17
     $            il,ih,i1s,i1e,i2s,i2e)                                10d2s17
             i10=i1s                                                    10d2s17
             i1n=nvirtc(isbv)                                           10d2s17
             nrow=(noc(is12)*(noc(is12)+1))/2                           10d2s17
             ijm=jmats(is)-1                                            10d2s17
             do i2=i2s,i2e                                              10d2s17
              if(i2.eq.i2e)i1n=i1e                                      10d2s17
              i2m=i2-1                                                  10d2s17
              do i1=i10,i1n                                             10d2s17
               i1m=i1-1                                                 10d2s17
               if(i1.ge.i2)then                                         10d2s17
                jvv=ivv+((i1m*(i1m+1))/2)+i2m                           10d2s17
                do i34=1,idoub(is12)                                     10d2s17
                 iadi=ijm+(i34*(i34+1))/2                                10d2s17
                 bc(jvv)=bc(jvv)+4d0*bc(iadi)                           10d2s17
                end do                                                   10d2s17
               end if                                                   10d2s17
               ijm=ijm+nrow                                             10d2s17
              end do                                                    10d2s17
              i10=1                                                     10d2s17
             end do                                                     10d2s17
            end if                                                      10d2s17
           end do                                                       10d2s17
c
c     part f
c
           do is=1,nsdlkk                                               10d2s17
            if(isblkk(3,is).eq.isbv.and.isblkk(4,is).eq.isbv)then       10d2s17
             is12=isblkk(1,is)                                          10d2s17
             call ilimts(nvirtc(isbv),nvirtc(isbv),mynprocg,mynowprog,  10d2s17
     $            il,ih,i1s,i1e,i2s,i2e)                                10d2s17
             i10=i1s                                                    10d2s17
             i1n=nvirtc(isbv)                                           10d2s17
             ikm=kmats(is)                                              10d2s17
             nrow=noc(is12)*noc(is12)                                   10d2s17
             do i2=i2s,i2e                                              10d2s17
              if(i2.eq.i2e)i1n=i1e                                      10d2s17
              i2m=i2-1                                                  10d2s17
              do i1=i10,i1n                                             10d2s17
               if(i1.ge.i2)then                                         10d2s17
                jvv=ivv+((i1*(i1-1))/2)+i2m                             10d2s17
                do i34=0,idoub(is12)-1                                   10d2s17
                 iadi=ikm+i34*(noc(is12)+1)                             10d2s17
                 bc(jvv)=bc(jvv)-2d0*bc(iadi)                           10d2s17
                end do                                                   10d2s17
               end if                                                   10d2s17
               ikm=ikm+nrow                                             10d2s17
              end do                                                    10d2s17
              i10=1                                                     10d2s17
             end do                                                     10d2s17
            end if                                                      10d2s17
           end do
c
c     parts a&c                                                         10d2s17
c
           do is=1,nsdlk                                                10d2s17
            if(isblk(3,is).eq.isao.and.isblk(4,is).eq.isao)then         10d2s17
             is12=isblk(1,is)                                           10d2s17
             call ilimts(noc(isao),noc(isao),mynprocg,mynowprog,il,ih,  10d2s17
     $            i1s,i1e,i2s,i2e)                                      10d2s17
             i10=i1s                                                    10d2s17
             i1n=noc(isao)                                              10d2s17
             nrow=(noc(is12)*(noc(is12)+1))/2                           11d17s17
             io=ioooo(is)-1                                             10d2s17
             do i2=i2s,i2e                                              10d2s17
              if(i2.eq.i2e)i1n=i1e                                      10d2s17
              i2m=i2-1-idoub(isao)                                      10d2s17
              iden=iden1(isao)+iacto(isao)*i2m                          10d2s17
              do i1=i10,i1n                                             10d2s17
               if(i1.gt.idoub(isao).and.i2.gt.idoub(isao))then          10d2s17
                i1m=i1-idoub(isao)-1
                do i34=1,idoub(is12)                                    10d2s17
                 iadi=io+(i34*(i34+1))/2                                10d2s17
                 fact=2d0*bc(iadi)                                      10d2s17
                 do i5=0,iacto(isao)-1                                  10d2s17
                  ix=max(i5,i1m)                                        10d2s17
                  in=min(i5,i1m)                                        10d2s17
                  iaap=((ix*(ix+1))/2)+in                               10d2s17
                  jaa=iaa+iaap                                          10d2s17
                  jaa2=iaa2+iaap                                        10d2s17
                  bc(jaa)=bc(jaa)+bc(iden+i5)*fact*bc(jaa2)             10d2s17
                 end do                                                 10d2s17
                end do                                                  10d2s17
               end if                                                   10d2s17
               io=io+nrow                                               10d2s17
              end do                                                    10d2s17
              i10=1                                                     10d2s17
             end do                                                     10d2s17
            end if                                                      10d2s17
c
c     parts b,d
c
            if(isblk(1,is).eq.isao.and.isblk(3,is).eq.isao)then         10d2s17
             is24=isblk(2,is)                                           10d2s17
             call ilimts(noc(isao),noc(is24),mynprocg,mynowprog,il,ih,  10d2s17
     $            i1s,i1e,i2s,i2e)                                      10d2s17
             i10=i1s                                                    10d2s17
             i1n=noc(isao)                                              10d2s17
             io=ioooo(is)                                               10d2s17
             if(isao.eq.is24)then                                       10d2s17
              nrow=(noc(isao)*(noc(isao)+1))/2                          10d2s17
              iswitch=0                                                 10d2s17
             else                                                       10d2s17
              nrow=noc(isao)*noc(is24)                                  10d2s17
              iswitch=1                                                 10d2s17
             end if                                                     10d2s17
             do i2=i2s,i2e                                              10d2s17
              i2m=i2-1                                                  10d2s17
              if(i2.eq.i2e)i1n=i1e                                      10d2s17
              do i1=i10,i1n                                             10d2s17
               if(i1.gt.idoub(isao).and.i2.le.idoub(is24))then          10d2s17
                i1m=i1-1-idoub(isao)                                    10d2s17
                iden=iden1(isao)+iacto(isao)*i1m                        10d2s17
                do i3=0,iacto(isao)-1                                   10d2s17
                 i3p=i3+idoub(isao)                                     10d2s17
                 ix=max(i3p,i2m)                                        10d2s17
                 in=min(i3p,i2m)                                        10d2s17
                 ieq=((ix*(ix+1))/2)+in                                 10d2s17
                 inot=i3p+noc(isao)*i2m                                 10d2s17
                 iadi=io+(inot-ieq)*iswitch+ieq                         10d2s17
                 do i5=0,iacto(isao)-1                                  10d2s17
                  ix=max(i5,i3)                                         10d2s17
                  in=min(i5,i3)                                         10d2s17
                  iaap=((ix*(ix+1))/2)+in                               10d2s17
                  jaa=iaa+iaap                                          10d2s17
                  jaa2=iaa2+iaap                                        10d2s17
                  bc(jaa)=bc(jaa)-bc(iden+i5)*bc(iadi)*bc(jaa2)         10d2s17
                 end do                                                 10d2s17
                end do                                                  10d2s17
               end if                                                   10d2s17
               io=io+nrow                                               10d2s17
              end do                                                    10d2s17
              i10=1                                                     10d2s17
             end do                                                     10d2s17
            else if(isblk(2,is).eq.isao.and.isblk(3,is).eq.isao)then    10d2s17
             is14=isblk(1,is)                                           10d2s17
             call ilimts(noc(isao),noc(is14),mynprocg,mynowprog,il,ih,  10d2s17
     $            i1s,i1e,i2s,i2e)                                      10d2s17
             i10=i1s                                                    10d2s17
             i1n=noc(isao)                                              10d2s17
             io=ioooo(is)                                               10d2s17
             nrow=noc(isao)*noc(is14)                                   10d2s17
             do i2=i2s,i2e                                              10d2s17
              i2m=i2-1                                                  10d2s17
              if(i2.eq.i2e)i1n=i1e                                      10d2s17
              do i1=i10,i1n                                             10d2s17
               if(i1.gt.idoub(isao).and.i2.le.idoub(is14))then          10d2s17
                i1m=i1-1-idoub(isao)                                    10d2s17
                iden=iden1(isao)+iacto(isao)*i1m                        10d2s17
                do i3=0,iacto(isao)-1                                   10d2s17
                 i3p=i3+idoub(isao)                                     10d2s17
                 iadi=io+i2m+noc(is14)*i3p                              10d2s17
                 do i5=0,iacto(isao)-1                                  10d2s17
                  ix=max(i5,i3)                                         10d2s17
                  in=min(i5,i3)                                         10d2s17
                  iaap=((ix*(ix+1))/2)+in                               10d2s17
                  jaa=iaa+iaap                                          10d2s17
                  jaa2=iaa2+iaap                                        10d2s17
                  bc(jaa)=bc(jaa)-bc(iden+i5)*bc(iadi)*bc(jaa2)         10d2s17
                 end do                                                 10d2s17
                end do                                                  10d2s17
               end if                                                   10d2s17
               io=io+nrow                                               10d2s17
              end do                                                    10d2s17
              i10=1                                                     10d2s17
             end do                                                     10d2s17
            else if(isblk(2,is).eq.isao.and.isblk(4,is).eq.isao)then    10d2s17
             is13=isblk(1,is)                                           10d2s17
             call ilimts(noc(is13),noc(isao),mynprocg,mynowprog,il,ih,  10d2s17
     $            i1s,i1e,i2s,i2e)                                      10d2s17
             i10=i1s                                                    10d2s17
             i1n=noc(is13)                                              10d2s17
             io=ioooo(is)                                               10d2s17
             nrow=noc(isao)*noc(is13)                                   10d2s17
             do i2=i2s,i2e                                              10d2s17
              i2m=i2-1-idoub(isao)                                      10d2s17
              iden=iden1(isao)+iacto(isao)*i2m                          10d2s17
              if(i2.eq.i2e)i1n=i1e                                      10d2s17
              do i1=i10,i1n                                             10d2s17
               if(i2.gt.idoub(isao).and.i1.le.idoub(is13))then          10d2s17
                i1m=i1-1                                                10d2s17
                do i3=0,iacto(isao)-1                                   10d2s17
                 i3p=i3+idoub(isao)                                     10d2s17
                 iadi=io+i1m+noc(is13)*i3p                              10d2s17
                 do i5=0,iacto(isao)-1                                  10d2s17
                  ix=max(i5,i3)                                         10d2s17
                  in=min(i5,i3)                                         10d2s17
                  iaap=((ix*(ix+1))/2)+in                               10d2s17
                  jaa=iaa+iaap                                          10d2s17
                  jaa2=iaa2+iaap                                        10d2s17
                  bc(jaa)=bc(jaa)-bc(iden+i5)*bc(iadi)*bc(jaa2)         10d2s17
                 end do                                                 10d2s17
                end do                                                  10d2s17
               end if                                                   10d2s17
               io=io+nrow                                               10d2s17
              end do                                                    10d2s17
              i10=1                                                     10d2s17
             end do                                                     10d2s17
            else if(isblk(1,is).eq.isao.and.isblk(4,is).eq.isao)then    10d2s17
             is23=isblk(2,is)                                           10d2s17
             call ilimts(noc(is23),noc(isao),mynprocg,mynowprog,il,ih,  10d2s17
     $            i1s,i1e,i2s,i2e)                                      10d2s17
             i10=i1s                                                    10d2s17
             i1n=noc(is23)                                              10d2s17
             io=ioooo(is)                                               10d2s17
             nrow=noc(isao)*noc(is23)                                   10d2s17
             do i2=i2s,i2e                                              10d2s17
              i2m=i2-1-idoub(isao)                                      10d2s17
              iden=iden1(isao)+iacto(isao)*i2m                          10d2s17
              if(i2.eq.i2e)i1n=i1e                                      10d2s17
              do i1=i10,i1n                                             10d2s17
               if(i2.gt.idoub(isao).and.i1.le.idoub(is23))then          10d2s17
                i1m=i1-1                                                10d2s17
                do i3=0,iacto(isao)-1                                   10d2s17
                 i3p=i3+idoub(isao)                                     10d2s17
                 iadi=io+i3p+noc(isao)*i1m                              10d2s17
                 do i5=0,iacto(isao)-1                                  10d2s17
                  ix=max(i5,i3)                                         10d2s17
                  in=min(i5,i3)                                         10d2s17
                  iaap=((ix*(ix+1))/2)+in                               10d2s17
                  jaa=iaa+iaap                                          10d2s17
                  jaa2=iaa2+iaap                                        10d2s17
                  bc(jaa)=bc(jaa)-bc(iden+i5)*bc(iadi)*bc(jaa2)         10d2s17
                 end do                                                 10d2s17
                end do                                                  10d2s17
               end if                                                   10d2s17
               io=io+nrow                                               10d2s17
              end do                                                    10d2s17
              i10=1                                                     10d2s17
             end do                                                     10d2s17
            end if                                                      10d2s17
           end do                                                       10d2s17
          end if                                                        10d2s17
          if(idoit(5).ne.0)then                                         9d25s17
          do is=1,nsdlk                                                 4d14s17
           if(isblk(1,is).eq.isbo)then                                  4d14s17
            js2=isblk(2,is)                                             4d14s17
            js3=isblk(3,is)                                             4d14s17
            js4=isblk(4,is)                                             4d14s17
            call ilimts(noc(js3),noc(js4),mynprocg,mynowprog,il,ih,     4d14s17
     $           i1s,i1e,i2s,i2e)                                       4d14s17
            i10=i1s                                                     4d14s17
            i1n=noc(js3)                                                4d14s17
            if(js3.eq.js4)then                                          4d14s17
             nrowd=(iacto(js2)*(iacto(js2)+1))/2                        4d14s17
             nrowo=(noc(js2)*(noc(js2)+1))/2
            else                                                        4d14s17
             nrowd=iacto(js2)*iacto(isbo)                               4d14s17
             nrowo=noc(js2)*noc(isbo)                                   4d14s17
            end if                                                      4d14s17
            i4o=ioooo(is)                                               4d14s17
            do i2=i2s,i2e                                               4d14s17
             if(i2.eq.i2e)i1n=i1e                                       4d14s17
             ip=i2-idoub(js4)-1                                         4d14s17
             do i1=i10,i1n                                              4d14s17
              io=i1-idoub(js3)-1                                        4d14s17
              if(ip.ge.0.and.io.ge.0)then                               7d19s17
               if(js3.eq.js4)then
                ix=max(io,ip)                                            4d14s17
                in=min(io,ip)                                            4d14s17
                jadd=jdenpt(is)+nrowd*(((ix*(ix+1))/2)+in)               4d14s17
                do n=0,iacto(js2)-1                                      4d14s17
                 do iq=0,iacto(isbo)-1                                   4d14s17
                  ix=max(iq,n)                                           4d14s17
                  in=min(iq,n)                                           4d14s17
                  iadd=jadd+((ix*(ix+1))/2)+in                           4d14s17
                  d1=bc(iadd)                                           10d20s17
                  ix=ix+idoub(isbo)                                      4d14s17
                  in=in+idoub(isbo)                                      4d14s17
                  i4o2=i4o+((ix*(ix+1))/2)+in                            4d14s17
                  eri2=bc(i4o2)                                         10d20s17
                  do ir=0,iacto(isbo)-1                                  4d14s17
                   ix=max(ir,n)                                          4d14s17
                   in=min(ir,n)                                          4d14s17
                   iadd=jadd+((ix*(ix+1))/2)+in                          4d14s17
                   ix=ix+idoub(isbo)                                     4d14s17
                   in=in+idoub(isbo)                                     4d14s17
                   i4o1=i4o+((ix*(ix+1))/2)+in                           4d14s17
                   ix=max(iq,ir)                                         4d14s17
                   in=min(iq,ir)                                         4d14s17
                   ix=((ix*(ix+1))/2)+in                                10d20s17
                   jaa=iaa+ix                                           10d20s17
                   jaa2=iaa2+ix                                         10d20s17
                   bc(jaa)=bc(jaa)+(d1*bc(i4o1)+bc(iadd)*eri2)          10d20s17
     $                  *bc(jaa2)                                       10d20s17
                  end do                                                 4d14s17
                 end do                                                  4d14s17
                end do                                                   4d14s17
               else                                                      4d14s17
                do n=0,iacto(js2)-1                                      4d14s17
                 jadd=jdenpt(is)+n*iacto(isbo)                           4d14s17
     $                +nrowd*(io+iacto(js3)*ip)                          4d14s17
                 j4o=i4o+idoub(isbo)+noc(isbo)*(n+idoub(js2))            4d14s17
                 do iq=0,iacto(isbo)-1                                   4d14s17
                  iadd=jadd+iq                                           4d14s17
c
c     get 1 rather than 2 because density is 4 times too large, but
c     we don't sum over both abcd and bacd
c
                  d1=bc(iadd)*0.5d0                                     10d20s17
                  i4o2=j4o+iq                                            4d14s17
                  eri2=bc(i4o2)*0.5d0                                   10d20s17
                  do ir=0,iacto(isbo)-1                                  4d14s17
                   iadd=jadd+ir                                          4d14s17
                   i4o1=j4o+ir                                           4d14s17
                   ix=max(iq,ir)                                         4d14s17
                   in=min(iq,ir)                                         4d14s17
                   ix=((ix*(ix+1))/2)+in                                10d20s17
                   jaa=iaa+ix                                           10d20s17
                   jaa2=iaa2+ix                                         10d20s17
                   bc(jaa)=bc(jaa)+(d1*bc(i4o1)+bc(iadd)*eri2)          10d20s17
     $                  *bc(jaa2)                                       10d20s17
                  end do                                                 4d14s17
                 end do                                                  4d14s17
                end do                                                   4d14s17
               end if                                                    4d14s17
              end if                                                    7d19s17
              i4o=i4o+nrowo                                             4d14s17
             end do                                                     4d14s17
             i10=1                                                      4d14s17
            end do                                                      4d14s17
           else if(isblk(2,is).eq.isbo)then                                  4d14s17
            js1=isblk(1,is)                                             4d14s17
            js3=isblk(3,is)                                             4d14s17
            js4=isblk(4,is)                                             4d14s17
            call ilimts(noc(js3),noc(js4),mynprocg,mynowprog,il,ih,     4d14s17
     $           i1s,i1e,i2s,i2e)                                       4d14s17
            i10=i1s                                                     4d14s17
            i1n=noc(js3)                                                4d14s17
            nrowd=iacto(js1)*iacto(isbo)                                4d14s17
            nrowo=noc(js1)*noc(isbo)                                    4d14s17
            i4o=ioooo(is)                                               4d14s17
            do i2=i2s,i2e                                               4d14s17
             if(i2.eq.i2e)i1n=i1e                                       4d14s17
             ip=i2-idoub(js4)-1                                         4d14s17
             do i1=i10,i1n                                              4d14s17
              io=i1-idoub(js3)-1                                        4d14s17
              if(io.ge.0.and.ip.ge.0)then                               7d19s17
               do n=0,iacto(js1)-1                                       4d14s17
                jadd=jdenpt(is)+n                                        4d14s17
     $               +nrowd*(io+iacto(js3)*ip)                           4d14s17
                j4o=i4o+n+idoub(js1)                                     4d14s17
                do iq=0,iacto(isbo)-1                                    4d14s17
                 iadd=jadd+iq*iacto(js1)                                 4d14s17
c
c     get 1 rather than 2 because density is 4 times too large
c     but we don't sum over both abcd and bacd
c
                 d1=bc(iadd)*0.5d0                                      10d20s17
                 i4o2=j4o+(iq+idoub(isbo))*noc(js1)                      4d14s17
                 eri2=bc(i4o2)*0.5d0                                    10d20s17
                 do ir=0,iacto(isbo)-1                                   4d14s17
                  iadd=jadd+ir*iacto(js1)                                4d14s17
                  i4o1=j4o+(ir+idoub(isbo))*noc(js1)                     4d14s17
                  ix=max(iq,ir)                                          4d14s17
                  in=min(iq,ir)                                          4d14s17
                  ix=((ix*(ix+1))/2)+in                                 10d20s17
                  jaa=iaa+ix                                            10d20s17
                  jaa2=iaa2+ix                                          10d20s17
                  bc(jaa)=bc(jaa)+(d1*bc(i4o1)+bc(iadd)*eri2)           10d20s17
     $                 *bc(jaa2)                                        10d20s17
                 end do                                                  4d14s17
                end do                                                   4d14s17
               end do                                                    4d14s17
              end if                                                    7d19s17
              i4o=i4o+nrowo                                             4d14s17
             end do                                                     4d14s17
             i10=1                                                      4d14s17
            end do                                                      4d14s17
           end if                                                       4d14s17
          end do                                                        4d14s17
          end if                                                        9d25s17
          call dws_gsumf(bc(iaa),naa+nvv)                               10d2s17
          nbas=nbasdws(isbo)*nvirtc(isbo)/nvirt(isbo)                   9d25s17
          if(idoit(2).ne.0)then                                         9d25s17
          jaa=iaa                                                       4d14s17
          do ir=0,iacto(isbo)-1                                         4d14s17
           idr=iden1(isbo)+iacto(isbo)*ir                               4d14s17
           ihr=ipt(isbo)+1+idoub(isbo)+nbas*(ir+idoub(isbo))            10d2s17
           do iq=0,ir                                                   4d14s17
            idq=iden1(isbo)+iacto(isbo)*iq                              10d2s17
            ihq=ipt(isbo)+1+idoub(isbo)+nbas*(iq+idoub(isbo))           4d25s17
            do m=0,iacto(isbo)-1                                        4d14s17
             bc(jaa)=bc(jaa)+bc(idr+m)*h0mo(ihq+m)+bc(idq+m)*h0mo(ihr+m)10d2s17
            end do                                                      4d14s17
            jaa=jaa+1                                                   4d14s17
           end do                                                       4d14s17
          end do                                                        4d14s17
          end if                                                        9d25s17
          i10=i1sh                                                      4d14s17
          i1n=nvirtc(isbv)                                              4d14s17
          ih=ihessc(isa,isb)                                             4d14s17
          do i2=i2sh,i2eh                                               4d14s17
           if(i2.eq.i2eh)i1n=i1eh                                       4d14s17
           do i1=i10,i1n                                                4d14s17
            ix=max(i1,i2)                                               10d2s17
            in=min(i1,i2)                                               10d2s17
            jvv=ivv+((ix*(ix-1))/2)+in-1                                10d2s17
            do i4=0,iacto(isbo)-1                                       10d2s17
             i4p=i4+idoub(isbo)                                         10d2s17
             iden=iden1(isbo)+iacto(isbo)*i4                            10d2s17
             jh=ih+idoub(isbo)+noc(isbo)*i4p                            10d2s17
             do i3=0,iacto(isbo)-1                                      10d2s17
              bc(jh+i3)=bc(jh+i3)+bc(jvv)*bc(iden+i3)                   10d2s17
             end do                                                     10d2s17
            end do                                                      10d2s17
            if(i1.eq.i2)then                                            4d14s17
             jaa=iaa                                                    4d14s17
             ihqr=ih-1+idoub(isbo)+noc(isbo)*idoub(isbo)                4d14s17
             do ir=0,iacto(isbo)-1                                      4d14s17
              do iq=0,ir-1                                              4d14s17
               ihrq=ih+ir+idoub(isbo)+noc(isbo)*(iq+idoub(isbo))        4d14s17
               bc(ihrq)=bc(ihrq)-bc(jaa)                                4d14s17
               ihqr=ih+iq+idoub(isbo)+noc(isbo)*(ir+idoub(isbo))        4d14s17
               bc(ihqr)=bc(ihqr)-bc(jaa)                                4d14s17
               jaa=jaa+1                                                4d14s17
              end do                                                    4d14s17
              ihrr=ihqr+1                                               4d14s17
              bc(ihrr)=bc(ihrr)-bc(jaa)                                 4d14s17
              jaa=jaa+1                                                 4d14s17
             end do                                                     4d14s17
            end if                                                      4d14s17
            ih=ih+nrowh                                                 4d14s17
           end do                                                       4d14s17
           i10=1                                                        4d14s17
          end do                                                        4d14s17
          ibcoff=iaa                                                    10d2s17
c
c     for rotations between doub and active
c
          idd=ibcoff                                                    4d25s17
          ndd=(idoub(isao)*(idoub(isao)+1))/2                           4d25s17
          iaa=idd+ndd                                                   4d25s17
          naa=(iacto(isav)*(iacto(isav)+1))/2                           4d25s17
          ibcoff=iaa+naa                                                4d25s17
          idd2=ibcoff                                                   9d5s17
          ibcoff=idd2+ndd                                               9d5s17
          iddf=ibcoff                                                   9d13s17
          ibcoff=iddf+ndd                                               9d13s17
          iaaf=ibcoff                                                   9d13s17
          ibcoff=iaaf+naa                                               9d13s17
          call enough('buildhesscas.  6',bc,ibc)
          do i=0,naa+ndd*2-1                                            9d5s17
           bc(idd+i)=0d0                                                4d25s17
          end do                                                        4d25s17
c
c     iddf is 2 for diagonals and 1 for off                             9d13s17
c
          do i=0,idoub(isao)-1                                          9d13s17
           ix=iddf+((i*(i+1))/2)                                        9d13s17
           do j=0,i-1                                                   9d13s17
            bc(ix+j)=1d0                                                9d13s17
           end do                                                       9d13s17
           bc(ix+i)=2d0                                                 9d13s17
          end do                                                        9d13s17
c
c     iaaf is 2 for diagonals and 1 for off                             9d13s17
c
          do i=0,iacto(isav)-1                                          9d13s17
           ix=iaaf+((i*(i+1))/2)                                        9d13s17
           do j=0,i-1                                                   9d13s17
            bc(ix+j)=1d0                                                9d13s17
           end do                                                       9d13s17
           bc(ix+i)=2d0                                                 9d13s17
          end do                                                        9d13s17
c
c     aaaa contribution to iaa
c
          if(idoit(5).ne.0)then                                         9d25s17
          do is=1,nsdlk                                                 9d15s17
           if(isblk(1,is).eq.isav.and.jdenpt(is).gt.0)then              9d15s17
            if(isblk(1,is).eq.isblk(2,is))then                          9d15s17
             if(isblk(1,is).eq.isblk(3,is))then                         9d15s17
              fmult=-2d0                                                9d15s17
             else                                                       9d15s17
c
c     this is my guess                                                  9d15s17
c
              fmult=-2d0                                                9d15s17
             end if                                                     9d15s17
            else                                                        9d15s17
             fmult=-1d0                                                 9d15s17
            end if                                                       9d15s17
            is2=isblk(2,is)
            is3=isblk(3,is)                                             9d15s17
            is4=isblk(4,is)                                             9d15s17
            call ilimts(noc(is3),noc(is4),mynprocg,mynowprog,il,ih,     9d15s17
     $           i1s,i1e,i2s,i2e)                                       9d15s17
            i10=i1s                                                     9d15s17
            i1n=noc(is3)                                                9d15s17
            i4o=ioooo(is)                                               9d15s17
            if(isav.eq.is2)then                                         9d15s17
             nrowi=(noc(isav)*(noc(isav)+1))/2                          9d15s17
             nrowd=(iacto(isav)*(iacto(isav)+1))/2                      9d15s17
            else                                                        9d15s17
             nrowi=noc(isav)*noc(is2)                                   9d15s17
             nrowd=iacto(isav)*iacto(is2)                               9d15s17
            end if                                                      9d15s17
            do i2=i2s,i2e                                               9d15s17
             i2m=i2-1-idoub(is4)                                        9d15s17
             if(i2.eq.i2e)i1n=i1e                                       9d15s17
             do i1=i10,i1n                                              9d15s17
              if(i1.gt.idoub(is3).and.i2.gt.idoub(is4))then             9d15s17
               i1m=i1-1-idoub(is3)                                      9d15s17
               ix=max(i1m,i2m)                                          9d15s17
               in=min(i1m,i2m)                                          9d15s17
               ieq=((ix*(ix+1))/2)+in                                   9d15s17
               inot=i1m+iacto(is3)*i2m                                  9d15s17
               iswitch=iabs(is3-is4)                                    9d15s17
               iswitch=iswitch/max(1,iswitch)                           9d15s17
               iadx=(inot-ieq)*iswitch+ieq                               9d15s17
               iden=jdenpt(is)+nrowd*iadx                                9d15s17
               do i4=0,iacto(is2)-1                                     9d15s17
                i4p=i4+idoub(is2)                                       9d15s17
                do i3=0,iacto(isav)-1                                   9d15s17
                 i3p=i3+idoub(isav)                                     9d15s17
                 ix=max(i3p,i4p)                                        9d15s17
                 in=min(i3p,i4p)                                        9d15s17
                 ieq=((ix*(ix+1))/2)+in                                 9d15s17
                 inot=i3p+noc(isav)*i4p                                 9d15s17
                 iswitch=iabs(isav-is2)                                 9d15s17
                 iswitch=iswitch/max(1,iswitch)                         9d15s17
                 iad=(inot-ieq)*iswitch+ieq+i4o                         9d15s17
                 fact=fmult*bc(iad)                                     9d15s17
                 do i5=0,iacto(isav)-1                                  9d15s17
                  ix=max(i4,i5)                                         9d15s17
                  in=min(i4,i5)                                         9d15s17
                  ieq=((ix*(ix+1))/2)+in                                9d15s17
                  inot=i5+iacto(isav)*i4                                9d15s17
                  iswitch=iabs(isav-is2)                                9d15s17
                  iswitch=iswitch/max(1,iswitch)                        9d15s17
                  iadd=iden+(inot-ieq)*iswitch+ieq                      9d15s17
                  ix=max(i3,i5)                                         9d15s17
                  in=min(i3,i5)                                         9d15s17
                  ix=((ix*(ix+1))/2)+in                                 9d15s17
                  jaa=iaa+ix                                            9d15s17
                  jaaf=iaaf+ix                                          9d15s17
                  bc(jaa)=bc(jaa)+fact*bc(iadd)*bc(jaaf)                9d15s17
                 end do                                                 9d15s17
                end do                                                  9d15s17
               end do                                                   9d15s17
              end if                                                    9d15s17
              i4o=i4o+nrowi                                             9d15s17
             end do                                                     9d15s17
             i10=1                                                      9d15s17
            end do                                                      9d15s17
           else if(isblk(2,is).eq.isav.and.jdenpt(is).gt.0)then              9d15s17
            fmult=-1d0                                                  9d15s17
            is1=isblk(1,is)
            is3=isblk(3,is)                                             9d15s17
            is4=isblk(4,is)                                             9d15s17
            call ilimts(noc(is3),noc(is4),mynprocg,mynowprog,il,ih,     9d15s17
     $           i1s,i1e,i2s,i2e)                                       9d15s17
            i10=i1s                                                     9d15s17
            i1n=noc(is3)                                                9d15s17
            i4o=ioooo(is)                                               9d15s17
            if(isav.eq.is1)then                                         9d15s17
             nrowi=(noc(isav)*(noc(isav)+1))/2                          9d15s17
             nrowd=(iacto(isav)*(iacto(isav)+1))/2                      9d15s17
            else                                                        9d15s17
             nrowi=noc(isav)*noc(is1)                                   9d15s17
             nrowd=iacto(isav)*iacto(is1)                               9d15s17
            end if                                                      9d15s17
            do i2=i2s,i2e                                               9d15s17
             i2m=i2-1-idoub(is4)                                        9d15s17
             if(i2.eq.i2e)i1n=i1e                                       9d15s17
             do i1=i10,i1n                                              9d15s17
              if(i1.gt.idoub(is3).and.i2.gt.idoub(is4))then             9d15s17
               i1m=i1-1-idoub(is3)                                      9d15s17
               ix=max(i1m,i2m)                                          9d15s17
               in=min(i1m,i2m)                                          9d15s17
               ieq=((ix*(ix+1))/2)+in                                   9d15s17
               inot=i1m+iacto(is3)*i2m                                  9d15s17
               iswitch=iabs(is3-is4)                                    9d15s17
               iswitch=iswitch/max(1,iswitch)                           9d15s17
               iadx=(inot-ieq)*iswitch+ieq                               9d15s17
               iden=jdenpt(is)+nrowd*iadx                                9d15s17
               do i4=0,iacto(is1)-1                                     9d15s17
                i4p=i4+idoub(is1)                                       9d15s17
                do i3=0,iacto(isav)-1                                   9d15s17
                 i3p=i3+idoub(isav)                                     9d15s17
                 ix=max(i3p,i4p)                                        9d15s17
                 in=min(i3p,i4p)                                        9d15s17
                 ieq=((ix*(ix+1))/2)+in                                 9d15s17
                 inot=i4p+noc(is1)*i3p                                  9d15s17
                 iswitch=iabs(isav-is1)                                 9d15s17
                 iswitch=iswitch/max(1,iswitch)                         9d15s17
                 iad=(inot-ieq)*iswitch+ieq+i4o                         9d15s17
                 fact=fmult*bc(iad)                                     9d15s17
                 do i5=0,iacto(isav)-1                                  9d15s17
                  ix=max(i4,i5)                                         9d15s17
                  in=min(i4,i5)                                         9d15s17
                  ieq=((ix*(ix+1))/2)+in                                9d15s17
                  inot=i4+iacto(is1)*i5                                 9d15s17
                  iswitch=iabs(isav-is1)                                9d15s17
                  iswitch=iswitch/max(1,iswitch)                        9d15s17
                  iadd=iden+(inot-ieq)*iswitch+ieq                      9d15s17
                  ix=max(i3,i5)                                         9d15s17
                  in=min(i3,i5)                                         9d15s17
                  ix=((ix*(ix+1))/2)+in                                 9d15s17
                  jaa=iaa+ix                                            9d15s17
                  jaaf=iaaf+ix                                          9d15s17
                  bc(jaa)=bc(jaa)+fact*bc(iadd)*bc(jaaf)                9d15s17
                 end do                                                 9d15s17
                end do                                                  9d15s17
               end do                                                   9d15s17
              end if                                                    9d15s17
              i4o=i4o+nrowi                                             9d15s17
             end do                                                     9d15s17
             i10=1                                                      9d15s17
            end do                                                      9d15s17
           end if                                                       9d15s17
          end do                                                        9d15s17
          end if                                                        9d25s17
          do is=1,nsdlk                                                 4d25s17
c
c     ddaa: dd' part for part i and j.
c
           if(idoit(4).ne.0)then                                        9d25s17
           if(isblk(1,is).eq.isao.and.isblk(2,is).eq.isbo.and.          9d5s17
     $        isblk(3,is).eq.isblk(4,is))then                           9d5s17
            is34=isblk(3,is)                                            9d5s17
            call ilimts(noc(is34),noc(is34),mynprocg,mynowprog,il,ih,   9d5s17
     $           i1s,i1e,i2s,i2e)                                       9d5s17
            i10=i1s                                                     9d5s17
            i1n=noc(is34)                                               9d5s17
            i4o=ioooo(is)                                               10d31s17
            nrow=(noc(isao)*(noc(isao)+1))/2                            9d5s17
            do i2=i2s,i2e                                               9d5s17
             if(i2.eq.i2e)i1n=i1e                                       9d5s17
             do i1=i10,i1n                                              9d5s17
              if(i1.le.idoub(is34).and.i1.eq.i2)then                    10d31s17
               do i4=0,idoub(isbo)-1                                    10d31s17
                do i3=0,idoub(isao)-1                                   10d31s17
                 in=min(i3,i4)                                          10d31s17
                 ix=max(i3,i4)                                          10d31s17
                 ix=((ix*(ix+1))/2)+in                                  10d31s17
                 jdd2=idd2+ix                                             9d13s17
                 iad=i4o+ix                                             10d31s17
                 bc(jdd2)=bc(jdd2)+2d0*bc(iad)*bc(iddf+ix)              10d31s17
                end do                                                  10d31s17
               end do                                                   10d31s17
              end if                                                    10d31s17
              i4o=i4o+nrow                                              10d31s17
             end do                                                     10d31s17
             i10=1                                                      10d31s17
            end do                                                      10d31s17
           end if                                                       10d31s17
           if(isblk(1,is).eq.isao.and.isblk(3,is).eq.isbo.and.          9d5s17
     $        isblk(2,is).eq.isblk(4,is))then                           9d5s17
            is24=isblk(2,is)                                            9d5s17
            call ilimts(noc(isbo),noc(is24),mynprocg,mynowprog,il,ih,   9d5s17
     $           i1s,i1e,i2s,i2e)                                       9d5s17
            i10=i1s                                                     9d5s17
            i1n=noc(isbo)                                               9d5s17
            i4o=ioooo(is)-1                                             9d5s17
            if(isao.eq.is24)then                                        9d5s17
             nrow=(noc(isao)*(noc(isao)+1))/2                           9d5s17
            else                                                        9d5s17
             nrow=noc(isao)*noc(is24)                                    9d5s17
            end if                                                      9d5s17
            do i2=i2s,i2e                                               9d5s17
             if(i2.eq.i2e)i1n=i1e                                       9d5s17
             do i1=i10,i1n                                              9d5s17
              if(i1.le.idoub(isbo).and.i2.le.idoub(is24))then           9d5s17
               if(isao.eq.is24)then                                     9d5s17
                do i3=1,idoub(isao)                                     9d5s17
                 ix=max(i1,i3)                                          9d5s17
                 in=min(i1,i3)                                          9d5s17
                 ix=((ix*(ix-1))/2)+in-1                                9d13s17
                 jdd2=idd2+ix                                           9d13s17
                 iaddf=iddf+ix                                          9d13s17
                 ix=max(i3,i2)                                          9d5s17
                 in=min(i3,i2)                                          9d5s17
                 iad=i4o+((ix*(ix-1))/2)+in                             9d5s17
                 bc(jdd2)=bc(jdd2)-bc(iaddf)*bc(iad)                    9d13s17
                end do                                                   9d5s17
               else                                                     9d5s17
                do i3=1,idoub(isao)                                     9d5s17
                 ix=max(i1,i3)                                          9d5s17
                 in=min(i1,i3)                                          9d5s17
                 ix=((ix*(ix-1))/2)+in-1                                9d13s17
                 jdd2=idd2+ix                                           9d13s17
                 iaddf=iddf+ix                                          9d13s17
                 iad=i4o+i3+noc(isao)*(i2-1)                            1d11s17
                 bc(jdd2)=bc(jdd2)-bc(iaddf)*bc(iad)                    9d13s17
                end do                                                   9d5s17
               end if                                                   9d5s17
              end if                                                    9d5s17
              i4o=i4o+nrow                                              9d5s17
             end do                                                     9d5s17
             i10=1                                                      9d5s17
            end do                                                      9d5s17
           else if(isblk(2,is).eq.isao.and.isblk(4,is).eq.isbo.and.     10d31s17
     $        isblk(1,is).eq.isblk(3,is))then                           10d31s17
            is13=isblk(1,is)                                            10d31s17
            call ilimts(noc(is13),noc(isbo),mynprocg,mynowprog,il,ih,   10d31s17
     $           i1s,i1e,i2s,i2e)                                       9d5s17
            i10=i1s                                                     9d5s17
            i1n=noc(is13)                                               10d31s17
            i4o=ioooo(is)-1                                             9d5s17
            nrow=noc(isao)*noc(is13)                                    10d31s17
            do i2=i2s,i2e                                               9d5s17
             if(i2.eq.i2e)i1n=i1e                                       9d5s17
             do i1=i10,i1n                                              9d5s17
              if(i1.le.idoub(is13).and.i2.le.idoub(isbo))then           10d31s17
                do i3=1,idoub(isao)                                     9d5s17
                 ix=max(i2,i3)                                          9d5s17
                 in=min(i2,i3)                                          9d5s17
                 ix=((ix*(ix-1))/2)+in-1                                9d13s17
                 jdd2=idd2+ix                                           9d13s17
                 iaddf=iddf+ix                                          9d13s17
                 iad=i4o+i1+noc(is13)*(i3-1)                            1d11s17
                 bc(jdd2)=bc(jdd2)-bc(iaddf)*bc(iad)                    9d13s17
                end do                                                   9d5s17
              end if                                                    9d5s17
              i4o=i4o+nrow                                              9d5s17
             end do                                                     9d5s17
             i10=1                                                      9d5s17
            end do                                                      9d5s17
           end if                                                       9d5s17
           end if                                                       9d25s17
c
c     dddd delta nn part
c
           if(isblk(1,is).eq.isao.and.isblk(2,is).eq.isao)then          4d25s17
            is34=isblk(3,is)                                            4d25s17
            call ilimts(noc(is34),noc(is34),mynprocg,mynowprog,il,ih,   4d25s17
     $           i1s,i1e,i2s,i2e)                                       4d25s17
            i10=i1s                                                     4d25s17
            i1n=noc(is34)                                               4d25s17
            i4o=ioooo(is)                                               4d25s17
            nrow=(noc(isao)*(noc(isao)+1))/2                            4d25s17
            do i2=i2s,i2e                                               4d25s17
             i2m=i2-1-idoub(is34)                                       9d6s17
             if(i2.eq.i2e)i1n=i1e                                       4d25s17
             do i1=i10,i1n                                              4d25s17
c
c     J part - note i34 is loop over i3 and i4, but since both idd and  4d25s17
c     i4o are stored as triangles, we can do them together.             4d25s17
c                                                                       4d25s17
              if(i1.eq.i2.and.i1.le.idoub(is34).and.idoit(3).ne.0)then  9d25s17
c     dddd part
               do i34=0,ndd-1                                           4d25s17
                bc(idd+i34)=bc(idd+i34)-8d0*bc(i4o+i34)                 4d25s17
               end do                                                   4d25s17
              end if                                                    9d6s17
              if(i1.gt.idoub(is34).and.i2.gt.idoub(is34).               9d25s17
     $             and.idoit(4).ne.0)then                               9d25s17
               i1m=i1-1-idoub(is34)                                     9d6s17
c     ddaa, part a
               iadd=iden1(is34)+i1m+iacto(is34)*i2m                     9d6s17
               fact=-4d0*bc(iadd)                                       9d6s17
               do i34=0,ndd-1                                           9d6s17
                bc(idd+i34)=bc(idd+i34)+fact*bc(i4o+i34)                9d6s17
               end do                                                   9d6s17
              end if                                                    4d25s17
              i4o=i4o+nrow                                              4d25s17
             end do                                                     4d25s17
             i10=1                                                      4d25s17
            end do                                                      4d25s17
           end if                                                       4d25s17
c
c     ddaa, part b
c
           if(idoit(4).ne.0)then                                        9d25s17
           igotddaab=0                                                   9d6s17
           if(isblk(1,is).eq.isao.and.isblk(3,is).eq.isao)then          9d6s17
            igotddaab=1                                                 9d6s17
            is24=isblk(2,is)                                            9d6s17
            call ilimts(noc(isao),noc(is24),mynprocg,mynowprog,il,ih,   4d25s17
     $           i1s,i1e,i2s,i2e)                                       4d25s17
            i10=i1s                                                     4d25s17
            i1n=noc(isao)                                               4d25s17
            i4o=ioooo(is)                                               4d25s17
            if(isao.eq.is24)then                                        9d6s17
             nrow=(noc(isao)*(noc(isao)+1))/2                            4d25s17
            else                                                        9d6s17
             nrow=noc(isao)*noc(is24)                                   9d6s17
            end if                                                      9d6s17
            do i2=i2s,i2e                                               9d6s17
             if(i2.eq.i2e)i1n=i1e                                       9d6s17
             i2m=i2-1-idoub(is24)                                       9d6s17
             do i1=i10,i1n                                              9d6s17
              if(i1.le.idoub(isao).and.i2.gt.idoub(is24))then           9d6s17
               i1m=i1-1                                                 9d6s17
               jadd=iden1(is24)+iacto(is24)*i2m                         9d6s17
               if(isao.eq.is24)then                                     9d6s17
                do i4=0,idoub(isao)-1                                    9d6s17
                 ix=max(i1m,i4)                                         9d6s17
                 in=min(i1m,i4)                                         9d6s17
                 ix=((ix*(ix+1))/2)+in                                  9d13s17
                 jdd=idd+ix                                             9d13s17
                 jddf=iddf+ix                                           9d13s17
                 do i3=0,iacto(is24)-1                                   9d6s17
                  i3p=i3+idoub(is24)                                    9d6s17
                  ix=max(i4,i3p)                                        9d6s17
                  in=min(i4,i3p)
                  iade=i4o+((ix*(ix+1))/2)+in                           9d6s17
                  bc(jdd)=bc(jdd)+bc(iade)*bc(jadd+i3)*bc(jddf)         9d13s17
                 end do                                                  9d6s17
                end do
               else                                                     9d6s17
                do i3=0,iacto(is24)-1                                   9d6s17
                 i3p=i3+idoub(is24)                                     9d6s17
                 fact=bc(jadd+i3)                                       9d8s17
                 jade=i4o+noc(isao)*i3p                                 9d6s17
                 do i4=0,idoub(isao)-1                                    9d6s17
                  ix=max(i1m,i4)                                        9d6s17
                  in=min(i1m,i4)                                        9d6s17
                  ix=((ix*(ix+1))/2)+in                                 9d13s17
                  jdd=idd+ix                                            9d13s17
                  jddf=iddf+ix                                          9d13s17
                  bc(jdd)=bc(jdd)+bc(jddf)*fact*bc(jade+i4)             9d13s17
                 end do                                                  9d6s17
                end do
               end if                                                   9d6s17
              end if                                                    9d6s17
              i4o=i4o+nrow                                              9d6s17
             end do                                                     9d6s17
             i10=1                                                      9d6s17
            end do                                                      9d6s17
           end if                                                       9d6s17
           if(isblk(2,is).eq.isao.and.isblk(3,is).eq.isao.and.          9d6s17
     $          igotddaab.eq.0)then                                     9d6s17
            igotddaab=1                                                 9d6s17
            is14=isblk(1,is)                                            9d6s17
            call ilimts(noc(isao),noc(is14),mynprocg,mynowprog,il,ih,   4d25s17
     $           i1s,i1e,i2s,i2e)                                       4d25s17
            i10=i1s                                                     4d25s17
            i1n=noc(isao)                                               4d25s17
            i4o=ioooo(is)                                               4d25s17
            nrow=noc(isao)*noc(is14)                                    9d6s17
            do i2=i2s,i2e                                               9d6s17
             if(i2.eq.i2e)i1n=i1e                                       9d6s17
             i2m=i2-1-idoub(is14)                                       9d6s17
             do i1=i10,i1n                                              9d6s17
              if(i1.le.idoub(isao).and.i2.gt.idoub(is14))then           9d6s17
               i1m=i1-1                                                 9d6s17
               jadd=iden1(is14)+iacto(is14)*i2m                         9d6s17
               do i4=0,idoub(isao)-1                                    9d6s17
                ix=max(i1m,i4)                                          9d6s17
                in=min(i1m,i4)                                          9d6s17
                ix=((ix*(ix+1))/2)+in                                   9d13s17
                jdd=idd+ix                                              9d13s17
                jddf=iddf+ix                                            9d13s17
                jade=i4o+idoub(is14)+noc(is14)*i4                       9d6s17
                do i3=0,iacto(is14)-1                                   9d6s17
                 bc(jdd)=bc(jdd)+bc(jddf)*bc(jadd+i3)*bc(jade+i3)       9d13s17
                end do                                                  9d6s17
               end do                                                   9d6s17
              end if                                                    9d6s17
              i4o=i4o+nrow                                              9d6s17
             end do                                                     9d6s17
             i10=1                                                      9d6s17
            end do                                                      9d6s17
           end if                                                       9d6s17
           if(isblk(2,is).eq.isao.and.isblk(4,is).eq.isao.and.          9d6s17
     $          igotddaab.eq.0)then                                     9d6s17
            igotddaab=1                                                 9d6s17
            is13=isblk(1,is)                                            9d6s17
            call ilimts(noc(is13),noc(isao),mynprocg,mynowprog,il,ih,   4d25s17
     $           i1s,i1e,i2s,i2e)                                       4d25s17
            i10=i1s                                                     4d25s17
            i1n=noc(is13)                                               4d25s17
            i4o=ioooo(is)                                               4d25s17
            nrow=noc(isao)*noc(is13)                                    9d6s17
            do i2=i2s,i2e                                               9d6s17
             if(i2.eq.i2e)i1n=i1e                                       9d6s17
             i2m=i2-1                                                   9d6s17
             do i1=i10,i1n                                              9d6s17
              if(i2.le.idoub(isao).and.i1.gt.idoub(is13))then           9d6s17
               i1m=i1-1-idoub(is13)                                     9d6s17
               jadd=iden1(is13)+iacto(is13)*i1m                         9d6s17
               do i4=0,idoub(isao)-1                                    9d6s17
                ix=max(i2m,i4)                                          9d6s17
                in=min(i2m,i4)                                          9d6s17
                ix=((ix*(ix+1))/2)+in                                   9d13s17
                jdd=idd+ix                                              9d13s17
                jddf=iddf+ix                                            9d13s17
                jade=i4o+idoub(is13)+noc(is13)*i4                       9d6s17
                do i3=0,iacto(is13)-1                                   9d6s17
                 bc(jdd)=bc(jdd)+bc(jadd+i3)*bc(jade+i3)*bc(jddf)       9d13s17
                end do                                                  9d6s17
               end do                                                   9d6s17
              end if                                                    9d6s17
              i4o=i4o+nrow                                              9d6s17
             end do                                                     9d6s17
             i10=1                                                      9d6s17
            end do                                                      9d6s17
           end if                                                       9d6s17
           if(isblk(1,is).eq.isao.and.isblk(4,is).eq.isao.and.          9d6s17
     $          igotddaab.eq.0)then                                     9d6s17
            igotddaab=1                                                 9d6s17
            is23=isblk(2,is)                                            9d6s17
            call ilimts(noc(is23),noc(isao),mynprocg,mynowprog,il,ih,   4d25s17
     $           i1s,i1e,i2s,i2e)                                       4d25s17
            i10=i1s                                                     4d25s17
            i1n=noc(is23)                                               4d25s17
            i4o=ioooo(is)                                               4d25s17
            nrow=noc(isao)*noc(is23)                                    9d6s17
            do i2=i2s,i2e                                               9d6s17
             if(i2.eq.i2e)i1n=i1e                                       9d6s17
             i2m=i2-1                                                   9d6s17
             do i1=i10,i1n                                              9d6s17
              if(i2.le.idoub(isao).and.i1.gt.idoub(is23))then           9d6s17
               i1m=i1-1-idoub(is23)                                     9d6s17
               do i3=0,iacto(is13)-1                                    9d6s17
                jade=i4o+noc(isao)*(i3+idoub(is13))                     9d6s17
                jadd=iden1(is23)+i3+iacto(is23)*i1m                     9d6s17
                fact=bc(jadd)                                           9d8s17
                do i4=0,idoub(isao)-1                                    9d6s17
                 ix=max(i2m,i4)                                          9d6s17
                 in=min(i2m,i4)                                          9d6s17
                 ix=((ix*(ix+1))/2)+in                                  9d13s17
                 jdd=idd+ix                                             9d13s17
                 jddf=iddf+ix                                           9d13s17
                 bc(jdd)=bc(jdd)+bc(jddf)*fact*bc(jade+i4)              9d13s17
                end do                                                  9d6s17
               end do                                                   9d6s17
              end if                                                    9d6s17
              i4o=i4o+nrow                                              9d6s17
             end do                                                     9d6s17
             i10=1                                                      9d6s17
            end do                                                      9d6s17
           end if                                                       9d6s17
           end if                                                       9d25s17
c
c     dddd delta dd J part
c
           if(idoit(3).ne.0)then                                        9d25s17
           if(isblk(1,is).eq.isav.and.isblk(2,is).eq.isav)then          4d25s17
            is34=isblk(3,is)                                            4d25s17
            call ilimts(noc(is34),noc(is34),mynprocg,mynowprog,il,ih,   4d25s17
     $           i1s,i1e,i2s,i2e)                                       4d25s17
            i10=i1s                                                     4d25s17
            i1n=noc(is34)                                               4d25s17
            i4o=ioooo(is)                                               4d25s17
            nrow=(noc(isav)*(noc(isav)+1))/2                            4d25s17
            do i2=i2s,i2e                                               4d25s17
             if(i2.eq.i2e)i1n=i1e                                       4d25s17
             do i1=i10,i1n                                              4d25s17
              if(i1.eq.i2.and.i1.le.idoub(is34))then                    4d25s17
               do i4=0,iacto(isav)-1                                    5d1s17
                i4p=i4+idoub(isav)                                      5d1s17
                do i3=0,i4                                              4d25s17
                 i3p=i3+idoub(isav)                                     5d1s17
                 iada=iaa+((i4*(i4+1))/2)+i3                            4d25s17
                 iado=i4o+((i4p*(i4p+1))/2)+i3p                         4d25s17
                 bc(iada)=bc(iada)+8d0*bc(iado)                         4d25s17
                end do                                                  4d25s17
               end do                                                   4d25s17
              end if                                                    4d25s17
              i4o=i4o+nrow                                              4d25s17
             end do                                                     4d25s17
             i10=1                                                      4d25s17
            end do                                                      4d25s17
           end if                                                       4d25s17
c
c     dddd delta nn, K part
c
           igot=0                                                       4d25s17
           if(isblk(1,is).eq.isao.and.isblk(3,is).eq.isao)then          4d25s17
            igot=1
            is24=isblk(2,is)                                            4d25s17
            call ilimts(noc(isao),noc(is24),mynprocg,mynowprog,il,ih,   4d25s17
     $           i1s,i1e,i2s,i2e)                                       4d25s17
            i10=i1s                                                     4d25s17
            i1n=noc(isao)                                               4d25s17
            i4o=ioooo(is)                                               4d25s17
            if(isao.eq.is24)then                                        4d25s17
             nrow=(noc(isao)*(noc(isao)+1))/2                           4d25s17
            else                                                        4d25s17
             nrow=noc(isao)*noc(is24)                                   4d25s17
            end if                                                      4d25s17
            do i2=i2s,i2e                                               4d25s17
             if(i2.eq.i2e)i1n=i1e                                       4d25s17
             i2m=i2-1                                                   4d25s17
             do i1=i10,i1n                                              4d25s17
              i1m=i1-1                                                  4d25s17
              if(i1.le.idoub(isao).and.i2.le.idoub(is24))then           4d25s17
               if(isao.eq.is24)then                                     4d25s17
                do i4=0,idoub(isao)-1                                    4d25s17
                 ix=max(i4,i1m)                                         4d25s17
                 in=min(i4,i1m)                                         4d25s17
                 ix=((ix*(ix+1))/2)+in                                  9d13s17
                 iadd=idd+ix
                 iaddf=iddf+ix                                          9d13s17
                 ix=max(i4,i2m)                                         4d25s17
                 in=min(i4,i2m)                                         4d25s17
                 iado=i4o+((ix*(ix+1))/2)+in                            4d25s17
                 bc(iadd)=bc(iadd)+2d0*bc(iado)*bc(iaddf)               9d13s17
                end do                                                   4d25s17
               else                                                     4d25s17
                do i4=0,idoub(isao)-1                                   4d25s17
                 ix=max(i4,i1m)                                         4d25s17
                 in=min(i4,i1m)                                         4d25s17
                 ix=((ix*(ix+1))/2)+in                                  9d13s17
                 iadd=idd+ix                                            9d13s17
                 iaddf=iddf+ix                                          9d13s17
                 iado=i4o+i4+noc(isao)*i2m                              4d25s17
c     my best guess
                 bc(iadd)=bc(iadd)+2d0*bc(iado)*bc(iaddf)               10d30s17
                end do                                                  4d25s17
               end if                                                   4d25s17
              end if                                                    4d25s17
              i4o=i4o+nrow                                              4d25s17
             end do                                                     4d25s17
             i10=1                                                      4d25s17
            end do                                                      4d25s17
           end if                                                       4d25s17
           if(isblk(1,is).eq.isao.and.isblk(4,is).eq.isao.and.          4d25s17
     $          igot.eq.0)then                                          4d25s17
            igot=1
            is23=isblk(2,is)                                            4d25s17
            call ilimts(noc(is23),noc(isao),mynprocg,mynowprog,il,ih,   4d25s17
     $           i1s,i1e,i2s,i2e)                                       4d25s17
            i10=i1s                                                     4d25s17
            i1n=noc(is23)                                               4d25s17
            i4o=ioooo(is)                                               4d25s17
            nrow=noc(isao)*noc(is23)                                    4d25s17
            do i2=i2s,i2e                                               4d25s17
             if(i2.eq.i2e)i1n=i1e                                       4d25s17
             i2m=i2-1                                                   4d25s17
             do i1=i10,i1n                                              4d25s17
              i1m=i1-1                                                  4d25s17
              if(i1.le.idoub(is23).and.i2.le.idoub(isao))then           4d25s17
               do i4=0,idoub(isao)-1                                    4d25s17
                ix=max(i4,i2m)                                          4d25s17
                in=min(i4,i2m)                                          4d25s17
                ix=((ix*(ix+1))/2)+in                                   9d13s17
                iadd=idd+ix                                             9d13s17
                iaddf=iddf+ix                                           9d13s17
                iado=i4o+i4+noc(isao)*i1m                               4d25s17
c     my best guess
                bc(iadd)=bc(iadd)+2d0*bc(iado)*bc(iaddf)                10d30s17
               end do                                                   4d25s17
              end if                                                    4d25s17
              i4o=i4o+nrow                                              4d25s17
             end do                                                     4d25s17
             i10=1                                                      4d25s17
            end do                                                      4d25s17
           end if                                                       4d25s17
           if(isblk(2,is).eq.isao.and.isblk(4,is).eq.isao.and.          4d25s17
     $          igot.eq.0)then                                          4d25s17
            igot=1
            is13=isblk(1,is)                                            4d25s17
            call ilimts(noc(is13),noc(isao),mynprocg,mynowprog,il,ih,   4d25s17
     $           i1s,i1e,i2s,i2e)                                       4d25s17
            i10=i1s                                                     4d25s17
            i1n=noc(is13)                                               4d25s17
            i4o=ioooo(is)                                               4d25s17
            nrow=noc(isao)*noc(is13)                                    4d25s17
            do i2=i2s,i2e                                               4d25s17
             if(i2.eq.i2e)i1n=i1e                                       4d25s17
             i2m=i2-1                                                   4d25s17
             do i1=i10,i1n                                              4d25s17
              i1m=i1-1                                                  4d25s17
              if(i1.le.idoub(is13).and.i2.le.idoub(isao))then           4d25s17
               do i4=0,idoub(isao)-1                                    4d25s17
                ix=max(i4,i2m)                                          4d25s17
                in=min(i4,i2m)                                          4d25s17
                ix=((ix*(ix+1))/2)+in                                   9d13s17
                iadd=idd+ix                                             9d13s17
                iaddf=iddf+ix                                           9d13s17
                iado=i4o+i1m+noc(is13)*i4                               4d25s17
                bc(iadd)=bc(iadd)+2d0*bc(iado)*bc(iaddf)                10d30s17
               end do                                                   4d25s17
              end if                                                    4d25s17
              i4o=i4o+nrow                                              4d25s17
             end do                                                     4d25s17
             i10=1                                                      4d25s17
            end do                                                      4d25s17
           end if                                                       4d25s17
           if(isblk(2,is).eq.isao.and.isblk(3,is).eq.isao.and.          4d25s17
     $          igot.eq.0)then                                          4d25s17
            igot=1
            is14=isblk(1,is)                                            4d25s17
            call ilimts(noc(isao),noc(is14),mynprocg,mynowprog,il,ih,   4d25s17
     $           i1s,i1e,i2s,i2e)                                       4d25s17
            i10=i1s                                                     4d25s17
            i1n=noc(isao)                                               4d25s17
            i4o=ioooo(is)                                               4d25s17
            nrow=noc(isao)*noc(is14)                                    4d25s17
            do i2=i2s,i2e                                               4d25s17
             if(i2.eq.i2e)i1n=i1e                                       4d25s17
             i2m=i2-1                                                   4d25s17
             do i1=i10,i1n                                              4d25s17
              i1m=i1-1                                                  4d25s17
              if(i1.le.idoub(isao).and.i2.le.idoub(is14))then           4d25s17
               do i4=0,idoub(isao)-1                                    4d25s17
                ix=max(i4,i1m)                                          4d25s17
                in=min(i4,i1m)                                          4d25s17
                ix=((ix*(ix+1))/2)+in                                   9d13s17
                iadd=idd+ix                                             9d13s17
                iado=i4o+i2m+noc(is14)*i4                               4d25s17
c     my best guess
                bc(iadd)=bc(iadd)+2d0*bc(iado)*bc(iaddf)                10d30s17
               end do                                                   4d25s17
              end if                                                    4d25s17
              i4o=i4o+nrow                                              4d25s17
             end do                                                     4d25s17
             i10=1                                                      4d25s17
            end do                                                      4d25s17
           end if                                                       4d25s17
           end if                                                       9d25s17
c
c     ddaa part g
c
           if(idoit(4).ne.0)then                                        9d25s17
           if(isblk(1,is).eq.isav.and.isblk(2,is).eq.isav)then          9d7s17
            is34=isblk(3,is)                                            9d7s17
            call ilimts(noc(is34),noc(is34),mynprocg,mynowprog,il,ih,   9d7s17
     $           i1s,i1e,i2s,i2e)                                       9d7s17
            i10=i1s                                                     9d7s17
            i1n=noc(is34)                                               9d7s17
            i4o=ioooo(is)                                               9d7s17
            nrow=(noc(isav)*(noc(isav)+1))/2                            9d7s17
            do i2=i2s,i2e                                               9d7s17
             if(i2.eq.i2e)i1n=i1e                                       9d7s17
             i2m=i2-1-idoub(is34)                                       9d7s17
             do i1=i10,i1n                                              9d7s17
              if(i1.gt.idoub(is34).and.i2.gt.idoub(is34))then           9d7s17
               i1m=i1-1-idoub(is34)                                     9d7s17
               iadd=iden1(is34)+i1m+iacto(is34)*i2m                     9d7s17
               fact=bc(iadd)*4d0                                        9d7s17
               jaa=iaa                                                  9d7s17
               do i4=0,iacto(isav)-1                                    9d7s17
                i4p=i4+idoub(isav)                                      9d7s17
                j4o=i4o+((i4p*(i4p+1))/2)+idoub(isav)                   9d7s17
                do i3=0,i4                                              9d7s17
                 bc(jaa)=bc(jaa)+fact*bc(j4o)                           9d7s17
                 jaa=jaa+1                                              9d7s17
                 j4o=j4o+1                                              9d7s17
                end do                                                  9d7s17
               end do                                                   9d7s17
              end if                                                    9d7s17
              i4o=i4o+nrow                                              9d7s17
             end do                                                     9d7s17
             i10=1                                                      9d7s17
            end do                                                      9d7s17
           end if                                                       9d7s17
c
c     ddaa part h
c
           if(isblk(2,is).eq.isav.and.isblk(4,is).eq.isav)then          9d7s17
            is13=isblk(3,is)                                            9d7s17
            call ilimts(noc(is13),noc(isav),mynprocg,mynowprog,il,ih,   9d7s17
     $           i1s,i1e,i2s,i2e)                                       9d7s17
            i10=i1s                                                     9d7s17
            i1n=noc(is13)                                               9d7s17
            i4o=ioooo(is)                                               9d7s17
            if(isav.eq.is13)then                                        9d7s17
             nrow=(noc(isav)*(noc(isav)+1))/2                            9d7s17
            else                                                        9d7s17
             nrow=noc(isav)*noc(is13)                                   9d7s17
            end if                                                      9d7s17
            do i2=i2s,i2e                                               9d7s17
             if(i2.eq.i2e)i1n=i1e                                       9d7s17
             i2m=i2-1-idoub(isav)                                       9d7s17
             do i1=i10,i1n                                              9d7s17
              if(i1.gt.idoub(is13).and.i2.gt.idoub(isav))then           9d7s17
               i1m=i1-1-idoub(is13)                                     9d7s17
               iadd=iden1(is13)+iacto(is13)*i1m                         9d7s17
               jaa=iaa+((i2m*(i2m+1))/2)
               if(is13.eq.isav)then                                     9d7s17
                do i4=0,min(iacto(isav)-1,i2m)                           9d7s17
                 i4p=i4+idoub(isav)                                      9d7s17
                 kaa=jaa+i4                                             9d7s17
                 do i3=0,iacto(is13)-1                                   9d7s17
                  i3p=i3+idoub(is13)                                    9d7s17
                  ix=max(i3p,i4p)                                        9d7s17
                  in=min(i3p,i4p)                                        9d7s17
                  iade=i4o+((ix*(ix+1))/2)+in                           9d7s17
                  bc(kaa)=bc(kaa)-2d0*bc(iadd+i3)*bc(iade)              9d7s17
                 end do                                                  9d7s17
                end do                                                   9d7s17
               else                                                     9d7s17
                do i4=0,min(iacto(isav)-1,i2m)                           9d7s17
                 i4p=i4+idoub(isav)                                      9d7s17
                 kaa=jaa+i4                                             9d7s17
                 iade=i4o+idoub(is13)+noc(is13)*i4p                     9d7s17
                 do i3=0,iacto(is13)-1                                   9d7s17
                  bc(kaa)=bc(kaa)-2d0*bc(iadd+i3)*bc(iade+i3)           9d7s17
                 end do                                                  9d7s17
                end do                                                   9d7s17
               end if                                                   9d7s17
              end if                                                    9d7s17
              i4o=i4o+nrow                                              9d7s17
             end do                                                     9d7s17
             i10=1                                                      9d7s17
            end do                                                      9d7s17
           else if(isblk(2,is).eq.isav.and.isblk(3,is).eq.isav)then     9d7s17
            is14=isblk(4,is)                                            9d7s17
            call ilimts(noc(isav),noc(is14),mynprocg,mynowprog,il,ih,   9d7s17
     $           i1s,i1e,i2s,i2e)                                       9d7s17
            i10=i1s                                                     9d7s17
            i1n=noc(isav)                                               9d7s17
            i4o=ioooo(is)                                               9d7s17
            nrow=noc(isav)*noc(is14)                                    9d7s17
            do i2=i2s,i2e                                               9d7s17
             if(i2.eq.i2e)i1n=i1e                                       9d7s17
             i2m=i2-1-idoub(is14)                                       9d7s17
             do i1=i10,i1n                                              9d7s17
              if(i2.gt.idoub(is14).and.i1.gt.idoub(isav))then           9d7s17
               i1m=i1-1-idoub(isav)                                     9d7s17
               iadd=iden1(is14)+iacto(is14)*i2m                         9d7s17
               jaa=iaa+((i1m*(i1m+1))/2)
               do i4=0,min(iacto(isav)-1,i1m)                           9d7s17
                i4p=i4+idoub(isav)                                      9d7s17
                kaa=jaa+i4                                              9d7s17
                iade=i4o+idoub(is14)+noc(is14)*i4p                      9d7s17
                do i3=0,iacto(is14)-1                                   9d7s17
                 bc(kaa)=bc(kaa)-2d0*bc(iadd+i3)*bc(iade+i3)            9d7s17
                end do                                                  9d7s17
               end do                                                   9d7s17
              end if                                                    9d7s17
              i4o=i4o+nrow                                              9d7s17
             end do                                                     9d7s17
             i10=1                                                      9d7s17
            end do                                                      9d7s17
           else if(isblk(1,is).eq.isav.and.isblk(3,is).eq.isav)then     9d7s17
            is24=isblk(4,is)                                            9d7s17
            call ilimts(noc(isav),noc(is24),mynprocg,mynowprog,il,ih,   9d7s17
     $           i1s,i1e,i2s,i2e)                                       9d7s17
            i10=i1s                                                     9d7s17
            i1n=noc(isav)                                               9d7s17
            i4o=ioooo(is)                                               9d7s17
            nrow=noc(isav)*noc(is24)                                    9d7s17
            do i2=i2s,i2e                                               9d7s17
             if(i2.eq.i2e)i1n=i1e                                       9d7s17
             i2m=i2-1-idoub(is24)                                       9d7s17
             do i1=i10,i1n                                              9d7s17
              if(i2.gt.idoub(is24).and.i1.gt.idoub(isav))then           9d7s17
               i1m=i1-1-idoub(isav)                                     9d7s17
               iadd=iden1(is24)+iacto(is24)*i2m                         9d7s17
               jaa=iaa+((i1m*(i1m+1))/2)                                9d7s17
               do i4=0,iacto(is24)-1                                    9d7s17
                i4p=i4+idoub(is24)                                      9d7s17
                fact=-2d0*bc(iadd+i4)                                   9d7s17
                iade=i4o+idoub(isav)+noc(isav)*i4p                      9d7s17
                do i3=0,min(iacto(isav)-1,i1m)                          9d7s17
                 bc(jaa+i3)=bc(jaa+i3)+fact*bc(iade+i3)                 9d7s17
                end do                                                  9d7s17
               end do                                                   9d7s17
              end if                                                    9d7s17
              i4o=i4o+nrow                                              9d7s17
             end do                                                     9d7s17
             i10=1                                                      9d7s17
            end do                                                      9d7s17
           else if(isblk(1,is).eq.isav.and.isblk(4,is).eq.isav)then     9d7s17
            is23=isblk(3,is)                                            9d7s17
            call ilimts(noc(is23),noc(isav),mynprocg,mynowprog,il,ih,   9d7s17
     $           i1s,i1e,i2s,i2e)                                       9d7s17
            i10=i1s                                                     9d7s17
            i1n=noc(is23)                                               9d7s17
            i4o=ioooo(is)                                               9d7s17
            nrow=noc(isav)*noc(is23)                                    9d7s17
            do i2=i2s,i2e                                               9d7s17
             if(i2.eq.i2e)i1n=i1e                                       9d7s17
             i2m=i2-1-idoub(isav)                                       9d7s17
             do i1=i10,i1n                                              9d7s17
              if(i2.gt.idoub(isav).and.i1.gt.idoub(is23))then           9d7s17
               i1m=i1-1-idoub(is23)                                     9d7s17
               iadd=iden1(is23)+iacto(is23)*i1m                         9d7s17
               jaa=iaa+((i2m*(i2m+1))/2)                                9d7s17
               do i4=0,iacto(is23)-1                                    9d7s17
                i4p=i4+idoub(is23)                                      9d7s17
                fact=-2d0*bc(iadd+i4)                                   9d7s17
                iade=i4o+idoub(isav)+noc(isav)*i4p                      9d7s17
                do i3=0,min(iacto(isav)-1,i2m)                          9d7s17
                 bc(jaa+i3)=bc(jaa+i3)+fact*bc(iade+i3)                 9d7s17
                end do                                                  9d7s17
               end do                                                   9d7s17
              end if                                                    9d7s17
              i4o=i4o+nrow                                              9d7s17
             end do                                                     9d7s17
             i10=1                                                      9d7s17
            end do                                                      9d7s17
           end if                                                       9d7s17
c
c     ddaa part c&e
c
           if(isblk(1,is).eq.isblk(2,is).and.isblk(3,is).eq.isav)then   9d7s17
            is12=isblk(1,is)                                            9d7s17
            call ilimts(noc(isav),noc(isav),mynprocg,mynowprog,il,ih,   9d7s17
     $           i1s,i1e,i2s,i2e)                                       9d7s17
            i10=i1s                                                     9d7s17
            i1n=noc(isav)                                               9d7s17
            i4o=ioooo(is)-1                                             9d7s17
            nrow=(noc(is12)*(noc(is12)+1))/2                            9d7s17
            do i2=i2s,i2e                                               9d7s17
             i2m=i2-1-idoub(isav)                                       9d7s17
             iadd=iden1(isav)+iacto(isav)*i2m                           9d7s17
             if(i2.eq.i2e)i1n=i1e                                       9d7s17
             do i1=i10,i1n                                              9d7s17
              if(i1.gt.idoub(isav).and.i2.gt.idoub(isav))then           9d7s17
               i1m=i1-1-idoub(isav)                                     9d7s17
               jaa=iaa+((i1m*(i1m+1))/2)
               do i34=1,idoub(is12)
                j4o=i4o+((i34*(i34+1))/2)
                fact=-2d0*bc(j4o)
                do i5=0,i1m
                 bc(jaa+i5)=bc(jaa+i5)+fact*bc(iadd+i5)
                end do
                do i5=i1m,iacto(isav)-1                                 9d13s17
                 kaa=iaa+((i5*(i5+1))/2)+i1m                            9d13s17
                 bc(kaa)=bc(kaa)+fact*bc(iadd+i5)                       9d13s17
                end do                                                  9d13s17
               end do
              end if                                                    9d7s17
              i4o=i4o+nrow                                              9d7s17
             end do                                                     9d7s17
             i10=1                                                      9d7s17
            end do                                                      9d7s17
           end if                                                       9d7s17
c
c     ddaa part d&f
c
           if(isblk(1,is).eq.isblk(3,is).and.isblk(2,is).eq.isav)then   9d7s17
            is13=isblk(1,is)                                            9d7s17
            call ilimts(noc(is13),noc(isav),mynprocg,mynowprog,il,ih,   9d7s17
     $           i1s,i1e,i2s,i2e)                                       9d7s17
            i10=i1s                                                     9d7s17
            i1n=noc(is13)                                               9d7s17
            i4o=ioooo(is)                                               9d7s17
            if(is13.eq.isav)then                                        9d7s17
             nrow=(noc(isav)*(noc(isav)+1))/2                           9d7s17
            else                                                        9d7s17
             nrow=noc(isav)*noc(is13)                                   9d7s17
            end if                                                      9d7s17
            do i2=i2s,i2e                                               9d7s17
             i2m=i2-1-idoub(isav)                                       9d7s17
             jadd=iden1(isav)+iacto(isav)*i2m                           9d7s17
             if(i2.eq.i2e)i1n=i1e                                       9d7s17
             do i1=i10,i1n                                              9d7s17
              if(i1.le.idoub(is13).and.i2.gt.idoub(isav))then           9d7s17
               i1m=i1-1                                                 9d7s17
               jaa=iaa                                                  9d7s17
               if(is13.eq.isav)then                                     9d7s17
                do i3=0,iacto(isav)-1                                   9d7s17
                 i3p=i3+idoub(isav)                                     9d7s17
                 ix=max(i1m,i3p)                                        9d7s17
                 in=min(i1m,i3p)                                        9d7s17
                 j4o2=i4o+((ix*(ix+1))/2)+in                            9d7s17
                 do i4=0,i3                                             9d7s17
                  i4p=i4+idoub(isav)                                    9d7s17
                  iadd=jadd+i3                                          9d7s17
                  iadd2=jadd+i4                                         9d7s17
                  ix=max(i1m,i4p)                                       9d7s17
                  in=min(i1m,i4p)                                       9d7s17
                  j4o=i4o+((ix*(ix+1))/2)+in                            9d7s17
                  bc(jaa)=bc(jaa)+bc(iadd)*bc(j4o)+bc(iadd2)*bc(j4o2)   9d7s17
                  jaa=jaa+1                                             9d7s17
                 end do                                                 9d7s17
                end do                                                  9d7s17
               else                                                     9d7s17
                do i3=0,iacto(isav)-1                                   9d7s17
                 i3p=i3+idoub(isav)                                     9d7s17
                 j4o2=i4o+i1m+noc(is13)*i3p                             9d7s17
                 do i4=0,i3                                             9d7s17
                  i4p=i4+idoub(isav)                                    9d7s17
                  iadd=jadd+i3                                          9d7s17
                  iadd2=jadd+i4                                         9d7s17
                  j4o=i4o+i1m+noc(is13)*i4p                             9d7s17
                  bc(jaa)=bc(jaa)+bc(iadd)*bc(j4o)+bc(iadd2)*bc(j4o2)   9d7s17
                  jaa=jaa+1                                             9d7s17
                 end do                                                 9d7s17
                end do                                                  9d7s17
               end if                                                   9d7s17
              end if                                                    9d7s17
              i4o=i4o+nrow                                              9d7s17
             end do                                                     9d7s17
             i10=1                                                      9d7s17
            end do                                                      9d7s17
           else if(isblk(1,is).eq.isblk(4,is).and.isblk(2,is).eq.isav)  9d7s17
     $           then                                                   9d7s17
            is14=isblk(1,is)                                            9d7s17
            call ilimts(noc(isav),noc(is14),mynprocg,mynowprog,il,ih,   9d7s17
     $           i1s,i1e,i2s,i2e)                                       9d7s17
            i10=i1s                                                     9d7s17
            i1n=noc(isav)                                               9d7s17
            i4o=ioooo(is)                                               9d7s17
            nrow=noc(isav)*noc(is14)                                    9d7s17
            do i2=i2s,i2e                                               9d7s17
             i2m=i2-1
             if(i2.eq.i2e)i1n=i1e                                       9d7s17
             do i1=i10,i1n                                              9d7s17
              if(i1.gt.idoub(isav).and.i2.le.idoub(is14))then           9d7s17
               i1m=i1-1-idoub(isav)                                     9d7s17
               jaa=iaa                                                  9d7s17
               jadd=iden1(isav)+iacto(isav)*i1m                           9d7s17
               do i3=0,iacto(isav)-1                                    9d7s17
                i3p=i3+idoub(isav)                                      9d7s17
                j4o2=i4o+i2m+noc(is14)*i3p                              9d7s17
                do i4=0,i3                                              9d7s17
                 i4p=i4+idoub(isav)                                     9d7s17
                 iadd=jadd+i3                                           9d7s17
                 iadd2=jadd+i4                                          9d7s17
                 j4o=i4o+i2m+noc(is14)*i4p                              9d7s17
                 bc(jaa)=bc(jaa)+bc(iadd)*bc(j4o)+bc(iadd2)*bc(j4o2)    9d7s17
                 jaa=jaa+1                                              9d7s17
                end do                                                  9d7s17
               end do                                                   9d7s17
              end if                                                    9d7s17
              i4o=i4o+nrow                                              9d7s17
             end do                                                     9d7s17
             i10=1                                                      9d7s17
            end do                                                      9d7s17
           else if(isblk(2,is).eq.isblk(4,is).and.isblk(1,is).eq.isav)  9d7s17
     $           then                                                   9d7s17
            is24=isblk(2,is)                                            9d7s17
            call ilimts(noc(isav),noc(is24),mynprocg,mynowprog,il,ih,   9d7s17
     $           i1s,i1e,i2s,i2e)                                       9d7s17
            i10=i1s                                                     9d7s17
            i1n=noc(isav)                                               9d7s17
            i4o=ioooo(is)                                               9d7s17
            nrow=noc(isav)*noc(is24)                                    9d7s17
            do i2=i2s,i2e                                               9d7s17
             i2m=i2-1
             if(i2.eq.i2e)i1n=i1e                                       9d7s17
             do i1=i10,i1n                                              9d7s17
              if(i1.gt.idoub(isav).and.i2.le.idoub(is24))then           9d7s17
               i1m=i1-1-idoub(isav)                                     9d7s17
               jaa=iaa                                                  9d7s17
               jadd=iden1(isav)+iacto(isav)*i1m                           9d7s17
               do i3=0,iacto(isav)-1                                    9d7s17
                i3p=i3+idoub(isav)                                      9d7s17
                j4o2=i4o+i3p+i2m*noc(isav)                              9d7s17
                do i4=0,i3                                              9d7s17
                 i4p=i4+idoub(isav)                                     9d7s17
                 iadd=jadd+i3                                           9d7s17
                 iadd2=jadd+i4                                          9d7s17
                 j4o=i4o+i4p+i2m*noc(isav)                              9d7s17
                 bc(jaa)=bc(jaa)+bc(iadd)*bc(j4o)+bc(iadd2)*bc(j4o2)    9d7s17
                 jaa=jaa+1                                              9d7s17
                end do                                                  9d7s17
               end do                                                   9d7s17
              end if                                                    9d7s17
              i4o=i4o+nrow                                              9d7s17
             end do                                                     9d7s17
             i10=1                                                      9d7s17
            end do                                                      9d7s17
           else if(isblk(2,is).eq.isblk(3,is).and.isblk(1,is).eq.isav)  9d7s17
     $           then                                                   9d7s17
            is23=isblk(2,is)                                            9d7s17
            call ilimts(noc(is23),noc(isav),mynprocg,mynowprog,il,ih,   9d7s17
     $           i1s,i1e,i2s,i2e)                                       9d7s17
            i10=i1s                                                     9d7s17
            i1n=noc(is23)                                               9d7s17
            i4o=ioooo(is)                                               9d7s17
            nrow=noc(isav)*noc(is23)                                    9d7s17
            do i2=i2s,i2e                                               9d7s17
             i2m=i2-1-idoub(isav)                                       9d7s17
             jadd=iden1(isav)+iacto(isav)*i2m                           9d7s17
             if(i2.eq.i2e)i1n=i1e                                       9d7s17
             do i1=i10,i1n                                              9d7s17
              if(i1.le.idoub(is23).and.i2.gt.idoub(isav))then           9d7s17
               i1m=i1-1
               jaa=iaa                                                  9d7s17
               do i3=0,iacto(isav)-1                                    9d7s17
                i3p=i3+idoub(isav)                                      9d7s17
                j4o2=i4o+i3p+i1m*noc(isav)                              9d7s17
                do i4=0,i3                                              9d7s17
                 i4p=i4+idoub(isav)                                     9d7s17
                 iadd=jadd+i3                                           9d7s17
                 iadd2=jadd+i4                                          9d7s17
                 j4o=i4o+i4p+i1m*noc(isav)                              9d7s17
                 bc(jaa)=bc(jaa)+bc(iadd)*bc(j4o)+bc(iadd2)*bc(j4o2)    9d7s17
                 jaa=jaa+1                                              9d7s17
                end do                                                  9d7s17
               end do                                                   9d7s17
              end if                                                    9d7s17
              i4o=i4o+nrow                                              9d7s17
             end do                                                     9d7s17
             i10=1                                                      9d7s17
            end do                                                      9d7s17
           end if                                                       9d7s17
           end if                                                       9d25s17
c
c     dddd delta dd, K part
c
           if(idoit(3).ne.0)then                                        9d25s17
           igot=0                                                       4d25s17
           if(isblk(1,is).eq.isav.and.isblk(3,is).eq.isav)then          4d25s17
            igot=1
            is24=isblk(2,is)                                            4d25s17
            call ilimts(noc(isav),noc(is24),mynprocg,mynowprog,il,ih,   4d25s17
     $           i1s,i1e,i2s,i2e)                                       4d25s17
            i10=i1s                                                     4d25s17
            i1n=noc(isav)                                               4d25s17
            i4o=ioooo(is)                                               4d25s17
            if(isav.eq.is24)then                                        4d25s17
             nrow=(noc(isav)*(noc(isav)+1))/2                           4d25s17
            else                                                        4d25s17
             nrow=noc(isav)*noc(is24)                                   4d25s17
            end if                                                      4d25s17
            do i2=i2s,i2e                                               4d25s17
             if(i2.eq.i2e)i1n=i1e                                       4d25s17
             i2m=i2-1                                                   4d25s17
             do i1=i10,i1n                                              4d25s17
              i1m=i1-1-idoub(isav)                                      4d25s17
              if(i1.gt.idoub(isav).and.i2.le.idoub(is24))then                    4d25s17
               if(isav.eq.is24)then                                     4d25s17
                do i4=0,iacto(isav)-1                                    4d25s17
                 i4p=i4+idoub(isav)                                     4d25s17
                 ix=max(i4,i1m)                                         4d25s17
                 in=min(i4,i1m)                                         4d25s17
                 ix=((ix*(ix+1))/2)+in                                  9d13s17
                 iada=iaa+ix                                            9d13s17
                 iadaf=iaaf+ix                                          9d13s17
                 ix=max(i4p,i2m)                                        4d25s17
                 in=min(i4p,i2m)                                        4d25s17
                 iado=i4o+((ix*(ix+1))/2)+in                            4d25s17
                 bc(iada)=bc(iada)-2d0*bc(iado)*bc(iadaf)               9d13s17
                end do                                                   4d25s17
               else                                                     4d25s17
                do i4=0,iacto(isav)-1                                   4d25s17
                 i4p=i4+idoub(isav)                                     4d25s17
                 ix=max(i4,i1m)                                         4d25s17
                 in=min(i4,i1m)                                         4d25s17
                 ix=((ix*(ix+1))/2)+in                                  9d13s17
                 iada=iaa+ix                                            9d13s17
                 iadaf=iaaf+ix                                          9d13s17
                 iado=i4o+i4p+noc(isav)*i2m                              4d25s17
                 bc(iada)=bc(iada)-2d0*bc(iado)*bc(iadaf)               9d13s17
                end do                                                  4d25s17
               end if                                                   4d25s17
              end if                                                    4d25s17
              i4o=i4o+nrow                                              4d25s17
             end do                                                     4d25s17
             i10=1                                                      4d25s17
            end do                                                      4d25s17
           end if                                                       4d25s17
           if(isblk(1,is).eq.isav.and.isblk(4,is).eq.isav.and.          4d25s17
     $          igot.eq.0)then                                          4d25s17
            igot=1
            is23=isblk(2,is)                                            4d25s17
            call ilimts(noc(is23),noc(isav),mynprocg,mynowprog,il,ih,   4d25s17
     $           i1s,i1e,i2s,i2e)                                       4d25s17
            i10=i1s                                                     4d25s17
            i1n=noc(is23)                                               4d25s17
            i4o=ioooo(is)                                               4d25s17
            nrow=noc(isav)*noc(is23)                                    4d25s17
            do i2=i2s,i2e                                               4d25s17
             if(i2.eq.i2e)i1n=i1e                                       4d25s17
             i2m=i2-1-idoub(isav)                                       4d25s17
             do i1=i10,i1n                                              4d25s17
              i1m=i1-1                                                  4d25s17
              if(i2.gt.idoub(isav).and.i1.le.idoub(is23))then           4d25s17
               do i4=0,iacto(isav)-1                                    4d25s17
                i4p=i4+idoub(isav)                                      4d25s17
                ix=max(i4,i2m)                                          4d25s17
                in=min(i4,i2m)                                          4d25s17
                ix=((ix*(ix+1))/2)+in                                   9d13s17
                iada=iaa+ix                                             9d13s17
                iadaf=iaaf+ix                                           9d13s17
                iado=i4o+i4p+noc(isav)*i1m                              4d25s17
                bc(iada)=bc(iada)-2d0*bc(iado)*bc(iadaf)                9d13s17
               end do                                                   4d25s17
              end if                                                    4d25s17
              i4o=i4o+nrow                                              4d25s17
             end do                                                     4d25s17
             i10=1                                                      4d25s17
            end do                                                      4d25s17
           end if                                                       4d25s17
           if(isblk(2,is).eq.isav.and.isblk(4,is).eq.isav.and.          4d25s17
     $          igot.eq.0)then                                          4d25s17
            igot=1
            is13=isblk(1,is)                                            4d25s17
            call ilimts(noc(is13),noc(isav),mynprocg,mynowprog,il,ih,   4d25s17
     $           i1s,i1e,i2s,i2e)                                       4d25s17
            i10=i1s                                                     4d25s17
            i1n=noc(is13)                                               4d25s17
            i4o=ioooo(is)                                               4d25s17
            nrow=noc(isav)*noc(is13)                                    4d25s17
            do i2=i2s,i2e                                               4d25s17
             if(i2.eq.i2e)i1n=i1e                                       4d25s17
             i2m=i2-1-idoub(isav)                                       4d25s17
             do i1=i10,i1n                                              4d25s17
              i1m=i1-1                                                  4d25s17
              if(i2.gt.idoub(isav).and.i1.le.idoub(is13))then           4d25s17
               do i4=0,iacto(isav)-1                                    4d25s17
                i4p=i4+idoub(isav)                                      4d25s17
                ix=max(i4,i2m)                                          4d25s17
                in=min(i4,i2m)                                          4d25s17
                ix=((ix*(ix+1))/2)+in                                   9d13s17
                iada=iaa+ix                                             9d13s17
                iadaf=iaaf+ix                                           9d13s17
                iado=i4o+i1m+noc(is13)*i4p                              4d25s17
                bc(iada)=bc(iada)-2d0*bc(iado)*bc(iadaf)                9d13s17
               end do                                                   4d25s17
              end if                                                    4d25s17
              i4o=i4o+nrow                                              4d25s17
             end do                                                     4d25s17
             i10=1                                                      4d25s17
            end do                                                      4d25s17
           end if                                                       4d25s17
           if(isblk(2,is).eq.isav.and.isblk(3,is).eq.isav.and.          4d25s17
     $          igot.eq.0)then                                          4d25s17
            igot=1
            is14=isblk(1,is)                                            4d25s17
            call ilimts(noc(isav),noc(is14),mynprocg,mynowprog,il,ih,   4d25s17
     $           i1s,i1e,i2s,i2e)                                       4d25s17
            i10=i1s                                                     4d25s17
            i1n=noc(isav)                                               4d25s17
            i4o=ioooo(is)                                               4d25s17
            nrow=noc(isav)*noc(is14)                                    4d25s17
            do i2=i2s,i2e                                               4d25s17
             if(i2.eq.i2e)i1n=i1e                                       4d25s17
             i2m=i2-1                                                   4d25s17
             do i1=i10,i1n                                              4d25s17
              i1m=i1-1-idoub(isav)                                      4d25s17
              if(i1.gt.idoub(isav).and.i2.le.idoub(is14))then           4d25s17
               do i4=0,iacto(isav)-1                                    4d25s17
                i4p=i4+idoub(isav)                                      4d25s17
                ix=max(i4,i1m)                                          4d25s17
                in=min(i4,i1m)                                          4d25s17
                ix=((ix*(ix+1))/2)+in                                   9d13s17
                iada=iaa+ix                                             9d13s17
                iadaf=iaaf+ix                                           9d13s17
                iado=i4o+i2m+noc(is14)*i4p                              4d25s17
                bc(iada)=bc(iada)-2d0*bc(iado)*bc(iadaf)                9d13s17
               end do                                                   4d25s17
              end if                                                    4d25s17
              i4o=i4o+nrow                                              4d25s17
             end do                                                     4d25s17
             i10=1                                                      4d25s17
            end do                                                      4d25s17
           end if                                                       4d25s17
           end if                                                       9d25s17
          end do                                                        4d25s17
c                                                                       7d21s17
c     1 part density parts of hessa                                     7d21s17
c
          igot=0                                                        7d21s17
          do is=1,nsdlk                                                 7d21s17
           if(idoub(isblk(1,is)).ne.0)then                              7d21s17
           end if                                                       7d21s17
          end do                                                        7d21s17
          nsum=ndd*2+naa                                                9d5s17
          call dws_gsumf(bc(idd),nsum)                                  9d5s17
c
c     delta dd h0 part
c
          if(idoit(1).ne.0)then                                         9d25s17
          jaa=iaa                                                       4d25s17
          nbas=nbasdws(isav)*nvirtc(isav)/nvirt(isav)                   9d25s17
          do i4=0,iacto(isav)-1                                         4d25s17
           i4p=i4+idoub(isav)                                           4d25s17
           do i3=0,i4                                                   4d25s17
            i3p=i3+idoub(isav)+1                                        4d25s17
            iad=ipt(isav)+i3p+nbas*i4p                                  4d25s17
            bc(jaa)=bc(jaa)+4d0*h0mo(iad)                               4d25s17
            jaa=jaa+1                                                   4d25s17
           end do                                                       4d25s17
          end do                                                        4d25s17
          end if                                                        9d25s17
c
c     aa contribution to iaa
c
          if(idoit(2).ne.0)then                                         9d25s17
          jaa=iaa                                                       9d14s17
          do i4=0,iacto(isav)-1                                         9d14s17
           iadd4=iden1(isav)+iacto(isav)*i4                             9d14s17
           i4p=i4+idoub(isav)                                           9d14s17
           ih3=ipt(isav)+1+nbas*i4p+idoub(isav)                         9d14s17
           do i3=0,i4                                                   9d14s17
            i3p=i3+idoub(isav)                                          9d14s17
            ih4=ipt(isav)+1+nbas*i3p+idoub(isav)                        9d14s17
            iadd3=iden1(isav)+iacto(isav)*i3                            9d14s17
            do i2=0,iacto(isav)-1                                       9d14s17
             bc(jaa)=bc(jaa)-bc(iadd4+i2)*h0mo(ih4+i2)                  9d14s17
     $            -bc(iadd3+i2)*h0mo(ih3+i2)                            9d14s17
            end do                                                      9d14s17
            jaa=jaa+1                                                   9d14s17
           end do                                                       9d14s17
          end do                                                        9d14s17
          end if                                                        9d25s17
c
          do i2=0,iacto(isav)-1                                         4d25s17
           do i1=0,iacto(isav)-1                                        4d25s17
            ix=max(i1,i2)                                               4d25s17
            in=min(i1,i2)                                               4d25s17
            iada=iaa+((ix*(ix+1))/2)+in                                 4d25s17
            do i34=0,idoub(isao)-1                                      4d25s17
             iad=ihessa(isa,isb)+i34*(idoub(isao)+1)                    4d25s17
     $            +nrowa*(i1+iacto(isav)*i2)                            4d25s17
             bc(iad)=bc(iad)+bc(iada)                                   4d25s17
            end do                                                      4d25s17
           end do                                                       4d25s17
          end do                                                        4d25s17
c
c     delta nn part
c
          if(idoit(1).ne.0)then                                         9d25s17
          jdd=idd                                                       4d25s17
          nbas=nbasdws(isao)*nvirtc(isao)/nvirt(isao)                   9d25s17
          do i4=0,idoub(isao)-1                                         4d25s17
           do i3=0,i4                                                   4d25s17
            iad=ipt(isao)+i3+1+nbas*i4                                  4d25s17
            bc(jdd)=bc(jdd)-4d0*h0mo(iad)                               4d25s17
            jdd=jdd+1                                                   4d25s17
           end do                                                       4d25s17
          end do                                                        4d25s17
          end if                                                        9d25s17
          do i12=0,iacto(isav)-1                                        4d25s17
           jad=ihessa(isa,isb)+nrowa*i12*(iacto(isav)+1)                 4d25s17
           do i4=0,idoub(isao)-1                                        4d25s17
            do i3=0,idoub(isao)-1                                       4d25s17
             ix=max(i3,i4)                                              4d25s17
             in=min(i3,i4)                                              4d25s17
             iadd=idd+((ix*(ix+1))/2)+in                                4d25s17
             iad=jad+i3+idoub(isao)*i4                                  4d25s17
             bc(iad)=bc(iad)+bc(iadd)                                   4d25s17
            end do                                                      4d25s17
           end do                                                       4d25s17
          end do                                                        4d25s17
c                                                                       9d14s17
c     aa contribution to idd2                                           9d14s17
c                                                                       9d14s17
          if(idoit(2).ne.0)then                                         9d25s17
          jdd=idd2                                                      9d14s17
          do i4=0,idoub(isao)-1                                         9d14s17
           do i3=0,i4                                                   9d14s17
            iad=ipt(isao)+i3+1+nbas*i4                                  9d14s17
            bc(jdd)=bc(jdd)+2d0*h0mo(iad)                               9d14s17
            jdd=jdd+1                                                   9d14s17
           end do                                                       9d14s17
          end do                                                        9d14s17
          end if                                                        9d25s17
c
c     ddaa parts i and j
c
          jad=ihessa(isa,isb)                                           9d5s17
          do i2=0,iacto(isav)-1                                         9d5s17
           do i1=0,iacto(isav)-1                                        9d5s17
            iad1=iden1(isav)+i1+iacto(isav)*i2                          9d5s17
            do i4=0,idoub(isao)-1
             do i3=0,idoub(isao)-1                                      9d5s17
              ix=max(i3,i4)                                             9d5s17
              in=min(i3,i4)                                             9d5s17
              iadd=idd2+((ix*(ix+1))/2)+in                              9d5s17
              bc(jad)=bc(jad)+bc(iadd)*bc(iad1)                         9d5s17
              jad=jad+1                                                 9d5s17
             end do                                                     9d5s17
            end do                                                      9d5s17
           end do                                                       9d5s17
          end do                                                        9d5s17
          ibcoff=idd                                                    4d25s17
c
c     for doub&act vs doub and virt
c
          iav=ibcoff                                                    4d25s17
          nwds=iacto(isa)*nvirtc(isb)                                   4d25s17
          if(nwds.gt.0)then                                             4d25s17
           ibcoff=iav+nwds                                              4d25s17
           call enough('buildhesscas.  7',bc,ibc)
           do i=0,nwds-1                                                4d25s17
            bc(iav+i)=0d0                                               4d25s17
           end do                                                       4d25s17
           do is=1,nsdlk1                                               4d25s17
c
c     aaaa nd,vd part
c
            if(idoit(5).ne.0)then                                       9d25s17
            if(isblk1(4,is).eq.isbv)then                                9d25s17
             is3=isblk1(3,is)                                           9d25s17
             is1=isblk1(1,is)                                           9d25s17
             is2=isblk1(2,is)                                           9d25s17
             iswitch=iabs(is1-is2)                                      9d25s17
             iswitch=iswitch/max(1,iswitch)                             9d25s17
             call ilimts(noc(is3),nvirtc(isbv),mynprocg,mynowprog,il,   9d25s17
     $            ih,i1s,i1e,i2s,i2e)                                   9d25s17
             if(is1.eq.is2)then                                         9d25s17
              nrow=(noc(is1)*(noc(is1)+1))/2                            9d25s17
             else                                                       9d25s17
              nrow=noc(is1)*noc(is2)                                    9d25s17
             end if                                                     9d25s17
             do isd=1,nsdlk                                             9d25s17
              if(isblk(1,isd).eq.isblk(2,isd))then                      9d25s17
               fmult=1d0                                                9d25s17
              else                                                      9d25s17
               fmult=0.5d0                                              9d25s17
              end if                                                    9d25s17
              if(isblk(1,isd).eq.isa.and.isblk(2,isd).eq.is3.and.       9d25s17
     $             isblk(3,isd).eq.is1.and.jdenpt(isd).gt.0)then        9d25s17
               jswitch=iabs(isa-is3)                                    9d25s17
               jswitch=jswitch/max(1,jswitch)                           9d25s17
               kswitch=iabs(is1-is2)                                    9d25s17
               kswitch=kswitch/max(1,kswitch)                           9d25s17
               if(isa.eq.is3)then                                       9d25s17
                nrowd=(iacto(isa)*(iacto(isa)+1))/2                     9d25s17
               else                                                     9d25s17
                nrowd=iacto(isa)*iacto(is3)                             9d25s17
               end if                                                   9d25s17
               i10=i1s                                                    9d25s17
               i1n=noc(is3)                                               9d25s17
               ix1=ionex(is)                                              9d25s17
               do i2=i2s,i2e                                              9d25s17
                if(i2.eq.i2e)i1n=i1e                                      9d25s17
                i2m=i2-1                                                9d25s17
                jav=iav+iacto(isa)*i2m                                  9d25s17
                do i1=i10,i1n                                             9d25s17
                 if(i1.gt.idoub(is3))then                               9d25s17
                  i1m=i1-1-idoub(is3)                                   9d25s17
                  do i4=0,iacto(is2)-1                                  9d25s17
                   i4p=i4+idoub(is2)                                    9d25s17
                   do i3=0,iacto(is1)-1                                 9d25s17
                    i3p=i3+idoub(is1)                                   9d25s17
                    ix=max(i3p,i4p)                                     9d25s17
                    in=min(i3p,i4p)                                     9d25s17
                    ieq=((ix*(ix+1))/2)+in                              9d25s17
                    inot=i3p+noc(is1)*i4p                               9d25s17
                    iad1=ix1+(inot-ieq)*iswitch+ieq                     9d25s17
                    fact=-bc(iad1)*fmult                                9d25s17
                    ix=max(i3,i4)                                       9d25s17
                    in=min(i3,i4)                                       9d25s17
                    ieq=((ix*(ix+1))/2)+in                              9d25s17
                    inot=i3+iacto(is1)*i4                               9d25s17
                    icol=(inot-ieq)*kswitch+ieq                         9d25s17
                    iden2=jdenpt(isd)+nrowd*icol                        9d25s17
                    do i5=0,iacto(isa)-1                                9d25s17
                     ix=max(i5,i1m)                                     9d25s17
                     in=min(i5,i1m)                                     9d25s17
                     ieq=((ix*(ix+1))/2)+in                             9d25s17
                     inot=i5+iacto(isa)*i1m                             9d25s17
                     jad=iden2+(inot-ieq)*jswitch+ieq                   9d25s17
                     bc(jav+i5)=bc(jav+i5)+fact*bc(jad)                 9d25s17
                    end do                                              9d25s17
                   end do                                               9d25s17
                  end do                                                9d25s17
                 end if                                                 9d25s17
                 ix1=ix1+nrow                                             9d25s17
                end do                                                    9d25s17
                i10=1                                                     9d25s17
               end do                                                     9d25s17
              else if(isblk(2,isd).eq.isa.and.isblk(1,isd).eq.is3.and.  9d25s17
     $             isblk(3,isd).eq.is1.and.jdenpt(isd).gt.0)then        9d25s17
               jswitch=iabs(isa-is3)                                    9d25s17
               jswitch=jswitch/max(1,jswitch)                           9d25s17
               kswitch=iabs(is1-is2)                                    9d25s17
               kswitch=kswitch/max(1,kswitch)                           9d25s17
               if(isa.eq.is3)then                                       9d25s17
                nrowd=(iacto(isa)*(iacto(isa)+1))/2                     9d25s17
               else                                                     9d25s17
                nrowd=iacto(isa)*iacto(is3)                             9d25s17
               end if                                                   9d25s17
               i10=i1s                                                    9d25s17
               i1n=noc(is3)                                               9d25s17
               ix1=ionex(is)                                              9d25s17
               do i2=i2s,i2e                                              9d25s17
                if(i2.eq.i2e)i1n=i1e                                      9d25s17
                i2m=i2-1                                                9d25s17
                jav=iav+iacto(isa)*i2m                                  9d25s17
                do i1=i10,i1n                                             9d25s17
                 if(i1.gt.idoub(is3))then                               9d25s17
                  i1m=i1-1-idoub(is3)                                   9d25s17
                  do i4=0,iacto(is2)-1                                  9d25s17
                   i4p=i4+idoub(is2)                                    9d25s17
                   do i3=0,iacto(is1)-1                                 9d25s17
                    i3p=i3+idoub(is1)                                   9d25s17
                    ix=max(i3p,i4p)                                     9d25s17
                    in=min(i3p,i4p)                                     9d25s17
                    ieq=((ix*(ix+1))/2)+in                              9d25s17
                    inot=i3p+noc(is1)*i4p                               9d25s17
                    iad1=ix1+(inot-ieq)*iswitch+ieq                     9d25s17
                    fact=-bc(iad1)*fmult                                9d25s17
                    ix=max(i3,i4)                                       9d25s17
                    in=min(i3,i4)                                       9d25s17
                    ieq=((ix*(ix+1))/2)+in                              9d25s17
                    inot=i3+iacto(is1)*i4                               9d25s17
                    icol=(inot-ieq)*kswitch+ieq                         9d25s17
                    iden2=jdenpt(isd)+nrowd*icol                        9d25s17
                    do i5=0,iacto(isa)-1                                9d25s17
                     inot=i1m+iacto(is3)*i5                             9d25s17
                     jad=iden2+inot                                     9d25s17
                     bc(jav+i5)=bc(jav+i5)+fact*bc(jad)                 9d25s17
                    end do                                              9d25s17
                   end do                                               9d25s17
                  end do                                                9d25s17
                 end if                                                 9d25s17
                 ix1=ix1+nrow                                             9d25s17
                end do                                                    9d25s17
                i10=1                                                     9d25s17
               end do                                                     9d25s17
              else if(isblk(2,isd).eq.isa.and.isblk(1,isd).eq.is3.and.  9d25s17
     $             isblk(3,isd).eq.is2.and.jdenpt(isd).gt.0)then        9d25s17
               jswitch=iabs(isa-is3)                                    9d25s17
               jswitch=jswitch/max(1,jswitch)                           9d25s17
               kswitch=iabs(is1-is2)                                    9d25s17
               kswitch=kswitch/max(1,kswitch)                           9d25s17
               if(isa.eq.is3)then                                       9d25s17
                nrowd=(iacto(isa)*(iacto(isa)+1))/2                     9d25s17
               else                                                     9d25s17
                nrowd=iacto(isa)*iacto(is3)                             9d25s17
               end if                                                   9d25s17
               i10=i1s                                                    9d25s17
               i1n=noc(is3)                                               9d25s17
               ix1=ionex(is)                                              9d25s17
               do i2=i2s,i2e                                              9d25s17
                if(i2.eq.i2e)i1n=i1e                                      9d25s17
                i2m=i2-1                                                9d25s17
                jav=iav+iacto(isa)*i2m                                  9d25s17
                do i1=i10,i1n                                             9d25s17
                 if(i1.gt.idoub(is3))then                               9d25s17
                  i1m=i1-1-idoub(is3)                                   9d25s17
                  do i4=0,iacto(is2)-1                                  9d25s17
                   i4p=i4+idoub(is2)                                    9d25s17
                   do i3=0,iacto(is1)-1                                 9d25s17
                    i3p=i3+idoub(is1)                                   9d25s17
                    ix=max(i3p,i4p)                                     9d25s17
                    in=min(i3p,i4p)                                     9d25s17
                    ieq=((ix*(ix+1))/2)+in                              9d25s17
                    inot=i3p+noc(is1)*i4p                               9d25s17
                    iad1=ix1+(inot-ieq)*iswitch+ieq                     9d25s17
                    fact=-bc(iad1)*fmult                                9d25s17
                    icol=i4+iacto(is2)*i3                               9d25s17
                    iden2=jdenpt(isd)+nrowd*icol                        9d25s17
                    do i5=0,iacto(isa)-1                                9d25s17
                     ix=max(i5,i1m)                                     9d25s17
                     in=min(i5,i1m)                                     9d25s17
                     ieq=((ix*(ix+1))/2)+in                             9d25s17
                     inot=i1m+iacto(is3)*i5                             9d25s17
                     jad=iden2+(inot-ieq)*jswitch+ieq                   9d25s17
                     bc(jav+i5)=bc(jav+i5)+fact*bc(jad)                 9d25s17
                    end do                                              9d25s17
                   end do                                               9d25s17
                  end do                                                9d25s17
                 end if                                                 9d25s17
                 ix1=ix1+nrow                                             9d25s17
                end do                                                    9d25s17
                i10=1                                                     9d25s17
               end do                                                     9d25s17
              else if(isblk(1,isd).eq.isa.and.isblk(2,isd).eq.is3.and.  9d25s17
     $             isblk(3,isd).eq.is2.and.jdenpt(isd).gt.0)then        9d25s17
               jswitch=iabs(isa-is3)                                    9d25s17
               jswitch=jswitch/max(1,jswitch)                           9d25s17
               kswitch=iabs(is1-is2)                                    9d25s17
               kswitch=kswitch/max(1,kswitch)                           9d25s17
               if(isa.eq.is3)then                                       9d25s17
                nrowd=(iacto(isa)*(iacto(isa)+1))/2                     9d25s17
               else                                                     9d25s17
                nrowd=iacto(isa)*iacto(is3)                             9d25s17
               end if                                                   9d25s17
               i10=i1s                                                    9d25s17
               i1n=noc(is3)                                               9d25s17
               ix1=ionex(is)                                              9d25s17
               do i2=i2s,i2e                                              9d25s17
                if(i2.eq.i2e)i1n=i1e                                      9d25s17
                i2m=i2-1                                                9d25s17
                jav=iav+iacto(isa)*i2m                                  9d25s17
                do i1=i10,i1n                                             9d25s17
                 if(i1.gt.idoub(is3))then                               9d25s17
                  i1m=i1-1-idoub(is3)                                   9d25s17
                  do i4=0,iacto(is2)-1                                  9d25s17
                   i4p=i4+idoub(is2)                                    9d25s17
                   do i3=0,iacto(is1)-1                                 9d25s17
                    i3p=i3+idoub(is1)                                   9d25s17
                    ix=max(i3p,i4p)                                     9d25s17
                    in=min(i3p,i4p)                                     9d25s17
                    ieq=((ix*(ix+1))/2)+in                              9d25s17
                    inot=i3p+noc(is1)*i4p                               9d25s17
                    iad1=ix1+(inot-ieq)*iswitch+ieq                     9d25s17
                    fact=-bc(iad1)*fmult                                9d25s17
                    icol=i3+iacto(is1)*i4                               9d25s17
                    iden2=jdenpt(isd)+nrowd*icol                        9d25s17
                    do i5=0,iacto(isa)-1                                9d25s17
                     ix=max(i5,i1m)                                     9d25s17
                     in=min(i5,i1m)                                     9d25s17
                     ieq=((ix*(ix+1))/2)+in                             9d25s17
                     inot=i1m+iacto(is3)*i5                             9d25s17
                     jad=iden2+(inot-ieq)*jswitch+ieq                   9d25s17
                     bc(jav+i5)=bc(jav+i5)+fact*bc(jad)                 9d25s17
                    end do                                              9d25s17
                   end do                                               9d25s17
                  end do                                                9d25s17
                 end if                                                 9d25s17
                 ix1=ix1+nrow                                             9d25s17
                end do                                                    9d25s17
                i10=1                                                     9d25s17
               end do                                                     9d25s17
              end if                                                    9d25s17
             end do                                                     9d25s17
            end if                                                      9d25s17
            end if                                                      9d25s17
c
c     ddaa nd,vd part a and c
c
            if(idoit(4).ne.0)then                                       9d25s17
            if(isblk1(4,is).eq.isbv.and.isblk1(3,is).eq.isbv)then       9d22s17
             is12=isblk1(1,is)                                          9d22s17
             call ilimts(noc(isbv),nvirtc(isbv),mynprocg,mynowprog,il,  9d22s17
     $            ih,i1s,i1e,i2s,i2e)                                   9d22s17
             i10=i1s                                                    9d22s17
             i1n=noc(isbv)                                              9d22s17
             ix1=ionex(is)-1                                            9d22s17
             nrow=(noc(is12)*(noc(is12)+1))/2                           9d22s17
             do i2=i2s,i2e                                              9d22s17
              if(i2.eq.i2e)i1n=i1e                                      9d22s17
              i2m=i2-1                                                  9d22s17
              do i1=i10,i1n                                             9d22s17
               if(i1.gt.idoub(isbv))then                                9d22s17
                i1m=i1-1-idoub(isbv)                                    9d22s17
                iden=iden1(isbv)+iacto(isbv)*i1m                        9d22s17
                jav=iav+iacto(isav)*i2m
                do i34=1,idoub(is12)                                    9d22s17
                 iadi=ix1+((i34*(i34+1))/2)                             9d22s17
                 fact=-bc(iadi)                                         9d22s17
                 do i4=0,iacto(isav)-1                                  9d22s17
                  bc(jav+i4)=bc(jav+i4)+fact*bc(iden+i4)                9d22s17
                 end do                                                 9d22s17
                end do
                jav=iav+i1m+iacto(isav)*i2m                             9d22s17
                jx1=ix1+1                                               9d22s17
                do i4=0,iacto(is12)-1                                   9d22s17
                 i4p=i4+idoub(is12)                                     9d22s17
                 iden=iden1(is12)+iacto(is12)*i4                        9d22s17
                 do i3=0,iacto(is12)-1                                  9d22s17
                  i3p=i3+idoub(is12)                                    9d22s17
                  ix=max(i3p,i4p)                                       9d22s17
                  in=min(i3p,i4p)                                       9d22s17
                  iadi=jx1+((ix*(ix+1))/2)+in                           9d22s17
                  bc(jav)=bc(jav)+2d0*bc(iadi)*bc(iden+i3)              9d22s17
                 end do                                                 9d22s17
                end do                                                  9d22s17
               end if                                                   9d22s17
               ix1=ix1+nrow                                             9d22s17
              end do                                                    9d22s17
              i10=1                                                     9d22s17
             end do                                                     9d22s17
            end if                                                      9d22s17
c
c     ddaa nd,vd part b and d
c
            if(isblk1(4,is).eq.isbv.and.                                9d22s17
     $           isblk1(3,is).eq.isblk1(1,is))then                      9d22s17
             is13=isblk1(1,is)                                          9d22s17
             call ilimts(noc(is13),nvirtc(isbv),mynprocg,mynowprog,il,  9d22s17
     $            ih,i1s,i1e,i2s,i2e)                                   9d22s17
             i10=i1s                                                    9d22s17
             i1n=noc(is13)                                              9d22s17
             if(is13.eq.isbv)then                                       9d22s17
              iswitch=0                                                 9d22s17
              nrow=(noc(is13)*(noc(is13)+1))/2                          9d22s17
             else                                                       9d22s17
              iswitch=1                                                 9d22s17
              nrow=noc(is13)*noc(isbv)                                  9d22s17
             end if                                                     9d22s17
             ix1=ionex(is)                                              9d22s17
             do i2=i2s,i2e                                              9d22s17
              if(i2.eq.i2e)i1n=i1e                                      9d22s17
              i2m=i2-1                                                  9d22s17
              do i1=i10,i1n                                             9d22s17
               jav=iav+iacto(isav)*i2m                                  9d22s17
               if(i1.le.idoub(is13))then                                9d22s17
                i1m=i1-1                                                9d22s17
                do i3=0,iacto(isav)-1                                   9d22s17
                 i3p=i3+idoub(isav)                                     9d22s17
                 ix=max(i3p,i1m)                                        9d22s17
                 in=min(i3p,i1m)                                        9d22s17
                 ieq=((ix*(ix+1))/2)+in                                 9d22s17
                 inot=i1m+noc(is13)*i3p                                 9d22s17
                 iadi=ix1+(inot-ieq)*iswitch+ieq                        9d22s17
                 fact=0.5d0*bc(iadi)                                    9d22s17
                 iden=iden1(isav)+iacto(isav)*i3                        9d22s17
                 do i4=0,iacto(isav)-1                                  9d22s17
                  bc(jav+i4)=bc(jav+i4)+bc(iden+i4)*fact                9d22s17
                 end do                                                 9d22s17
                end do                                                  9d22s17
               else                                                     9d22s17
                i1m=i1-1-idoub(is13)                                    9d22s17
                iden=iden1(is13)+iacto(is13)*i1m                        9d22s17
                do i4=0,iacto(isbv)-1                                   9d22s17
                 i4p=i4+idoub(isbv)                                     9d22s17
                 kav=jav+i4                                             9d22s17
                 do i3=0,iacto(is13)-1                                  9d22s17
                  i3p=i3+idoub(is13)                                    9d22s17
                  ix=max(i4p,i3p)                                       9d22s17
                  in=min(i4p,i3p)                                       9d22s17
                  ieq=((ix*(ix+1))/2)+in                                9d22s17
                  inot=i3p+noc(is13)*i4p                                9d22s17
                  iadi=ix1+(inot-ieq)*iswitch+ieq                       9d22s17
                  bc(kav)=bc(kav)-bc(iadi)*bc(iden+i3)                  9d22s17
                 end do                                                 9d22s17
                end do                                                  9d22s17
               end if                                                   9d22s17
               ix1=ix1+nrow                                             9d22s17
              end do                                                    9d22s17
              i10=1                                                     9d22s17
             end do                                                     9d22s17
            else if(isblk1(4,is).eq.isbv.and.                           9d22s17
     $           isblk1(3,is).eq.isblk1(2,is))then                      9d22s17
             is23=isblk1(2,is)                                          9d22s17
             call ilimts(noc(is23),nvirtc(isbv),mynprocg,mynowprog,il,  9d22s17
     $            ih,i1s,i1e,i2s,i2e)                                   9d22s17
             i10=i1s                                                    9d22s17
             i1n=noc(is23)                                              9d22s17
             nrow=noc(is23)*noc(isbv)                                   9d22s17
             ix1=ionex(is)                                              9d22s17
             do i2=i2s,i2e                                              9d22s17
              if(i2.eq.i2e)i1n=i1e                                      9d22s17
              i2m=i2-1                                                  9d22s17
              do i1=i10,i1n                                             9d22s17
               jav=iav+iacto(isav)*i2m                                  9d22s17
               if(i1.le.idoub(is23))then                                9d22s17
                i1m=i1-1                                                9d22s17
                do i3=0,iacto(isav)-1                                   9d22s17
                 i3p=i3+idoub(isav)                                     9d22s17
                 inot=i3p+noc(isbv)*i1m                                 9d22s17
                 iadi=ix1+inot                                          9d22s17
                 fact=bc(iadi)*0.5d0                                    9d22s17
                 iden=iden1(isav)+iacto(isav)*i3                        9d22s17
                 do i4=0,iacto(isav)-1                                  9d22s17
                  bc(jav+i4)=bc(jav+i4)+bc(iden+i4)*fact                9d22s17
                 end do                                                 9d22s17
                end do                                                  9d22s17
               else                                                     9d22s17
                i1m=i1-1-idoub(is23)                                    9d22s17
                iden=iden1(is23)+iacto(is23)*i1m                        9d22s17
                do i4=0,iacto(isbv)-1                                   9d22s17
                 i4p=i4+idoub(isbv)                                     9d22s17
                 kav=jav+i4                                             9d22s17
                 do i3=0,iacto(is23)-1                                  9d22s17
                  i3p=i3+idoub(is23)                                    9d22s17
                  inot=i4p+noc(isbv)*i3p                                9d22s17
                  iadi=ix1+inot                                         9d22s17
                  bc(kav)=bc(kav)-bc(iadi)*bc(iden+i3)                  9d22s17
                 end do                                                 9d22s17
                end do                                                  9d22s17
               end if                                                   9d22s17
               ix1=ix1+nrow                                             9d22s17
              end do                                                    9d22s17
              i10=1                                                     9d22s17
             end do                                                     9d22s17
            end if                                                      9d22s17
            end if                                                      9d25s17
c
c     dddd: "K" part
c
            if(idoit(3).ne.0)then                                       9d25s17
            if(isblk1(1,is).eq.isblk1(2,is).and.isblk1(3,is).eq.isa     4d25s17
     $           .and.isblk1(4,is).eq.isb)then                          4d25s17
             is12=isblk1(1,is)                                          4d25s17
             call ilimts(noc(isa),nvirtc(isb),mynprocg,mynowprog,il,ih, 4d25s17
     $           i1s,i1e,i2s,i2e)                                       4d25s17
             i10=i1s                                                    4d25s17
             i1n=noc(isa)                                               4d25s17
             ixo=ionex(is)-1                                            6d19s17
             nrow=(noc(is12)*(noc(is12)+1))/2                           4d25s17
             do i2=i2s,i2e                                              4d25s17
              i2m=i2-1                                                  4d25s17
              if(i2.eq.i2e)i1n=i1e                                      4d25s17
              do i1=i10,i1n                                             4d25s17
               i1m=i1-1-idoub(isa)                                      4d25s17
               if(i1.gt.idoub(isa))then                                 4d25s17
                do i34=0,idoub(is12)-1                                  4d25s17
                 i34p=i34+1                                             4d25s17
                 iadx=ixo+(i34p*(i34p+1))/2                             4d25s17
                 jav=iav+i1m+iacto(isa)*i2m                             4d25s17
                 bc(jav)=bc(jav)+4d0*bc(iadx)                           4d25s17
                end do                                                  4d25s17
               end if                                                   4d25s17
               ixo=ixo+nrow                                             4d25s17
              end do                                                    4d25s17
              i10=1                                                     4d25s17
             end do                                                     4d25s17
            end if                                                      4d25s17
c
c     dddd: "J" part
c
            igot=0                                                      4d25s17
            if(isblk1(1,is).eq.isa.and.isblk1(2,is).eq.isblk1(3,is)     4d25s17
     $           .and.isblk1(4,is).eq.isb)then                          4d25s17
             igot=1                                                     4d25s17
             is23=isblk1(2,is)                                           4d25s17
             call ilimts(noc(is23),nvirtc(isb),mynprocg,mynowprog,il,ih,4d25s17
     $            i1s,i1e,i2s,i2e)                                      4d25s17
             i10=i1s                                                    4d25s17
             i1n=noc(is23)                                              4d25s17
             if(isa.eq.is23)then                                        4d25s17
              nrow=(noc(isa)*(noc(isa)+1))/2                            4d25s17
             else                                                       4d25s17
              nrow=noc(isa)*noc(is23)                                   4d25s17
             end if                                                     4d25s17
             ixo=ionex(is)                                              4d25s17
             do i2=i2s,i2e                                              4d25s17
              i2m=i2-1                                                  4d25s17
              if(i2.eq.i2e)i1n=i1e                                      4d25s17
              do i1=i10,i1n                                             4d25s17
               i1m=i1-1                                                 4d25s17
               if(i1.le.idoub(is23))then                                4d25s17
                if(isa.eq.is23)then                                     4d25s17
                 i4=i1m
                  do i3=0,iacto(isa)-1                                   4d25s17
                   i3p=i3+idoub(isa)                                    4d25s17
                   iadx=ixo+((i3p*(i3p+1))/2)+i4                        4d25s17
                   jav=iav+i3+iacto(isa)*i2m                            4d25s17
                   bc(jav)=bc(jav)-2d0*bc(iadx)                         6d19s17
                  end do                                                4d25s17
                else                                                    4d25s17
                 i4=i1m                                                 7d14s17
                  do i3=0,iacto(isa)-1                                   4d25s17
                   i3p=i3+idoub(isa)                                    4d25s17
                   iadx=ixo+i3p+noc(isa)*i4                             4d25s17
                   jav=iav+i3+iacto(isa)*i2m                            4d25s17
                   bc(jav)=bc(jav)-2d0*bc(iadx)                         6d19s17
                  end do                                                4d25s17
                end if                                                  4d25s17
               end if                                                   4d25s17
               ixo=ixo+nrow                                             4d25s17
              end do                                                    4d25s17
              i10=1                                                     4d25s17
             end do                                                     4d25s17
            end if                                                      4d25s17
            if(isblk1(2,is).eq.isa.and.isblk1(1,is).eq.isblk1(3,is)     4d25s17
     $            .and.isblk1(4,is).eq.isb.and.igot.eq.0)then           4d25s17
             igot=1                                                     4d25s17
             is13=isblk1(1,is)                                           4d25s17
             call ilimts(noc(is13),nvirtc(isb),mynprocg,mynowprog,il,ih,4d25s17
     $            i1s,i1e,i2s,i2e)                                      4d25s17
             i10=i1s                                                    4d25s17
             i1n=noc(is13)                                              4d25s17
             nrow=noc(isa)*noc(is13)                                    4d25s17
             ixo=ionex(is)                                              4d25s17
             do i2=i2s,i2e                                              4d25s17
              i2m=i2-1                                                  4d25s17
              if(i2.eq.i2e)i1n=i1e                                      4d25s17
              do i1=i10,i1n                                             4d25s17
               i1m=i1-1                                                 4d25s17
               if(i1.le.idoub(is13))then                                4d25s17
                i4=i1m                                                  7d14s17
                 do i3=0,iacto(isa)-1                                   4d25s17
                  i3p=i3+idoub(isa)                                     4d25s17
                  iadx=ixo+i4+noc(is13)*i3p                             4d25s17
                  jav=iav+i3+iacto(isa)*i2m                             4d25s17
                  bc(jav)=bc(jav)-2d0*bc(iadx)                          11d7s17
                 end do                                                 4d25s17
               end if                                                   4d25s17
               ixo=ixo+nrow                                             4d25s17
              end do                                                    4d25s17
              i10=1                                                     4d25s17
             end do                                                     4d25s17
            end if                                                      4d25s17
            end if                                                      9d25s17
           end do                                                       4d25s17
           call dws_gsumf(bc(iav),nwds)                                 4d25s17
           nbas=nbasdws(isa)*nvirtc(isa)/nvirt(isa)                     9d25s17
           jav=iav                                                      4d25s17
           if(idoit(1).ne.0)then                                        9d25s17
           do i4=0,nvirtc(isb)-1                                        4d25s17
            i4p=i4+noc(isb)                                             6d19s17
            do i3=0,iacto(isa)-1                                        4d25s17
             i3p=i3+1+idoub(isa)                                        4d25s17
             iadh=ipt(isa)+i3p+nbas*i4p                                 6d19s17
             jav=iav+i3+iacto(isa)*i4                                   6d19s17
             bc(jav)=bc(jav)+2d0*h0mo(iadh)                             4d25s17
            end do                                                      4d25s17
           end do                                                       4d25s17
           end if                                                       9d25s17
c
c     aa nd,vd part
c
           if(idoit(2).ne.0)then                                        9d25s17
            nbas=nbasdws(isa)*(nvirtc(isa)/nvirt(isa))                  9d25s17
            do i3=0,iacto(isa)-1                                        9d25s17
             iadh=ipt(isa)+noc(isa)+nbas*(i3+idoub(isa))+1              9d25s17
             iden=iden1(isa)+iacto(isa)*i3                              9d25s17
             jav=iav                                                    9d25s17
             do i2=0,nvirtc(isa)-1                                      9d25s17
              fact=0.5d0*h0mo(iadh+i2)                                  9d25s17
              do i1=0,iacto(isa)-1                                      9d25s17
               bc(jav+i1)=bc(jav+i1)-bc(iden+i1)*fact                   9d25s17
              end do                                                    9d25s17
              jav=jav+iacto(isa)                                        9d25s17
             end do                                                     9d25s17
            end do                                                      9d25s17
           end if                                                       9d25s17
c
           i10=i1sb                                                     4d25s17
           i1n=noc(isbo)                                                4d25s17
           ihb=ihessb(isa,isb)                                          4d25s17
           do i2=i2sb,i2eb                                              4d25s17
            i2m=i2-1                                                    4d25s17
            if(i2.eq.i2eb)i1n=i1eb                                      4d25s17
            do i1=i10,i1n                                               4d25s17
             i1m=i1-1                                                   4d25s17
             if(i1.le.idoub(isbo))then                                  4d25s17
              do i4=0,iacto(isa)-1                                      4d25s17
               jav=iav+i4+iacto(isa)*i2m                                4d25s17
               jhb=ihb+i4+iacto(isa)*i1m                                6d19s17
               bc(jhb)=bc(jhb)+2d0*bc(jav)                              7d19s17
              end do                                                     4d25s17
             end if                                                     4d27s17
             ihb=ihb+nrowb                                              4d25s17
            end do                                                      4d25s17
            i10=1                                                       4d25s17
           end do                                                       4d25s17
           ibcoff=iav                                                   4d25s17
          end if                                                        4d25s17
          ibcoff=iav                                                    4d25s17
         end if                                                         4d13s17
c
c     doub&active vs active virt rotations
c
         if(isbo.eq.isa)then                                            4d25s17
          ivd=ibcoff                                                    4d27s17
          nwds=nvirtc(isb)*idoub(isao)                                  4d27s17
          if(nwds.gt.0)then                                             4d27s17
           ivd2=ivd+nwds                                                9d22s17
           ibcoff=ivd2+nwds                                             9d22s17
           call enough('buildhesscas.  8',bc,ibc)
           do i=0,nwds*2-1                                              9d22s17
            bc(ivd+i)=0d0                                               4d27s17
           end do                                                       4d27s17
           do is=1,nsdlk1                                               4d27s17
c
c     ddaa nd, vn part c and part a
c
            if(idoit(4).ne.0)then                                       9d25s17
            if(isblk1(3,is).eq.isao.and.                                9d22s17
     $           isblk1(1,is).eq.isblk1(2,is))then                      9d22s17
             is12=isblk1(1,is)                                          9d22s17
             call ilimts(noc(isao),nvirtc(isao),mynprocg,mynowprog,il,  9d22s17
     $            ih,i1s,i1e,i2s,i2e)                                   9d22s17
             i10=i1s                                                    9d22s17
             i1n=noc(isao)                                              9d22s17
             nrow=(noc(is12)*(noc(is12)+1))/2                           9d22s17
             ix1=ionex(is)-1                                            9d22s17
             do i2=i2s,i2e                                              9d22s17
              if(i2.eq.i2e)i1n=i1e                                      9d22s17
              i2m=i2-1                                                  9d22s17
              do i1=i10,i1n                                             9d22s17
               if(i1.le.idoub(isao))then                                9d22s17
                i1m=i1-1                                                 9d22s17
                jvd2=ivd2+i1m+idoub(isao)*i2m                            9d22s17
                do i34=1,idoub(is12)                                    9d22s17
                 jx1=ix1+(i34*(i34+1))/2                                9d22s17
                 bc(jvd2)=bc(jvd2)-4d0*bc(jx1)                          9d22s17
                end do                                                  9d22s17
                jvd=ivd+i2m+nvirtc(isao)*i1m                            9d22s17
                jx1=ix1+1                                               9d22s17
                do i4=0,iacto(is12)-1                                   9d22s17
                 i4p=i4+idoub(is12)                                     9d22s17
                 iden=iden1(is12)+iacto(is12)*i4                        9d22s17
                 do i3=0,iacto(is12)-1                                  9d22s17
                  i3p=i3+idoub(is12)                                    9d22s17
                  ix=max(i3p,i4p)                                       9d22s17
                  in=min(i3p,i4p)                                       9d22s17
                  iadi=jx1+((ix*(ix+1))/2)+in                           9d22s17
                  bc(jvd)=bc(jvd)+bc(iden+i3)*bc(iadi)*2d0              9d22s17
                 end do                                                 9d22s17
                end do                                                  9d22s17
               end if                                                   9d22s17
               ix1=ix1+nrow                                             9d22s17
              end do                                                    9d22s17
              i10=1                                                     9d22s17
             end do                                                     9d22s17
            end if                                                      9d22s17
c
c     ddaa nd, vn part d and part b
c
            if(isblk1(4,is).eq.isbv.and.                                9d22s17
     $           isblk1(2,is).eq.isblk1(3,is))then                      9d22s17
             is23=isblk1(2,is)                                          11d9s17
             call ilimts(noc(is23),nvirtc(isbv),mynprocg,mynowprog,il,  9d22s17
     $            ih,i1s,i1e,i2s,i2e)                                   9d22s17
             i10=i1s                                                    9d22s17
             i1n=noc(is23)                                              9d22s17
             if(isao.eq.is23)then                                       9d22s17
              nrow=(noc(is23)*(noc(is23)+1))/2                           9d22s17
             else                                                       9d22s17
              nrow=noc(is23)*noc(isao)                                  9d22s17
             end if                                                     9d22s17
             iswitch=iabs(is23-isao)                                    9d22s17
             iswitch=iswitch/max(1,iswitch)                             9d22s17
             ix1=ionex(is)                                              9d22s17
             do i2=i2s,i2e                                              9d22s17
              if(i2.eq.i2e)i1n=i1e                                      9d22s17
              i2m=i2-1                                                  9d22s17
              jvd2=ivd2+idoub(isao)*i2m                                 9d22s17
              do i1=i10,i1n                                             9d22s17
               if(i1.le.idoub(is23))then                                9d22s17
                i1m=i1-1                                                 9d22s17
                do i4=0,idoub(isao)-1                                   9d22s17
                 ix=max(i1m,i4)                                         9d22s17
                 in=min(i1m,i4)                                         9d22s17
                 ieq=((ix*(ix+1))/2)+in                                 9d22s17
                 inot=i4+noc(isao)*i1m                                  9d22s17
                 iadi=ix1+(inot-ieq)*iswitch+ieq                        9d22s17
                 bc(jvd2+i4)=bc(jvd2+i4)+2d0*bc(iadi)                   9d22s17
                end do                                                  9d22s17
               else                                                     9d22s17
                i1m=i1-1-idoub(is23)                                    9d22s17
                iden=iden1(is23)+iacto(is23)*i1m                        9d22s17
                kvd=ivd+i2m                                             9d22s17
                do i4=0,idoub(isao)-1                                   9d22s17
                 jvd=kvd+nvirtc(isbv)*i4                                9d22s17
                 do i3=0,iacto(is23)-1                                  9d22s17
                  i3p=i3+idoub(is23)                                    9d22s17
                  ix=max(i3p,i4)                                        9d22s17
                  in=min(i3p,i4)                                        9d22s17
                  ieq=((ix*(ix+1))/2)+in                                9d22s17
                  inot=i4+noc(isao)*i3p                                 9d22s17
                  iadi=ix1+(inot-ieq)*iswitch+ieq                       9d22s17
                  bc(jvd)=bc(jvd)-bc(iadi)*bc(iden+i3)                  9d22s17
                 end do                                                 9d22s17
                end do                                                  9d22s17
               end if                                                   9d22s17
               ix1=ix1+nrow                                             9d22s17
              end do                                                    9d22s17
              i10=1                                                     9d22s17
             end do                                                     9d22s17
            else if(isblk1(4,is).eq.isbv.and.                                9d22s17
     $           isblk1(1,is).eq.isblk1(3,is))then                      9d22s17
             is13=isblk1(1,is)                                          9d22s17
             call ilimts(noc(is13),nvirtc(isbv),mynprocg,mynowprog,il,  9d22s17
     $            ih,i1s,i1e,i2s,i2e)                                   9d22s17
             i10=i1s                                                    9d22s17
             i1n=noc(is13)                                              9d22s17
             nrow=noc(is13)*noc(isao)                                   9d22s17
             ix1=ionex(is)                                              9d22s17
             do i2=i2s,i2e                                              9d22s17
              if(i2.eq.i2e)i1n=i1e                                      9d22s17
              i2m=i2-1                                                  9d22s17
              jvd2=ivd2+idoub(isao)*i2m                                 9d22s17
              do i1=i10,i1n                                             9d22s17
               if(i1.le.idoub(is13))then                                9d22s17
                i1m=i1-1                                                 9d22s17
                do i4=0,idoub(isao)-1                                   9d22s17
                 inot=i1m+noc(is13)*i4                                  9d22s17
                 iadi=ix1+inot                                          9d22s17
                 bc(jvd2+i4)=bc(jvd2+i4)+2d0*bc(iadi)                   9d22s17
                end do                                                  9d22s17
               else                                                     9d22s17
                i1m=i1-1-idoub(is13)                                    9d22s17
                iden=iden1(is13)+iacto(is13)*i1m                        9d22s17
                kvd=ivd+i2m                                             9d22s17
                do i4=0,idoub(isao)-1                                   9d22s17
                 jvd=kvd+nvirtc(isbv)*i4                                9d22s17
                 do i3=0,iacto(is13)-1                                  9d22s17
                  i3p=i3+idoub(is13)                                    9d22s17
                  inot=i3p+noc(is13)*i4                                 9d22s17
                  iadi=ix1+inot                                         9d22s17
                  bc(jvd)=bc(jvd)-bc(iadi)*bc(iden+i3)                  9d22s17
                 end do                                                 9d22s17
                end do                                                  9d22s17
               end if                                                   9d22s17
               ix1=ix1+nrow                                             9d22s17
              end do                                                    9d22s17
              i10=1                                                     9d22s17
             end do                                                     9d22s17
            end if                                                      9d22s17
            end if                                                      9d25s17
c
c     dddd: "K" part
c
            if(idoit(3).ne.0)then                                       9d25s17
            if(isblk1(3,is).eq.isao.and.isblk1(4,is).eq.isb)then        4d27s17
             is12=isblk1(1,is)                                          4d27s17
             call ilimts(noc(isao),nvirtc(isb),mynprocg,mynowprog,il,ih,4d27s17
     $            i1s,i1e,i2s,i2e)                                      4d27s17
             i10=i1s                                                    4d27s17
             i1n=noc(isao)                                              4d27s17
             nrow=(noc(is12)*(noc(is12)+1))/2                           4d27s17
             iox=ionex(is)-1                                            7d18s17
             do i2=i2s,i2e                                              4d27s17
              i2m=i2-1                                                  4d27s17
              if(i2.eq.i2e)i1n=i1e                                      4d27s17
              do i1=i10,i1n                                             4d27s17
               i1m=i1-1                                                 4d27s17
               if(i1.le.idoub(isao))then                                4d27s17
                do i34=0,idoub(is12)-1                                  7d18s17
                 i34p=i34+1                                             4d27s17
                 iado=iox+(i34p*(i34p+1))/2                             4d27s17
                 jvd=ivd+i2m+nvirtc(isb)*i1m                            4d27s17
                 bc(jvd)=bc(jvd)+4d0*bc(iado)                           4d27s17
                end do                                                  4d27s17
               end if                                                   4d27s17
               iox=iox+nrow                                             4d27s17
              end do                                                    4d27s17
              i10=1                                                     4d27s17
             end do                                                     4d27s17
            end if                                                      4d27s17
c
c     "J" part
c
            igot=0
            if(isblk1(1,is).eq.isao.and.isblk1(4,is).eq.isb)then        4d27s17
             igot=1                                                     4d27s17
             is23=isblk1(2,is)                                          4d27s17
             call ilimts(noc(is23),nvirtc(isb),mynprocg,mynowprog,il,ih,4d27s17
     $            i1s,i1e,i2s,i2e)                                      4d27s17
             i10=i1s                                                    4d27s17
             i1n=noc(is23)                                              4d27s17
             if(isao.eq.is23)then                                       4d27s17
              nrow=(noc(is23)*(noc(is23)+1))/2                          4d27s17
             else                                                       4d27s17
              nrow=noc(is23)*noc(isao)                                  4d27s17
             end if                                                     4d27s17
             iox=ionex(is)                                              4d27s17
             do i2=i2s,i2e                                              4d27s17
              i2m=i2-1                                                  4d27s17
              if(i2.eq.i2e)i1n=i1e                                      4d27s17
              do i1=i10,i1n                                             4d27s17
               i1m=i1-1                                                 4d27s17
               if(i1.le.idoub(is23))then                                4d27s17
                i4=i1m                                                  4d27s17
                if(isao.eq.is23)then                                    4d27s17
                 do i3=0,idoub(isao)-1                                   4d27s17
                  ix=max(i3,i4)                                         4d27s17
                  in=min(i3,i4)                                         4d27s17
                  iado=iox+((ix*(ix+1))/2)+in                           4d27s17
                  jvd=ivd+i2m+nvirtc(isb)*i3                            4d27s17
                  bc(jvd)=bc(jvd)-2d0*bc(iado)                           4d27s17
                 end do                                                  4d27s17
                else                                                    4d27s17
                 do i3=0,idoub(isao)-1                                   4d27s17
                  iado=iox+i3+noc(isao)*i4                              4d27s17
                  jvd=ivd+i2m+nvirtc(isb)*i3                            4d27s17
                  bc(jvd)=bc(jvd)-2d0*bc(iado)                          11d7s17
                 end do                                                  4d27s17
                end if                                                  4d27s17
               end if                                                   4d27s17
               iox=iox+nrow                                             4d27s17
              end do                                                    4d27s17
              i10=1                                                     4d27s17
             end do                                                     4d27s17
            end if                                                      4d27s17
            if(isblk1(2,is).eq.isao.and.isblk1(4,is).eq.isb.and.        4d27s17
     $           igot.eq.0)then                                         4d27s17
             igot=1                                                     4d27s17
             is13=isblk1(1,is)                                          4d27s17
             call ilimts(noc(is13),nvirtc(isb),mynprocg,mynowprog,il,ih,4d27s17
     $            i1s,i1e,i2s,i2e)                                      4d27s17
             i10=i1s                                                    4d27s17
             i1n=noc(is13)                                              4d27s17
             nrow=noc(is13)*noc(isao)                                   4d27s17
             iox=ionex(is)                                              4d27s17
             do i2=i2s,i2e                                              4d27s17
              i2m=i2-1                                                  4d27s17
              if(i2.eq.i2e)i1n=i1e                                      4d27s17
              do i1=i10,i1n                                             4d27s17
               i1m=i1-1                                                 4d27s17
               if(i1.le.idoub(is13))then                                4d27s17
                i4=i1m                                                  4d27s17
                do i3=0,idoub(isao)-1                                   4d27s17
                 iado=iox+i4+noc(is13)*i3                               4d27s17
                 jvd=ivd+i2m+nvirtc(isb)*i3                             4d27s17
                 bc(jvd)=bc(jvd)-2d0*bc(iado)                           11d7s17
                end do                                                  4d27s17
               end if                                                   4d27s17
               iox=iox+nrow                                             4d27s17
              end do                                                    4d27s17
              i10=1                                                     4d27s17
             end do                                                     4d27s17
            end if                                                      4d27s17
            end if                                                      9d25s17
           end do                                                       4d27s17
           call dws_gsumf(bc(ivd),nwds*2)                               9d22s17
           nbas=nbasdws(isb)*nvirtc(isb)/nvirt(isb)                     9d25s17
           if(idoit(1).ne.0)then                                        9d25s17
           jvd=ivd                                                      4d27s17
           do i4=0,idoub(isao)-1                                        4d27s17
            do i3=0,nvirtc(isb)-1                                       4d27s17
             i3p=i3+noc(isb)+1                                          4d27s17
             iadh=ipt(isb)+i3p+nbas*i4                                  4d27s17
             bc(jvd)=bc(jvd)+2d0*h0mo(iadh)                             4d27s17
             jvd=jvd+1                                                  4d27s17
            end do                                                      4d27s17
           end do                                                       4d27s17
           end if                                                       9d25s17
c
           i10=i1sb                                                     4d27s17
           i1n=noc(isbo)                                                4d27s17
           ihb=ihessb(isa,isb)                                          4d27s17
           do i2=i2sb,i2eb                                              4d27s17
            i2m=i2-1                                                    4d27s17
            if(i2.eq.i2eb)i1n=i1eb                                      4d27s17
            do i1=i10,i1n                                               4d27s17
             i1m=i1-1-idoub(isbo)                                       4d27s17
             if(i1.gt.idoub(isbo))then                                  5d12s22
              do i4=0,idoub(isao)-1                                     4d27s17
               jvd=ivd+i2m+nvirtc(isb)*i4                               4d27s17
               jhb=ihb+i1m+iacto(isbo)*i4                               4d27s17
               bc(jhb)=bc(jhb)+bc(jvd)                                  4d27s17
              end do                                                    4d27s17
             end if                                                     4d27s17
             ihb=ihb+nrowb                                              4d27s17
            end do                                                      4d27s17
            i10=1                                                       4d27s17
           end do                                                       4d27s17
c
c     aa nd,vn part
c
           if(idoit(2).ne.0)then                                        9d25s17
           jvd2=ivd2                                                    9d25s17
           nbas=nbasdws(isao)*nvirtc(isao)/nvirt(isao)                   9d25s17
           iad=ipt(isao)+nbas*noc(isao)+1                               9d25s17
           do i1=0,nvirtc(isao)-1                                       9d25s17
            do i2=0,idoub(isao)-1                                        9d25s17
             bc(jvd2+i2)=bc(jvd2+i2)-2d0*h0mo(iad+i2)                   9d25s17
            end do                                                      9d25s17
            jvd2=jvd2+idoub(isao)                                       9d25s17
            iad=iad+nbas                                                9d25s17
           end do                                                       9d25s17
           end if                                                       9d25s17
c
c     ddaa nd,vn part c&d
c
           i10=i1sb                                                     9d22s17
           i1n=noc(isbo)                                                9d22s17
           ihb=ihessb(isa,isb)                                          9d22s17
           do i2=i2sb,i2eb                                              9d22s17
            if(i2.eq.i2eb)i1n=i1eb                                      9d22s17
            i2m=i2-1                                                    9d22s17
            do i1=i10,i1n                                               9d22s17
             if(i1.gt.idoub(isbo))then                                  9d22s17
              i1m=i1-1-idoub(isbo)                                      9d22s17
              iden=iden1(isbo)+iacto(isbo)*i1m                          9d22s17
              jhb=ihb                                                   9d22s17
              do i4=0,idoub(isao)-1                                     9d22s17
               jvd2=ivd2+i4+idoub(isao)*i2m                             9d22s17
               do i3=0,iacto(isav)-1                                    9d22s17
                bc(jhb)=bc(jhb)+bc(iden+i3)*bc(jvd2)                    9d22s17
                jhb=jhb+1                                               9d22s17
               end do                                                   9d22s17
              end do                                                    9d22s17
             end if                                                     9d22s17
             ihb=ihb+nrowb                                              9d22s17
            end do                                                      9d22s17
            i10=1                                                       9d22s17
           end do                                                       9d22s17
           ibcoff=ivd                                                   4d27s17
          end if                                                        4d27s17
         end if                                                         4d25s17
         if(nrowa*ncola.gt.0.and.idwsdeb.gt.10.and.lsym)then            11d30s17
          write(6,*)('doub-active part of hess (eshift)'),mynowprog                4d25s17
          write(6,*)('symmetry type '),isbo,isao,isa,isb,ihessa(isa,isb)
          call prntm2(bc(ihessa(isa,isb)),nrowa,ncola,nrowa)             4d25s17
         end if
         if(nrowb*nhereb.gt.0.and.idwsdeb.gt.0)then
          rms=0d0                                                       11d9s17
          do i=0,nrowb*nhereb-1                                         11d9s17
           rms=rms+bc(ihessb(isa,isb)+i)**2                             11d9s17
          end do                                                        11d9s17
          rms=sqrt(rms/dfloat(nrowb*nhereb))                            11d9s17
          if(rms.gt.1d-5)then                                           11d9s17
          write(6,*)('doub-active and all-virt part of hess (eshift)')   4d27s17
          write(6,*)('isa,isao,isbo,isb '),isa,isao,isbo,isb,
     $         ihessb(isa,isb)
          call prntm2(bc(ihessb(isa,isb)),nrowb,nhereb,nrowb)             4d25s17
          end if                                                        11d9s17
         end if
         if(nsymb.eq.1.and.idwsdeb.gt.10)then
          write(6,*)('hessa: ')
          if(lsym)call printa(bc(ihessa(1,1)),idoub,0,idoub,0,iacto,    11d30s17
     $         idoub,                                                   11d30s17
     $         iacto,idoub,bc(ibcoff))
          write(6,*)('hessb: ')
          ihb=ibcoff
          nwds=nrowb*noc(1)*nvirtc(1)
          ibcoff=ihb+nwds
          call enough('buildhesscas.  9',bc,ibc)
          do i1234=0,nwds-1
           bc(ihb+i1234)=0d0
          end do
          do i12=0,nhereb-1
           i12p=i12+ilb-1                                               4d11s17
           iad1=ihb+nrowb*i12p                                          4d11s17
           iad2=ihessb(isa,isb)+nrowb*i12                                4d11s17
           do i34=0,nrowb-1
            bc(iad1+i34)=bc(iad2+i34)                                   4d11s17
           end do
          end do
          call dws_gsumf(bc(ihb),nwds)
          write(6,*)('hessb over all processors ')
          call prntm2(bc(ihb),nrowb,noc(1)*nvirtc(1),nrowb)
          call printa(bc(ihb),iacto,idoub,idoub,0,noc,0,nvirtc,noc,
     $         bc(ibcoff))
          ibcoff=ihb
         end if
         if(idoit(4).ne.0.and.lsym)then                                 11d30s17
           if(isbv.eq.isav)then                                         10d18s17
            prefact=1d0                                                 10d18s17
           else                                                         10d18s17
            prefact=0.5d0                                               10d18s17
           end if                                                       10d18s17
           do is=1,nsdlkk                                               9d28s17
c
c     ddaa vd, vn parts e&g
c     these terms are doubled, because due to distribution of kmats,
c     we can not get the proper symmetric result. Thus we will average
c     the off diagonal matrix elements.
c
c     there appears to be a problem when isa ne isb
c     rows: isao,isbo  cols: isbv,isav
c
            if(isblkk(1,is).eq.isao.and.isblkk(3,is).eq.isbv.and.       9d28s17
     $           isblkk(4,is).eq.isav)then                                9d28s17
             i10=i1sh                                                   9d28s17
             i1n=nvirtc(isbv)                                           9d28s17
             ikm=kmats(is)                                              9d28s17
             iadh=ihessc(isa,isb)                                       9d28s17
             do i2=i2sh,i2eh                                            9d28s17
              if(i2.eq.i2eh)i1n=i1eh                                    9d28s17
              do i1=i10,i1n                                             9d28s17
               if(i1.ge.i2.or.isbv.ne.isav)then                         10d2s17
                fmul=16d0*prefact                                       10d18s17
                if(i1.eq.i2.and.isbv.eq.isav)fmul=8d0                   11d27s17
                do i4=0,iacto(isbo)-1                                    9d28s17
                 i4p=i4+idoub(isbo)                                      9d28s17
                 iden=iden1(isbo)+iacto(isbo)*i4                         9d28s17
                 iadi=ikm+noc(isao)*i4p                                  9d28s17
                 do i5=0,iacto(isbo)-1                                   9d28s17
                  fact=fmul*bc(iden+i5)                                  10d10s17
                  jadh=iadh+noc(isao)*(i5+idoub(isbo))                   9d28s17
                  do i3=0,idoub(isao)-1                                  9d28s17
                   bc(jadh+i3)=bc(jadh+i3)+fact*bc(iadi+i3)              9d28s17
                  end do                                                 9d28s17
                 end do                                                  9d28s17
                end do                                                   9d28s17
               end if                                                   10d2s17
               ikm=ikm+nrowh                                            9d28s17
               iadh=iadh+nrowh                                          9d28s17
              end do                                                    9d28s17
              i10=1                                                     9d28s17
             end do                                                     9d28s17
            end if                                                      9d28s17
            if(isblkk(2,is).eq.isao.and.isblkk(3,is).eq.isbv.and.       10d18s17
     $           isblkk(4,is).eq.isav)then                                9d28s17
             i10=i1sh                                                   9d28s17
             i1n=nvirtc(isbv)                                           9d28s17
             ikm=kmats(is)                                              9d28s17
             iadh=ihessc(isa,isb)                                       9d28s17
             do i2=i2sh,i2eh                                            9d28s17
              if(i2.eq.i2eh)i1n=i1eh                                    9d28s17
              do i1=i10,i1n                                             9d28s17
               if(i1.ge.i2.or.isbv.ne.isav)then                         10d2s17
                fmul=-4d0*prefact                                       10d18s17
                if(i1.eq.i2.and.isbv.eq.isav)fmul=-2d0                  11d27s17
                do i4=0,idoub(isao)-1                                   10d18s17
                 jadh=iadh+i4+noc(isao)*idoub(isbo)                     10d18s17
                 iadi=ikm+idoub(isbo)+noc(isbo)*i4                      10d18s17
                 do i3=0,iacto(isbo)-1                                  10d18s17
                  fact=fmul*bc(iadi+i3)                                 10d18s17
                  iden=iden1(isbo)+iacto(isbo)*i3                       10d18s17
                  do i5=0,iacto(isbo)-1                                 10d18s17
                   bc(jadh+i5*noc(isao))=bc(jadh+i5*noc(isao))          10d18s17
     $                  +fact*bc(iden+i5)                               10d18s17
                  end do                                                 9d28s17
                 end do                                                  9d28s17
                end do                                                   9d28s17
               end if                                                   10d2s17
               ikm=ikm+nrowh                                            9d28s17
               iadh=iadh+nrowh                                          9d28s17
              end do                                                    9d28s17
              i10=1                                                     9d28s17
             end do                                                     9d28s17
            end if                                                      9d28s17
c
c     ddaa vn, vd part e
c
            if(isblkk(1,is).eq.isao.and.isblkk(3,is).eq.isbv.and.       9d28s17
     $           isblkk(4,is).eq.isav)then                                9d28s17
             i10=i1sh                                                   9d28s17
             i1n=nvirtc(isbv)                                           9d28s17
             ikm=kmats(is)                                              9d28s17
             iadh=ihessc(isa,isb)                                       9d28s17
             do i2=i2sh,i2eh                                            9d28s17
              if(i2.eq.i2eh)i1n=i1eh                                    9d28s17
              do i1=i10,i1n                                             9d28s17
               if(i1.ge.i2.or.isbv.ne.isav)then                         10d2s17
                fmul=16d0*prefact                                       10d18s17
                if(i1.eq.i2.and.isbv.eq.isav)fmul=8d0                   11d27s17
                do i4=0,idoub(isbo)-1                                    9d28s17
                 iadi=ikm+idoub(isao)+noc(isao)*i4                       10d10s17
                 do i5=0,iacto(isao)-1                                   9d28s17
                  iden=iden1(isao)+iacto(isao)*i5                        10d10s17
                  jadh=iadh+idoub(isao)+i5+noc(isao)*i4                  10d10s17
                  do i3=0,iacto(isao)-1                                  10d10s17
                   bc(jadh)=bc(jadh)+fmul*bc(iden+i3)*bc(iadi+i3)        10d10s17
                  end do                                                 9d28s17
                 end do                                                  9d28s17
                end do                                                   9d28s17
               end if                                                   10d2s17
               ikm=ikm+nrowh                                            9d28s17
               iadh=iadh+nrowh                                          9d28s17
              end do                                                    9d28s17
              i10=1                                                     9d28s17
             end do                                                     9d28s17
            end if                                                      9d28s17
c
c     ddaa vn, vd part g, well, sort of. Since this involves transpose
c     of K, distribution causes us to actually do ddaa vd, vn part g
c
            if(isblkk(2,is).eq.isao.and.isblkk(3,is).eq.isbv.and.       10d18s17
     $           isblkk(4,is).eq.isav)then                                9d28s17
             i10=i1sh                                                   9d28s17
             i1n=nvirtc(isbv)                                           9d28s17
             ikm=kmats(is)                                              9d28s17
             iadh=ihessc(isa,isb)                                       9d28s17
             do i2=i2sh,i2eh                                            9d28s17
              if(i2.eq.i2eh)i1n=i1eh                                    9d28s17
              do i1=i10,i1n                                             9d28s17
               if(i1.ge.i2.or.isbv.ne.isav)then                         10d2s17
                fmul=-4d0*prefact                                       10d18s17
                if(i1.eq.i2.and.isbv.eq.isav)fmul=-2d0                  11d27s17
                do i4=0,iacto(isao)-1                                   10d18s17
                 jadh=iadh+i4+idoub(isao)                               10d18s17
                 iden=iden1(isao)+iacto(isao)*i4                        10d18s17
                 do i5=0,idoub(isbo)-1                                  10d18s17
                  iadi=ikm+i5+noc(isbo)*idoub(isao)                     10d18s17
                  kadh=jadh+noc(isao)*i5                                10d18s17
                  do i3=0,iacto(isao)-1                                 10d18s17
                   bc(kadh)=bc(kadh)+fmul*bc(iadi+i3*noc(isbo))         10d18s17
     $                  *bc(iden+i3)                                    10d18s17
                  end do                                                 9d28s17
                 end do                                                  9d28s17
                end do                                                   9d28s17
               end if                                                   10d2s17
               ikm=ikm+nrowh                                            9d28s17
               iadh=iadh+nrowh                                          9d28s17
              end do                                                    9d28s17
              i10=1                                                     9d28s17
             end do                                                     9d28s17
            end if                                                      9d28s17
           end do                                                       9d28s17
          do is=1,nsdlk                                                 9d28s17
c
c     ddaa vd, vn part f
c
           if(isblk(3,is).eq.isbv.and.isblk(4,is).eq.isav.and.          9d28s17
     $          isblk(2,is).eq.isao)then                                  9d28s17
            iswitch=iabs(isao-isbo)                                     9d28s17
            iswitch=iswitch/max(1,iswitch)                              9d28s17
            i10=i1sh                                                    9d28s17
            i1n=nvirtc(isbv)                                            9d28s17
            ijm=jmats(is)                                               9d28s17
            iadh=ihessc(isa,isb)                                        9d29s17
            if(isao.eq.isbo)then                                        9d29s17
             nrow=(noc(isao)*(noc(isao)+1))/2                           9d29s17
            else                                                        9d29s17
             nrow=noc(isao)*noc(isbo)                                   9d29s17
            end if                                                      9d29s17
            do i2=i2sh,i2eh                                             9d28s17
             if(i2.eq.i2eh)i1n=i1eh                                     9d28s17
             do i1=i10,i1n                                              9d28s17
              do i4=0,idoub(isao)-1                                     9d29s17
               do i3=0,iacto(isbo)-1                                    9d29s17
                i3p=i3+idoub(isbo)
                ix=max(i4,i3p)                                          9d29s17
                in=min(i4,i3p)                                          9d29s17
                ieq=((ix*(ix+1))/2)+in                                  9d29s17
                inot=i3p+noc(isbo)*i4                                   9d29s17
                iadi=ijm+(inot-ieq)*iswitch+ieq                         9d29s17
                fact=-2d0*bc(iadi)                                      10d18s17
                iden=iden1(isbo)+iacto(isbo)*i3                         9d29s17
                do i5=0,iacto(isbo)-1                                   9d29s17
                 i5p=i5+idoub(isbo)                                     9d29s17
                 jadh=iadh+i4+noc(isao)*i5p                             9d29s17
                 bc(jadh)=bc(jadh)+fact*bc(iden+i5)                     9d29s17
                end do                                                  9d29s17
               end do                                                   9d29s17
              end do                                                    9d29s17
              ijm=ijm+nrow
              iadh=iadh+nrowh                                           9d29s17
             end do                                                     9d29s17
             i10=1                                                      9d29s17
            end do                                                      9d28s17
           else if(isblk(3,is).eq.isbv.and.isblk(4,is).eq.isav.and.     9d29s17
     $        isblk(1,is).eq.isao)then                                  9d28s17
            i10=i1sh                                                    9d28s17
            i1n=nvirtc(isbv)                                            9d28s17
            ijm=jmats(is)                                               9d28s17
            iadh=ihessc(isa,isb)                                        9d29s17
            nrow=noc(isao)*noc(isbo)                                    9d29s17
            do i2=i2sh,i2eh                                             9d28s17
             if(i2.eq.i2eh)i1n=i1eh                                     9d28s17
             do i1=i10,i1n                                              9d28s17
              do i4=0,idoub(isao)-1                                     9d29s17
               do i3=0,iacto(isbo)-1                                    9d29s17
                i3p=i3+idoub(isbo)
                inot=i4+noc(isao)*i3p                                   9d29s17
                iadi=ijm+inot                                           9d29s17
                fact=-2d0*bc(iadi)                                      10d18s17
                iden=iden1(isbo)+iacto(isbo)*i3                         9d29s17
                do i5=0,iacto(isbo)-1                                   9d29s17
                 i5p=i5+idoub(isbo)                                     9d29s17
                 jadh=iadh+i4+noc(isao)*i5p                             9d29s17
                 bc(jadh)=bc(jadh)+fact*bc(iden+i5)                     9d29s17
                end do                                                  9d29s17
               end do                                                   9d29s17
              end do                                                    9d29s17
              ijm=ijm+nrow
              iadh=iadh+nrowh                                           9d29s17
             end do                                                     9d29s17
             i10=1                                                      9d29s17
            end do                                                      9d28s17
           end if                                                       9d28s17
c
c     ddaa vn, vd part f
c
           if(isblk(3,is).eq.isbv.and.isblk(4,is).eq.isav.and.          9d28s17
     $        isblk(2,is).eq.isao)then                                  9d28s17
            iswitch=iabs(isao-isbo)                                     9d28s17
            iswitch=iswitch/max(1,iswitch)                              9d28s17
            i10=i1sh                                                    9d28s17
            i1n=nvirtc(isbv)                                            9d28s17
            ijm=jmats(is)                                               9d28s17
            iadh=ihessc(isa,isb)                                        9d29s17
            if(isao.eq.isbo)then                                        9d29s17
             nrow=(noc(isao)*(noc(isao)+1))/2                           9d29s17
             fmul=-2d0                                                  12d26s17
            else                                                        9d29s17
             nrow=noc(isao)*noc(isbo)                                   9d29s17
             fmul=-2d0                                                  12d26s17
            end if                                                      9d29s17
            do i2=i2sh,i2eh                                             9d28s17
             if(i2.eq.i2eh)i1n=i1eh                                     9d28s17
             do i1=i10,i1n                                              9d28s17
              do i4=0,iacto(isao)-1                                     10d10s17
               i4p=i4+idoub(isao)                                       10d10s17
               do i3=0,idoub(isbo)-1                                    9d29s17
                ix=max(i4p,i3)                                          10d10s17
                in=min(i4p,i3)                                          10d10s17
                ieq=((ix*(ix+1))/2)+in                                  9d29s17
                inot=i3+noc(isbo)*i4p                                   10d10s17
                iadi=ijm+(inot-ieq)*iswitch+ieq                         9d29s17
                fact=fmul*bc(iadi)                                      12d26s17
                iden=iden1(isao)+iacto(isao)*i4                         9d29s17
                do i5=0,iacto(isao)-1                                   9d29s17
                 i5p=i5+idoub(isao)                                     9d29s17
                 jadh=iadh+i5p+noc(isao)*i3                             10d10s17
                 bc(jadh)=bc(jadh)+fact*bc(iden+i5)                     9d29s17
                end do                                                  9d29s17
               end do                                                   9d29s17
              end do                                                    9d29s17
              ijm=ijm+nrow
              iadh=iadh+nrowh                                           9d29s17
             end do                                                     9d29s17
             i10=1                                                      9d29s17
            end do                                                      9d28s17
           else if(isblk(3,is).eq.isbv.and.isblk(4,is).eq.isav.and.     9d29s17
     $        isblk(1,is).eq.isao)then                                  9d28s17
            i10=i1sh                                                    9d28s17
            i1n=nvirtc(isbv)                                            9d28s17
            ijm=jmats(is)                                               9d28s17
            iadh=ihessc(isa,isb)                                        9d29s17
            nrow=noc(isao)*noc(isbo)                                    9d29s17
            do i2=i2sh,i2eh                                             9d28s17
             if(i2.eq.i2eh)i1n=i1eh                                     9d28s17
             do i1=i10,i1n                                              9d28s17
              do i4=0,iacto(isao)-1                                     10d10s17
               i4p=i4+idoub(isao)                                       10d10s17
               do i3=0,idoub(isbo)-1                                    10d10s17
                inot=i4p+noc(isao)*i3                                   10d10s17
                iadi=ijm+inot                                           9d29s17
                fact=-2d0*bc(iadi)                                      10d18s17
                iden=iden1(isao)+iacto(isao)*i4                         9d29s17
                do i5=0,iacto(isao)-1                                   9d29s17
                 i5p=i5+idoub(isao)                                     9d29s17
                 jadh=iadh+i5p+noc(isao)*i3                             10d10s17
                 bc(jadh)=bc(jadh)+fact*bc(iden+i5)                     9d29s17
                end do                                                  9d29s17
               end do                                                   9d29s17
              end do                                                    9d29s17
              ijm=ijm+nrow
              iadh=iadh+nrowh                                           9d29s17
             end do                                                     9d29s17
             i10=1                                                      9d29s17
            end do                                                      9d28s17
           end if                                                       9d28s17
          end do                                                        9d28s17
         end if                                                         9d28s17
c
c     doub-virt rotations                                               4d27s17
c
         if(isa.eq.isb.and.idoub(isao).ne.0)then                        4d27s17
          noo=(idoub(isao)*(idoub(isao)+1))/2
          noox=(noc(isao)*(noc(isao)+1))/2                              7d19s17
          nvv=(nvirtc(isav)*(nvirtc(isav)+1))/2                         3d29s16
          nda=idoub(isao)*iacto(isao)                                   9d27s17
          ioop=ibcoff
          ivv=ioop+noo
          ida=ivv+nvv                                                   9d27s17
          ibcoff=ida+nda                                                9d27s17
           kg1=ida
          call enough('buildhesscas. 10',bc,ibc)
          do i=0,noo+nvv+nda-1                                          9d27s17
           bc(ioop+i)=0d0
          end do
          if(idoit(4).ne.0)then                                         9d28s17
c
c     ddaa contribution
c     vd,vd
c
           do is=1,nsdlk                                                9d28s17
c
c     j part
c
            if(isblk(1,is).eq.isao.and.isblk(2,is).eq.isao)then         9d28s17
             is34=isblk(3,is)                                           9d28s17
             call ilimts(noc(is34),noc(is34),mynprocg,mynowprog,il,     9d27s17
     $            ih,i1s,i1e,i2s,i2e)                                    9d27s17
             i10=i1s                                                    9d27s17
             i1n=noc(is34)                                              9d27s17
             io=ioooo(is)                                               9d27s17
             nrow=(noc(isao)*(noc(isao)+1))/2                           9d27s17
             do i2=i2s,i2e                                              9d27s17
              if(i2.eq.i2e)i1n=i1e                                      9d27s17
              i2m=i2-1-idoub(is34)                                      9d28s17
              do i1=i10,i1n                                             9d28s17
               if(i2.gt.idoub(is34).and.i1.gt.idoub(is34))then          9d28s17
                i1m=i1-1-idoub(is34)                                    9d28s17
                iden=iden1(is34)+i1m+iacto(is34)*i2m                    9d28s17
                fact=-bc(iden)*4d0                                      9d28s17
                joop=ioop                                               9d28s17
                jo=io                                                   9d28s17
                do i4=0,idoub(isao)-1                                   9d28s17
                 do i3=0,i4                                             9d28s17
                  bc(joop+i3)=bc(joop+i3)+fact*bc(jo+i3)                9d28s17
                 end do                                                 9d28s17
                 joop=joop+i4+1                                         9d28s17
                 jo=jo+i4+1                                             9d28s17
                end do                                                  9d28s17
               end if                                                   9d28s17
               io=io+nrow                                               9d28s17
              end do                                                    9d28s17
              i10=1                                                     9d28s17
             end do                                                     9d28s17
            end if                                                      9d28s17
c
c     k part
c
            if(isblk(2,is).eq.isao.and.isblk(4,is).eq.isao)then         9d28s17
             is13=isblk(1,is)                                           9d28s17
             iswitch=iabs(is13-isao)                                     9d28s17
             iswitch=iswitch/max(1,iswitch)                             9d28s17
             call ilimts(noc(is13),noc(isao),mynprocg,mynowprog,il,     9d27s17
     $            ih,i1s,i1e,i2s,i2e)                                    9d27s17
             i10=i1s                                                    9d27s17
             i1n=noc(is13)                                              9d27s17
             io=ioooo(is)                                               9d27s17
             if(is13.eq.isao)then                                       9d28s17
              nrow=(noc(isao)*(noc(isao)+1))/2                           9d27s17
             else                                                       9d28s17
              nrow=noc(isao)*noc(is13)                                  9d28s17
             end if                                                     9d28s17
             do i2=i2s,i2e                                              9d27s17
              if(i2.eq.i2e)i1n=i1e                                      9d28s17
              i2m=i2-1                                                  9d28s17
              joop=ioop+((i2m*(i2m+1))/2)                               9d28s17
              do i1=i10,i1n                                             9d28s17
               if(i1.gt.idoub(is13).and.i2.le.idoub(isao))then          9d28s17
                i1m=i1-1-idoub(is13)                                    9d28s17
                iden=iden1(is13)+iacto(is13)*i1m                        9d28s17
                do i4=0,i2m                                             9d28s17
                 do i3=0,iacto(is13)-1                                  9d28s17
                  i3p=i3+idoub(is13)                                    9d28s17
                  ix=max(i3p,i4)                                        9d28s17
                  in=min(i3p,i4)                                        9d28s17
                  ieq=((ix*(ix+1))/2)+in                                9d28s17
                  inot=i3p+noc(is13)*i4                                 9d28s17
                  iadi=io+(inot-ieq)*iswitch+ieq                        9d28s17
                  bc(joop+i4)=bc(joop+i4)+2d0*bc(iden+i3)*bc(iadi)      9d28s17
                 end do                                                 9d28s17
                end do                                                  9d28s17
               end if                                                   9d28s17
               io=io+nrow                                               9d28s17
              end do                                                    9d28s17
              i10=1                                                     9d28s17
             end do                                                     9d28s17
            else if(isblk(2,is).eq.isao.and.isblk(3,is).eq.isao)then    9d28s17
             is14=isblk(1,is)                                           9d28s17
             call ilimts(noc(isao),noc(is14),mynprocg,mynowprog,il,     9d27s17
     $            ih,i1s,i1e,i2s,i2e)                                    9d27s17
             i10=i1s                                                    9d27s17
             i1n=noc(isao)                                              9d27s17
             io=ioooo(is)                                               9d27s17
             nrow=noc(isao)*noc(is14)                                   9d28s17
             do i2=i2s,i2e                                              9d27s17
              if(i2.eq.i2e)i1n=i1e                                      9d28s17
              i2m=i2-1-idoub(is14)                                      9d28s17
              iden=iden1(is14)+iacto(is14)*i2m                          9d28s17
              do i1=i10,i1n                                             9d28s17
               if(i2.gt.idoub(is14).and.i1.le.idoub(isao))then          9d28s17
                i1m=i1-1                                                9d28s17
                joop=ioop+((i1m*(i1m+1))/2)                             9d28s17
                do i4=0,i1m                                             9d28s17
                 do i3=0,iacto(is14)-1                                  9d28s17
                  i3p=i3+idoub(is14)                                    9d28s17
                  iadi=io+i3p+noc(is14)*i4                              9d28s17
                  bc(joop+i4)=bc(joop+i4)+2d0*bc(iden+i3)*bc(iadi)      9d28s17
                 end do                                                 9d28s17
                end do                                                  9d28s17
               end if                                                   9d28s17
               io=io+nrow                                               9d28s17
              end do                                                    9d28s17
              i10=1                                                     9d28s17
             end do                                                     9d28s17
            else if(isblk(1,is).eq.isao.and.isblk(3,is).eq.isao)then    9d28s17
             is24=isblk(2,is)                                           9d28s17
             call ilimts(noc(isao),noc(is24),mynprocg,mynowprog,il,     9d27s17
     $            ih,i1s,i1e,i2s,i2e)                                    9d27s17
             i10=i1s                                                    9d27s17
             i1n=noc(isao)                                              9d27s17
             io=ioooo(is)                                               9d27s17
             nrow=noc(isao)*noc(is24)                                   9d28s17
             do i2=i2s,i2e                                              9d27s17
              if(i2.eq.i2e)i1n=i1e                                      9d28s17
              i2m=i2-1-idoub(is24)                                      9d28s17
              iden=iden1(is24)+iacto(is24)*i2m                          9d28s17
              do i1=i10,i1n                                             9d28s17
               if(i2.gt.idoub(is24).and.i1.le.idoub(isao))then          9d28s17
                i1m=i1-1                                                9d28s17
                joop=ioop+((i1m*(i1m+1))/2)                               9d28s17
                do i4=0,i1m                                             9d28s17
                 do i3=0,iacto(is24)-1                                  9d28s17
                  i3p=i3+idoub(is24)                                    9d28s17
                  iadi=io+i4+noc(isao)*i3p                              9d28s17
                  bc(joop+i4)=bc(joop+i4)+2d0*bc(iden+i3)*bc(iadi)      9d28s17
                 end do                                                 9d28s17
                end do                                                  9d28s17
               end if                                                   9d28s17
               io=io+nrow                                               9d28s17
              end do                                                    9d28s17
              i10=1                                                     9d28s17
             end do                                                     9d28s17
            else if(isblk(1,is).eq.isao.and.isblk(4,is).eq.isao)then    9d28s17
             is23=isblk(2,is)                                           9d28s17
             call ilimts(noc(is23),noc(isao),mynprocg,mynowprog,il,     9d27s17
     $            ih,i1s,i1e,i2s,i2e)                                    9d27s17
             i10=i1s                                                    9d27s17
             i1n=noc(is23)                                              9d27s17
             io=ioooo(is)                                               9d27s17
             nrow=noc(isao)*noc(is23)                                   9d28s17
             do i2=i2s,i2e                                              9d27s17
              if(i2.eq.i2e)i1n=i1e                                      9d28s17
              i2m=i2-1                                                  9d28s17
              joop=ioop+((i2m*(i2m+1))/2)                               9d28s17
              do i1=i10,i1n                                             9d28s17
               if(i2.le.idoub(isao).and.i1.gt.idoub(is23))then          9d28s17
                i1m=i1-1-idoub(is23)                                    9d28s17
                iden=iden1(is23)+iacto(is23)*i1m                        9d28s17
                do i4=0,i2m                                             9d28s17
                 do i3=0,iacto(is23)-1                                  9d28s17
                  i3p=i3+idoub(is23)                                    9d28s17
                  iadi=io+i4+noc(isao)*i3p                              9d28s17
                  bc(joop+i4)=bc(joop+i4)+2d0*bc(iden+i3)*bc(iadi)      9d28s17
                 end do                                                 9d28s17
                end do                                                  9d28s17
               end if                                                   9d28s17
               io=io+nrow                                               9d28s17
              end do                                                    9d28s17
              i10=1                                                     9d28s17
             end do                                                     9d28s17
            end if                                                      9d28s17
           end do                                                       9d28s17
c
c     now terms depending on 2 virtual indecies
c
c
c     jmats part
c
           do is=1,nsdlk                                                9d28s17
            if(isblk(3,is).eq.isav.and.isblk(4,is).eq.isav)then         9d28s17
             is12=isblk(1,is)                                           9d28s17
             call ilimts(nvirtc(isav),nvirtc(isav),mynprocg,mynowprog,  9d28s17
     $            il,ih,i1s,i1e,i2s,i2e)                                9d28s17
             i10=i1s                                                    9d27s17
             i1n=nvirtc(isav)                                           9d28s17
             ijm=jmats(is)                                              9d28s17
             nrow=(noc(is12)*(noc(is12)+1))/2                           9d27s17
             do i2=i2s,i2e                                              9d28s17
              if(i2.eq.i2e)i1n=i1e                                      9d28s17
              i2m=i2-1                                                  9d28s17
              do i1=i10,i1n                                             9d28s17
               if(i1.ge.i2)then                                         9d28s17
                jvv=ivv+((i1*(i1-1))/2)+i2m                             9d28s17
                do i4=0,iacto(is12)-1                                   9d28s17
                 i4p=i4+idoub(is12)                                     9d28s17
                 iden=iden1(is12)+iacto(is12)*i4                        9d28s17
                 do i3=0,iacto(is12)-1                                  9d28s17
                  i3p=i3+idoub(is12)                                    9d28s17
                  ix=max(i3p,i4p)                                       9d28s17
                  in=min(i3p,i4p)                                       9d28s17
                  iadi=ijm+((ix*(ix+1))/2)+in                           9d28s17
                  bc(jvv)=bc(jvv)+4d0*bc(iadi)*bc(iden+i3)              9d28s17
                 end do                                                 9d28s17
                end do                                                  9d28s17
               end if                                                   9d28s17
               ijm=ijm+nrow                                             9d28s17
              end do                                                    9d28s17
              i10=1                                                     9d28s17
             end do                                                     9d28s17
            end if                                                      9d28s17
           end do                                                       9d28s17
           do is=1,nsdlkk                                               9d28s17
c
c     kmats part
c
            if(isblkk(3,is).eq.isav.and.isblkk(4,is).eq.isav)then       9d28s17
             is12=isblkk(1,is)                                          9d28s17
             call ilimts(nvirtc(isav),nvirtc(isav),mynprocg,mynowprog,  9d28s17
     $            il,ih,i1s,i1e,i2s,i2e)                                9d28s17
             i10=i1s                                                    9d27s17
             i1n=nvirtc(isav)                                           9d28s17
             ikm=kmats(is)                                              9d28s17
             nrow=noc(is12)*noc(is12)                                   9d28s17
             do i2=i2s,i2e                                              9d28s17
              i2m=i2-1                                                  9d28s17
              if(i2.eq.i2e)i1n=i1e                                      10d2s17
              do i1=i10,i1n                                             9d28s17
               if(i1.ge.i2)then                                         9d28s17
                jvv=ivv+((i1*(i1-1))/2)+i2m                             9d28s17
                do i4=0,iacto(is12)-1                                   9d28s17
                 i4p=i4+idoub(is12)                                     9d28s17
                 iden=iden1(is12)+iacto(is12)*i4                        9d28s17
                 iadi=ikm+idoub(is12)+noc(is12)*i4p                     9d28s17
                 do i3=0,iacto(is12)-1                                  9d28s17
                  bc(jvv)=bc(jvv)-2d0*bc(iadi+i3)*bc(iden+i3)           9d28s17
                 end do                                                 9d28s17
                end do                                                  9d28s17
               end if                                                   9d28s17
               ikm=ikm+nrow                                             9d28s17
              end do                                                    9d28s17
              i10=1                                                     9d28s17
             end do                                                     9d28s17
            end if                                                      9d28s17
           end do                                                       9d28s17
c
c     ddaa contribution to vd,va
c
           do is=1,nsdlk                                                9d28s17
c
c     part a
c
            if(isblk(1,is).eq.isao.and.isblk(2,is).eq.isao)then         9d28s17
             is34=isblk(3,is)                                           9d28s17
             call ilimts(noc(is34),noc(is34),mynprocg,mynowprog,il,     9d27s17
     $            ih,i1s,i1e,i2s,i2e)                                    9d27s17
             i10=i1s                                                    9d27s17
             i1n=noc(is34)                                              9d27s17
             io=ioooo(is)                                               9d27s17
             nrow=(noc(isao)*(noc(isao)+1))/2                           9d27s17
             do i2=i2s,i2e                                              9d27s17
              if(i2.eq.i2e)i1n=i1e                                      9d27s17
              i2m=i2-1-idoub(is34)                                      9d28s17
              do i1=i10,i1n                                             9d28s17
               if(i1.gt.idoub(is34).and.i2.gt.idoub(is34))then          9d28s17
                i1m=i1-1-idoub(is34)                                    9d28s17
                iden=iden1(is34)+i1m+iacto(is34)*i2m                    9d28s17
                fact=-2d0*bc(iden)                                      9d28s17
                jda=ida                                                 9d28s17
                do i4=0,idoub(isao)-1                                   9d28s17
                 do i3=0,iacto(isao)-1                                   9d29s17
                  i3p=i3+idoub(isao)                                    9d28s17
                  iadi=io+((i3p*(i3p+1))/2)+i4                          9d28s17
                  bc(jda+i3)=bc(jda+i3)+fact*bc(iadi)                   9d28s17
                 end do                                                 9d28s17
                 jda=jda+iacto(isao)                                    9d29s17
                end do                                                  9d28s17
               end if                                                   9d28s17
               io=io+nrow                                               9d28s17
              end do                                                    9d28s17
              i10=1                                                     9d28s17
             end do                                                     9d28s17
            end if                                                      9d28s17
c
c     part b
c
            if(isblk(2,is).eq.isao.and.isblk(4,is).eq.isao)then         10d2s17
             is13=isblk(1,is)                                           10d2s17
             call ilimts(noc(is13),noc(isao),mynprocg,mynowprog,il,     9d27s17
     $            ih,i1s,i1e,i2s,i2e)                                    9d27s17
             i10=i1s                                                    9d27s17
             i1n=noc(is13)                                              9d27s17
             io=ioooo(is)                                               9d27s17
             if(is13.eq.isao)then                                       10d2s17
              nrow=(noc(isao)*(noc(isao)+1))/2                           9d27s17
              iswitch=0                                                 10d2s17
             else                                                       10d2s17
              nrow=noc(isao)*noc(is13)                                  10d2s17
              iswitch=1                                                 10d2s17
             end if                                                     10d2s17
             do i2=i2s,i2e                                              9d27s17
              if(i2.eq.i2e)i1n=i1e                                      9d27s17
              i2m=i2-1                                                  10d2s17
              do i1=i10,i1n                                             10d2s17
               if(i1.gt.idoub(is13).and.i2.le.idoub(isao))then          10d2s17
                i1m=i1-1-idoub(is13)                                    10d2s17
                iden=iden1(is13)+iacto(is13)*i1m                        10d2s17
                jda=ida+iacto(isao)*i2m                                 10d2s17
                do i4=0,iacto(isao)-1                                   10d2s17
                 i4p=i4+idoub(isao)                                     10d2s17
                 kda=jda+i4                                             10d2s17
                 do i3=0,iacto(is13)-1                                  10d2s17
                  i3p=i3+idoub(is13)                                    10d2s17
                  ix=max(i3p,i4p)                                       10d2s17
                  in=min(i3p,i4p)                                       10d2s17
                  ieq=((ix*(ix+1))/2)+in                                10d2s17
                  inot=i3p+noc(is13)*i4p                                10d2s17
                  iadi=io+(inot-ieq)*iswitch+ieq                        10d2s17
                  bc(kda)=bc(kda)+bc(iadi)*bc(iden+i3)                  10d2s17
                 end do                                                 10d2s17
                end do                                                  10d2s17
               end if                                                   10d2s17
               io=io+nrow                                               10d2s17
              end do                                                    10d2s17
              i10=1                                                     10d2s17
             end do                                                     10d2s17
            else if(isblk(1,is).eq.isao.and.isblk(4,is).eq.isao)then    10d2s17
             is23=isblk(2,is)                                           10d2s17
             call ilimts(noc(is23),noc(isao),mynprocg,mynowprog,il,     9d27s17
     $            ih,i1s,i1e,i2s,i2e)                                   10d2s17
             i10=i1s                                                    9d27s17
             i1n=noc(is23)                                              9d27s17
             io=ioooo(is)                                               9d27s17
             nrow=noc(isao)*noc(is23)                                   10d2s17
             do i2=i2s,i2e                                              9d27s17
              if(i2.eq.i2e)i1n=i1e                                      9d27s17
              i2m=i2-1                                                  10d2s17
              do i1=i10,i1n                                             10d2s17
               if(i1.gt.idoub(is23).and.i2.le.idoub(isao))then          10d2s17
                i1m=i1-1-idoub(is23)                                    10d2s17
                iden=iden1(is23)+iacto(is23)*i1m                        10d2s17
                jda=ida+iacto(isao)*i2m                                 10d2s17
                do i4=0,iacto(isao)-1                                   10d2s17
                 i4p=i4+idoub(isao)                                     10d2s17
                 kda=jda+i4                                             10d2s17
                 do i3=0,iacto(is23)-1                                  10d2s17
                  i3p=i3+idoub(is23)                                    10d2s17
                  iadi=io+i4p+noc(isao)*i3p                             10d2s17
                  bc(kda)=bc(kda)+bc(iadi)*bc(iden+i3)                  10d2s17
                 end do                                                 10d2s17
                end do                                                  10d2s17
               end if                                                   10d2s17
               io=io+nrow                                               10d2s17
              end do                                                    10d2s17
              i10=1                                                     10d2s17
             end do                                                     10d2s17
            else if(isblk(1,is).eq.isao.and.isblk(3,is).eq.isao)then    10d2s17
             is24=isblk(2,is)                                           10d2s17
             call ilimts(noc(isao),noc(is24),mynprocg,mynowprog,il,     9d27s17
     $            ih,i1s,i1e,i2s,i2e)                                   10d2s17
             i10=i1s                                                    9d27s17
             i1n=noc(isao)                                              9d27s17
             io=ioooo(is)                                               9d27s17
             nrow=noc(isao)*noc(is24)                                   10d2s17
             do i2=i2s,i2e                                              9d27s17
              if(i2.eq.i2e)i1n=i1e                                      9d27s17
              i2m=i2-1-idoub(is24)                                      10d2s17
              iden=iden1(is24)+iacto(is24)*i2m                          10d2s17
              do i1=i10,i1n                                             10d2s17
               if(i2.gt.idoub(is24).and.i1.le.idoub(isao))then          10d2s17
                i1m=i1-1                                                10d2s17
                jda=ida+iacto(isao)*i1m                                 10d2s17
                do i4=0,iacto(isao)-1                                   10d2s17
                 i4p=i4+idoub(isao)                                     10d2s17
                 kda=jda+i4                                             10d2s17
                 do i3=0,iacto(is24)-1                                  10d2s17
                  i3p=i3+idoub(is24)                                    10d2s17
                  iadi=io+i4p+noc(isao)*i3p                             10d2s17
                  bc(kda)=bc(kda)+bc(iadi)*bc(iden+i3)                  10d2s17
                 end do                                                 10d2s17
                end do                                                  10d2s17
               end if                                                   10d2s17
               io=io+nrow                                               10d2s17
              end do                                                    10d2s17
              i10=1                                                     10d2s17
             end do                                                     10d2s17
            else if(isblk(2,is).eq.isao.and.isblk(3,is).eq.isao)then    10d2s17
             is14=isblk(1,is)                                           10d2s17
             call ilimts(noc(isao),noc(is14),mynprocg,mynowprog,il,     9d27s17
     $            ih,i1s,i1e,i2s,i2e)                                   10d2s17
             i10=i1s                                                    9d27s17
             i1n=noc(isao)                                              9d27s17
             io=ioooo(is)                                               9d27s17
             nrow=noc(isao)*noc(is14)                                   10d2s17
             do i2=i2s,i2e                                              9d27s17
              if(i2.eq.i2e)i1n=i1e                                      9d27s17
              i2m=i2-1-idoub(is14)                                      10d2s17
              iden=iden1(is24)+iacto(is14)*i2m                          10d2s17
              do i1=i10,i1n                                             10d2s17
               if(i2.gt.idoub(is14).and.i1.le.idoub(isao))then          10d2s17
                i1m=i1-1                                                10d2s17
                jda=ida+iacto(isao)*i1m                                 10d2s17
                do i4=0,iacto(isao)-1                                   10d2s17
                 i4p=i4+idoub(isao)                                     10d2s17
                 kda=jda+i4                                             10d2s17
                 do i3=0,iacto(is14)-1                                  10d2s17
                  i3p=i3+idoub(is14)                                    10d2s17
                  iadi=io+i3p+noc(is14)*i4p                             10d2s17
                  bc(kda)=bc(kda)+bc(iadi)*bc(iden+i3)                  10d2s17
                 end do                                                 10d2s17
                end do                                                  10d2s17
               end if                                                   10d2s17
               io=io+nrow                                               10d2s17
              end do                                                    10d2s17
              i10=1                                                     10d2s17
             end do                                                     10d2s17
            end if                                                      10d2s17
c
c     part c
c
            if(isblk(3,is).eq.isao.and.isblk(4,is).eq.isao)then         10d2s17
             is12=isblk(1,is)                                           10d2s17
             call ilimts(noc(isao),noc(isao),mynprocg,mynowprog,il,     9d27s17
     $            ih,i1s,i1e,i2s,i2e)                                   10d2s17
             i10=i1s                                                    9d27s17
             i1n=noc(isao)                                              9d27s17
             io=ioooo(is)-1                                             10d2s17
             nrow=(noc(is12)*(noc(is12)+1))/2                           10d2s17
             do i2=i2s,i2e                                              10d2s17
              i2m=i2-1-idoub(isao)                                      10d2s17
              iden=iden1(isao)+iacto(isao)*i2m                          10d2s17
              if(i2.eq.i2e)i1n=i1e                                      10d2s17
              do i1=i10,i1n                                             10d2s17
               if(i1.le.idoub(isao).and.i2.gt.idoub(isao))then          10d2s17
                i1m=i1-1                                                10d2s17
                jda=ida+iacto(isao)*i1m                                 10d2s17
                do i34=1,idoub(is12)                                    10d2s17
                 iadi=io+(i34*(i34+1))/2                                10d2s17
                 fact=-bc(iadi)*2d0                                     10d2s17
                 do i5=0,iacto(isao)-1                                  10d2s17
                  bc(jda+i5)=bc(jda+i5)+bc(iden+i5)*fact                10d2s17
                 end do                                                 10d2s17
                end do                                                  10d2s17
               end if                                                   10d2s17
               io=io+nrow                                               10d2s17
              end do                                                    10d2s17
              i10=1                                                     10d2s17
             end do                                                     10d2s17
            end if                                                      10d2s17
c
c     part d
c
            if(isblk(1,is).eq.isao.and.isblk(3,is).eq.isao)then         10d2s17
             is24=isblk(2,is)                                           10d2s17
             call ilimts(noc(isao),noc(is24),mynprocg,mynowprog,il,     9d27s17
     $            ih,i1s,i1e,i2s,i2e)                                   10d2s17
             i10=i1s                                                    9d27s17
             i1n=noc(isao)                                              9d27s17
             if(is24.eq.isao)then                                       10d2s17
              nrow=(noc(isao)*(noc(isao)+1))/2                          10d2s17
              iswitch=0                                                 10d2s17
             else                                                       10d2s17
              nrow=noc(isao)*noc(is24)                                  10d2s17
              iswitch=1                                                 10d2s17
             end if                                                     10d2s17
             io=ioooo(is)                                               10d2s17
             do i2=i2s,i2e                                              10d2s17
              i2m=i2-1                                                  10d2s17
              if(i2.eq.i2e)i1n=i1e                                      10d2s17
              do i1=i10,i1n                                             10d2s17
               if(i1.gt.idoub(isao).and.i2.le.idoub(is24))then          10d2s17
                i1m=i1-1-idoub(isao)                                    10d2s17
                iden=iden1(isao)+iacto(isao)*i1m                        10d2s17
                do i3=0,idoub(isao)-1                                   10d2s17
                 ix=max(i3,i2m)                                         10d2s17
                 in=min(i3,i2m)                                         10d2s17
                 ieq=((ix*(ix+1))/2)+in                                 10d2s17
                 inot=i3+noc(isao)*i2m                                  10d2s17
                 iadi=io+(inot-ieq)*iswitch+ieq                         10d2s17
                 jda=ida+iacto(isao)*i3                                 10d2s17
                 do i5=0,iacto(isao)-1                                  10d2s17
                  bc(jda+i5)=bc(jda+i5)+bc(iadi)*bc(iden+i5)            10d2s17
                 end do                                                 10d2s17
                end do                                                  10d2s17
               end if                                                   10d2s17
               io=io+nrow                                               10d2s17
              end do                                                    10d2s17
              i10=1                                                     10d2s17
             end do                                                     10d2s17
            else if(isblk(2,is).eq.isao.and.isblk(3,is).eq.isao)then    11d20s17
             is14=isblk(1,is)                                           11d20s17
             call ilimts(noc(isao),noc(is14),mynprocg,mynowprog,il,     11d20s17
     $            ih,i1s,i1e,i2s,i2e)                                   10d2s17
             i10=i1s                                                    9d27s17
             i1n=noc(isao)                                              9d27s17
             nrow=noc(isao)*noc(is14)                                   11d20s17
             io=ioooo(is)                                               10d2s17
             do i2=i2s,i2e                                              10d2s17
              i2m=i2-1                                                  10d2s17
              if(i2.eq.i2e)i1n=i1e                                      10d2s17
              do i1=i10,i1n                                             10d2s17
               if(i1.gt.idoub(isao).and.i2.le.idoub(is14))then          11d20s17
                i1m=i1-1-idoub(isao)                                    10d2s17
                iden=iden1(isao)+iacto(isao)*i1m                        10d2s17
                do i3=0,idoub(isao)-1                                   10d2s17
                 iadi=io+i2m+noc(is14)*i3                               11d27s17
                 jda=ida+iacto(isao)*i3                                 10d2s17
                 do i5=0,iacto(isao)-1                                  10d2s17
                  bc(jda+i5)=bc(jda+i5)+bc(iadi)*bc(iden+i5)            10d2s17
                 end do                                                 10d2s17
                end do                                                  10d2s17
               end if                                                   10d2s17
               io=io+nrow                                               10d2s17
              end do                                                    10d2s17
              i10=1                                                     10d2s17
             end do                                                     10d2s17
            else if(isblk(2,is).eq.isao.and.isblk(4,is).eq.isao)then    11d27s17
             is13=isblk(1,is)                                           11d27s17
             call ilimts(noc(is13),noc(isao),mynprocg,mynowprog,il,     11d27s17
     $            ih,i1s,i1e,i2s,i2e)                                   10d2s17
             i10=i1s                                                    9d27s17
             i1n=noc(is13)                                              11d27s17
             nrow=noc(isao)*noc(is13)                                   11d27s17
             io=ioooo(is)                                               10d2s17
             do i2=i2s,i2e                                              10d2s17
              i2m=i2-1-idoub(isao)                                      11d27s17
              iden=iden1(isao)+iacto(isao)*i2m                          11d27s17
              if(i2.eq.i2e)i1n=i1e                                      10d2s17
              do i1=i10,i1n                                             10d2s17
               if(i2.gt.idoub(isao).and.i1.le.idoub(is13))then          11d27s17
                i1m=i1-1                                                11d27s17
                do i3=0,idoub(isao)-1                                   10d2s17
                 iadi=io+i1m+noc(is13)*i3                               11d27s17
                 jda=ida+iacto(isao)*i3                                 10d2s17
                 do i5=0,iacto(isao)-1                                  10d2s17
                  bc(jda+i5)=bc(jda+i5)+bc(iadi)*bc(iden+i5)            10d2s17
                 end do                                                 10d2s17
                end do                                                  10d2s17
               end if                                                   10d2s17
               io=io+nrow                                               10d2s17
              end do                                                    10d2s17
              i10=1                                                     10d2s17
             end do                                                     10d2s17
            else if(isblk(1,is).eq.isao.and.isblk(4,is).eq.isao)then    11d27s17
             is23=isblk(2,is)                                           11d27s17
             call ilimts(noc(is23),noc(isao),mynprocg,mynowprog,il,     11d27s17
     $            ih,i1s,i1e,i2s,i2e)                                   10d2s17
             i10=i1s                                                    9d27s17
             i1n=noc(is23)                                              11d27s17
             nrow=noc(isao)*noc(is23)                                   11d27s17
             io=ioooo(is)                                               10d2s17
             do i2=i2s,i2e                                              10d2s17
              i2m=i2-1-idoub(isao)                                      11d27s17
              iden=iden1(isao)+iacto(isao)*i2m                          11d27s17
              if(i2.eq.i2e)i1n=i1e                                      10d2s17
              do i1=i10,i1n                                             10d2s17
               if(i2.gt.idoub(isao).and.i1.le.idoub(is23))then          11d27s17
                i1m=i1-1                                                11d27s17
                do i3=0,idoub(isao)-1                                   10d2s17
                 iadi=io+i3+noc(isao)*i1m                               11d27s17
                 jda=ida+iacto(isao)*i3                                 10d2s17
                 do i5=0,iacto(isao)-1                                  10d2s17
                  bc(jda+i5)=bc(jda+i5)+bc(iadi)*bc(iden+i5)            10d2s17
                 end do                                                 10d2s17
                end do                                                  10d2s17
               end if                                                   10d2s17
               io=io+nrow                                               10d2s17
              end do                                                    10d2s17
              i10=1                                                     10d2s17
             end do                                                     10d2s17
            end if                                                      10d2s17
           end do                                                       9d28s17
          end if                                                        9d28s17
          if(idoit(3).ne.0)then                                         9d27s17
           do is=1,nsdlk                                                9d27s17
c
c     dddd j contribution to hessc from dv av.
c
            if(isblk(1,is).eq.isao.and.isblk(2,is).eq.isao)then         9d27s17
             is34=isblk(3,is)                                           9d27s17
             call ilimts(noc(is34),noc(is34),mynprocg,mynowprog,il,     9d27s17
     $            ih,i1s,i1e,i2s,i2e)                                    9d27s17
             i10=i1s                                                    9d27s17
             i1n=noc(is34)                                              9d27s17
             io=ioooo(is)                                               9d27s17
             nrow=(noc(isao)*(noc(isao)+1))/2                           9d27s17
             do i2=i2s,i2e                                              9d27s17
              if(i2.eq.i2e)i1n=i1e                                      9d27s17
              i2m=i2-1                                                  9d27s17
              do i1=i10,i1n
               if(i1.le.idoub(is34).and.i1.eq.i2)then                   10d2s17
                jda=ida                                                 9d27s17
                do i4=0,idoub(isao)-1                                   9d27s17
                 do i3=0,iacto(isao)-1                                  9d27s17
                  i3p=i3+idoub(isao)                                    9d27s17
                  iadi=io+((i3p*(i3p+1))/2)+i4                          9d27s17
                  bc(jda+i3)=bc(jda+i3)-4d0*bc(iadi)                    9d27s17
                 end do                                                 9d27s17
                 jda=jda+iacto(isao)                                    9d27s17
                end do                                                  9d27s17
               end if                                                   9d27s17
               io=io+nrow                                               9d27s17
              end do                                                    9d27s17
              i10=1                                                     9d27s17
             end do                                                     9d27s17
            end if                                                      9d27s17
c
c     dddd k contribution to hessc from dv av.
c
            if(isblk(1,is).eq.isao.and.isblk(3,is).eq.isao)then         9d27s17
             is24=isblk(2,is)                                           9d27s17
             iswitch=iabs(is24-isao)                                    9d27s17
             iswitch=iswitch/max(1,iswitch)                             9d27s17
             call ilimts(noc(isao),noc(is24),mynprocg,mynowprog,il,     9d27s17
     $            ih,i1s,i1e,i2s,i2e)                                    9d27s17
             i10=i1s                                                    9d27s17
             i1n=noc(isao)                                              9d27s17
             io=ioooo(is)                                               9d27s17
             if(is24.eq.isao)then                                       9d27s17
              nrow=(noc(isao)*(noc(isao)+1))/2                           9d27s17
             else                                                       9d27s17
              nrow=noc(isao)*noc(is24)                                  9d27s17
             end if                                                     9d27s17
             do i2=i2s,i2e                                              9d27s17
              i2m=i2-1                                                  9d27s17
              if(i2.eq.i2e)i1n=i1e                                      9d27s17
              do i1=i10,i1n                                             9d27s17
               if(i1.le.idoub(isao).and.i2.le.idoub(is24))then          9d27s17
                i1m=i1-1                                                9d27s17
                jda=ida+iacto(isao)*i1m                                 9d27s17
                do i3=0,iacto(isao)-1                                   9d27s17
                 i3p=i3+idoub(isao)                                     9d27s17
                 ix=max(i3p,i2m)                                        9d27s17
                 in=min(i3p,i2m)                                        9d27s17
                 ieq=((ix*(ix+1))/2)+in                                 9d27s17
                 inot=i3p+noc(isao)*i2m                                 9d27s17
                 iadi=io+(inot-ieq)*iswitch+ieq                         9d27s17
                 bc(jda+i3)=bc(jda+i3)+2d0*bc(iadi)                     9d27s17
                end do                                                  9d27s17
               end if                                                   9d27s17
               io=io+nrow                                               9d27s17
              end do                                                    9d27s17
              i10=1                                                     9d27s17
             end do                                                     9d27s17
            else if(isblk(2,is).eq.isao.and.isblk(3,is).eq.isao)then    9d27s17
             is14=isblk(1,is)                                           9d27s17
             call ilimts(noc(isao),noc(is14),mynprocg,mynowprog,il,     9d27s17
     $            ih,i1s,i1e,i2s,i2e)                                    9d27s17
             i10=i1s                                                    9d27s17
             i1n=noc(isao)                                              9d27s17
             io=ioooo(is)                                               9d27s17
             nrow=noc(isao)*noc(is14)                                   9d27s17
             do i2=i2s,i2e                                              9d27s17
              i2m=i2-1                                                  9d27s17
              if(i2.eq.i2e)i1n=i1e                                      9d27s17
              do i1=i10,i1n                                             9d27s17
               if(i1.le.idoub(isao).and.i2.le.idoub(is14))then          11d9s17
                i1m=i1-1                                                9d27s17
                jda=ida+iacto(isao)*i1m                                 9d27s17
                do i3=0,iacto(isao)-1                                   9d27s17
                 i3p=i3+idoub(isao)                                     9d27s17
                 iadi=io+i2m+noc(is14)*i3p                              11d9s17
                 bc(jda+i3)=bc(jda+i3)+2d0*bc(iadi)                     9d27s17
                end do                                                  9d27s17
               end if                                                   9d27s17
               io=io+nrow                                               9d27s17
              end do                                                    9d27s17
              i10=1                                                     9d27s17
             end do                                                     9d27s17
            else if(isblk(2,is).eq.isao.and.isblk(4,is).eq.isao)then    9d27s17
             is13=isblk(1,is)                                           9d27s17
             call ilimts(noc(is13),noc(isao),mynprocg,mynowprog,il,     9d27s17
     $            ih,i1s,i1e,i2s,i2e)                                    9d27s17
             i10=i1s                                                    9d27s17
             i1n=noc(is13)                                              9d27s17
             io=ioooo(is)                                               9d27s17
             nrow=noc(isao)*noc(is13)                                   9d27s17
             do i2=i2s,i2e                                              9d27s17
              i2m=i2-1                                                  9d27s17
              if(i2.eq.i2e)i1n=i1e                                      9d27s17
              do i1=i10,i1n                                             9d27s17
               if(i1.le.idoub(is13).and.i2.le.idoub(isao))then          9d27s17
                i1m=i1-1                                                9d27s17
                jda=ida+iacto(isao)*i2m                                 9d27s17
                do i3=0,iacto(isao)-1                                   9d27s17
                 i3p=i3+idoub(isao)                                     9d27s17
                 iadi=io+i1m+noc(is13)*i3p                              11d9s17
                 bc(jda+i3)=bc(jda+i3)+2d0*bc(iadi)                     9d27s17
                end do                                                  9d27s17
               end if                                                   9d27s17
               io=io+nrow                                               9d27s17
              end do                                                    9d27s17
              i10=1                                                     9d27s17
             end do                                                     9d27s17
            else if(isblk(1,is).eq.isao.and.isblk(4,is).eq.isao)then    9d27s17
             is23=isblk(2,is)                                           11d9s17
             call ilimts(noc(is23),noc(isao),mynprocg,mynowprog,il,     9d27s17
     $            ih,i1s,i1e,i2s,i2e)                                    9d27s17
             i10=i1s                                                    9d27s17
             i1n=noc(is23)                                              11d9s17
             io=ioooo(is)                                               9d27s17
             nrow=noc(isao)*noc(is23)                                   9d27s17
             do i2=i2s,i2e                                              9d27s17
              i2m=i2-1                                                  9d27s17
              if(i2.eq.i2e)i1n=i1e                                      9d27s17
              do i1=i10,i1n                                             9d27s17
               if(i1.le.idoub(is23).and.i2.le.idoub(isao))then          11d9s17
                i1m=i1-1                                                9d27s17
                jda=ida+iacto(isao)*i2m                                 9d27s17
                do i3=0,iacto(isao)-1                                   9d27s17
                 i3p=i3+idoub(isao)                                     9d27s17
                 iadi=io+i3p+noc(isao)*i1m                              9d27s17
                 bc(jda+i3)=bc(jda+i3)+2d0*bc(iadi)                     9d27s17
                end do                                                  9d27s17
               end if                                                   9d27s17
               io=io+nrow                                               9d27s17
              end do                                                    9d27s17
              i10=1                                                     9d27s17
             end do                                                     9d27s17
            end if                                                      9d27s17
           end do                                                       9d27s17
          end if                                                        9d27s17
          if(idoit(3).ne.0)then
           do is=1,nsdlk
           if(isblk(1,is).eq.isao.and.isblk(2,is).eq.isao)then
            iscv=isblk(3,is)
            isco=multh(iscv,ipsym)
            call ilimts(noc(iscv),noc(iscv),mynprocg,mynowprog,il,
     $           ih,i1s,i1e,i2s,i2e)
            ii=ioooo(is)
            i10=i1s
            i1n=noc(iscv)
            do i2=i2s,i2e
             if(i2.eq.i2e)i1n=i1e
             do i1=i10,i1n
              if(i1.eq.i2)then
               if(i1.le.idoub(iscv))then                                4d11s17
               do i34=0,noo-1
                bc(ioop+i34)=bc(ioop+i34)-8d0*bc(ii+i34)
               end do
               end if                                                   4d11s17
              end if
              ii=ii+noox                                                7d19s17
             end do
             i10=1
            end do
           end if
          end do
          if(idwsdeb.gt.10)then
          write(6,*)('oop after i1=i2 sum ')
          write(6,*)(bc(ioop+i),i=0,noo-1)
          end if
          do isc=1,nsymb
           if(noc(isc).gt.0)then
            do is=1,nsdlk
             if(isblk(1,is).eq.isc.and.isblk(2,is).eq.isao.and.
     $          isblk(3,is).eq.isc)then                                 3d29s16
              call ilimts(noc(isc),noc(isao),mynprocg,mynowprog,il,     3d29s16
     $           ih,i1s,i1e,i2s,i2e)
              ii=ioooo(is)
              if(isc.eq.isao)then
               noo=(noc(isao)*(noc(isao)+1))/2
              else
               noo=noc(isao)*noc(isc)
              end if
              i10=i1s
              i1n=noc(isc)                                              3d29s16
              do i2=i2s,i2e
               if(i2.eq.i2e)i1n=i1e
               i2m=i2-1
               do i1=i10,i1n
                if(i2.le.idoub(isao).and.i1.le.idoub(isc))then          4d27s17
                 i1m=i1-1
                 if(isao.eq.isc)then
                  do il=i2m,idoub(isao)-1                                4d11s17
                   iad=ioop+((il*(il+1))/2)+i2m
                   ix=max(il,i1m)
                   in=min(il,i1m)
                   iado=ii+((ix*(ix+1))/2)+in
                   bc(iad)=bc(iad)+4d0*bc(iado)
                  end do
                 else
                  do il=i2m,idoub(isao)-1                               4d27s17
                   iad=ioop+((il*(il+1))/2)+i2m
                   iado=ii+i1m+noc(isc)*il
                   bc(iad)=bc(iad)+4d0*bc(iado)
                  end do
                 end if
                end if                                                  4d27s17
                ii=ii+noo
               end do
               i10=1
              end do
              go to 4
             end if
             if(isblk(1,is).eq.isao.and.isblk(2,is).eq.isc.and.
     $          isblk(3,is).eq.isao)then                                3d29s16
              call ilimts(noc(isao),noc(isc),mynprocg,mynowprog,il,     3d29s16
     $           ih,i1s,i1e,i2s,i2e)
              ii=ioooo(is)
              if(isc.eq.isao)then
               noo=(noc(isao)*(noc(isao)+1))/2
              else
               noo=noc(isao)*noc(isc)
              end if
              i10=i1s
              i1n=noc(isao)                                             3d29s16
              do i2=i2s,i2e
               if(i2.eq.i2e)i1n=i1e
               i2m=i2-1
               do i1=i10,i1n
                if(i2.le.idoub(isc).and.i1.le.idoub(isao))then          4d27s17
                 i1m=i1-1
                 if(isao.eq.isc)then
                  do il=i1m,idoub(isao)-1                               4d27s17
                   iad=ioop+((il*(il+1))/2)+i1m
                   ix=max(il,i2m)
                   in=min(il,i2m)
                   iado=ii+((ix*(ix+1))/2)+in
                   bc(iad)=bc(iad)+4d0*bc(iado)
                  end do
                 else
                  do il=i1m,idoub(isao)-1
                   iad=ioop+((il*(il+1))/2)+i1m
                   iado=ii+il+noc(isao)*i2m
                   bc(iad)=bc(iad)+4d0*bc(iado)
                  end do
                 end if
                end if                                                  4d27s17
                ii=ii+noo
               end do
               i10=1
              end do
              go to 4
             end if
             if(isblk(1,is).eq.isao.and.isblk(2,is).eq.isc.and.
     $          isblk(3,is).eq.isc)then                                 3d29s16
              call ilimts(noc(isc),noc(isao),mynprocg,mynowprog,il,     3d29s16
     $           ih,i1s,i1e,i2s,i2e)
              ii=ioooo(is)
              if(isc.eq.isao)then
               noo=(noc(isao)*(noc(isao)+1))/2
              else
               noo=noc(isao)*noc(isc)
              end if
              i10=i1s
              i1n=noc(isc)                                              3d29s16
              do i2=i2s,i2e
               if(i2.eq.i2e)i1n=i1e
               i2m=i2-1
               do i1=i10,i1n
                if(i2.le.idoub(isao).and.i1.le.idoub(isc))then          4d27s17
                 i1m=i1-1
                 if(isao.eq.isc)then
                  do il=i2m,idoub(isao)-1
                   iad=ioop+((il*(il+1))/2)+i2m
                   ix=max(il,i1m)
                   in=min(il,i1m)
                   iado=ii+((ix*(ix+1))/2)+in
                   bc(iad)=bc(iad)+4d0*bc(iado)
                  end do
                 else
                  do il=i2m,idoub(isao)-1
                   iad=ioop+((il*(il+1))/2)+i2m
                   iado=ii+il+noc(isao)*i2m
                   bc(iad)=bc(iad)+4d0*bc(iado)
                  end do
                 end if
                end if                                                  4d27s17
                ii=ii+noo
               end do
               i10=1
              end do
              go to 4
             end if
             if(isblk(1,is).eq.isc.and.isblk(2,is).eq.isao.and.          3d29s16
     $          isblk(3,is).eq.isao)then                                3d29s16
              call ilimts(noc(isao),noc(isc),mynprocg,mynowprog,il,     3d29s16
     $           ih,i1s,i1e,i2s,i2e)
              ii=ioooo(is)
              noo=noc(isao)*noc(isc)
              i10=i1s
              i1n=noc(isao)                                             3d29s16
              do i2=i2s,i2e
               if(i2.eq.i2e)i1n=i1e
               i2m=i2-1
               do i1=i10,i1n
                if(i2.le.idoub(isc).and.i1.le.idoub(isao))then          4d27s17
                 i1m=i1-1
                 do il=i1m,idoub(isao)-1                                4d27s17
                  iad=ioop+((il*(il+1))/2)+i1m
                  iado=ii+i2m+il*noc(isc)
                  bc(iad)=bc(iad)+4d0*bc(iado)
                 end do
                end if                                                  4d27s17
                ii=ii+noo
               end do
               i10=1
              end do
              go to 4
             end if
            end do
         write(6,*)('in buildhess')
            write(6,*)('missing isc: '),isc,isao,isc,isao               3d29s16
            do is=1,nsdlk
             write(6,*)(isblk(j,is),j=1,4)
            end do
            call dws_sync
            call dws_finalize
            stop
    4       continue
           end if
          end do
          end if                                                        9d25s17
          if(idoit(5).ne.0)then                                         10d19s17
c
c     aaaa vd vn contribution                                           10d19s17
c
           do isd=1,nsdlk                                               10d19s17
            if(jdenpt(isd).gt.0.and.isblk(1,isd).eq.isbo)then           10d19s17
             is2=isblk(2,isd)                                           10d19s17
             is3=isblk(3,isd)                                           10d19s17
             is4=isblk(4,isd)                                           10d19s17
             if(isbo.eq.is2)then                                        10d19s17
              jswitch=0                                                 10d19s17
              nrowd=(iacto(isbo)*(iacto(isbo)+1))/2                     10d19s17
             else                                                       10d19s17
              jswitch=1                                                 10d19s17
              nrowd=iacto(isbo)*iacto(is2)                              10d19s17
             end if                                                     10d19s17
             do isi=1,nsdlk                                              10d19s17
              if(isblk(1,isi).eq.isao.and.is2.eq.isblk(2,isi).and.      10d19s17
     $           is3.eq.isblk(3,isi))then                               10d19s17
               call ilimts(noc(is3),noc(is4),mynprocg,mynowprog,il,ih,  10d19s17
     $             i1s,i1e,i2s,i2e)                                     10d19s17
               if(isao.eq.is2)then                                      10d19s17
                iswitch=0                                               10d19s17
                nrowi=(noc(isao)*(noc(isao)+1))/2                       10d19s17
               else                                                     10d19s17
                iswitch=1                                               10d19s17
                nrowi=noc(isao)*noc(is2)                                10d19s17
               end if                                                   10d19s17
               i10=i1s                                                  10d19s17
               i1n=noc(is3)                                             10d19s17
               i4o=ioooo(isi)                                           10d19s17
               do i2=i2s,i2e                                            10d19s17
                if(i2.eq.i2e)i1n=i1e                                    10d19s17
                i2m=i2-1-idoub(is4)                                     10d19s17
                do i1=i10,i1n                                           10d19s17
                 if(i1.gt.idoub(is3).and.i2.gt.idoub(is4))then          10d19s17
                  if(isblk(1,isd).ne.isblk(2,isd))then                  11d27s17
                   factr=-1d0                                           11d27s17
                  else                                                  11d27s17
                   factr=-2d0                                           11d27s17
                  end if                                                11d27s17
                  i1m=i1-1-idoub(is3)                                   10d19s17
                  ix=max(i1m,i2m)                                          10d19s17
                  in=min(i1m,i2m)                                          10d19s17
                  ieq=((ix*(ix+1))/2)+in                                10d19s17
                  inot=i1m+iacto(is3)*i2m                               10d19s17
                  iden=jdenpt(isd)+nrowd*((inot-ieq)*jswitch+ieq)       10d19s17
                  do i4=0,iacto(is2)-1                                   10d19s17
                   i4p=i4+idoub(is2)                                     10d19s17
                   do i3=0,idoub(isao)-1                                 10d19s17
                    ix=max(i3,i4p)                                       10d19s17
                    in=min(i3,i4p)                                       10d19s17
                    ieq=((ix*(ix+1))/2)+in                               10d19s17
                    inot=i3+noc(isao)*i4p                                10d19s17
                    iadi=i4o+(inot-ieq)*iswitch+ieq                      10d19s17
                    do i5=0,iacto(isbo)-1                                10d19s17
                     jda=ida+i5+iacto(isbo)*i3                          10d19s17
                     ix=max(i5,i4)                                      10d19s17
                     in=min(i5,i4)                                      10d19s17
                     ieq=((ix*(ix+1))/2)+in                             10d19s17
                     inot=i5+iacto(isbo)*i4                             10d19s17
                     jden=iden+(inot-ieq)*jswitch+ieq                   10d19s17
                     bc(jda)=bc(jda)+factr*bc(iadi)*bc(jden)            11d27s17
                    end do                                               10d19s17
                   end do                                                10d19s17
                  end do                                                 10d19s17
                 end if                                                 10d19s17
                 i4o=i4o+nrowi                                          10d19s17
                end do                                                  10d19s17
                i10=1                                                   10d19s17
               end do                                                   10d19s17
              end if                                                    10d19s17
             end do                                                      10d19s17
            else if(jdenpt(isd).gt.0.and.isblk(2,isd).eq.isbo)then      10d20s17
             is1=isblk(1,isd)                                           10d19s17
             is3=isblk(3,isd)                                           10d19s17
             is4=isblk(4,isd)                                           10d19s17
             nrowd=iacto(isbo)*iacto(is1)                               10d20s17
             do isi=1,nsdlk                                             10d20s17
              if(isblk(1,isi).eq.is1.and.isblk(2,isi).eq.isao.and.      10d20s17
     $             isblk(3,isi).eq.is3)then                             10d20s17
               call ilimts(noc(is3),noc(is4),mynprocg,mynowprog,il,ih,  10d19s17
     $             i1s,i1e,i2s,i2e)                                     10d19s17
               nrowi=noc(isao)*noc(is1)                                 10d20s17
               i10=i1s                                                  10d19s17
               i1n=noc(is3)                                             10d19s17
               i4o=ioooo(isi)                                           10d19s17
               do i2=i2s,i2e                                            10d19s17
                if(i2.eq.i2e)i1n=i1e                                    10d19s17
                i2m=i2-1-idoub(is4)                                     10d19s17
                do i1=i10,i1n                                           10d19s17
                 if(i1.gt.idoub(is3).and.i2.gt.idoub(is4))then          10d19s17
                  i1m=i1-1-idoub(is3)                                   10d19s17
                  inot=i1m+iacto(is3)*i2m                               10d19s17
                  iden=jdenpt(isd)+nrowd*inot                           10d20s17
                  do i4=0,iacto(is1)-1                                  10d20s17
                   i4p=i4+idoub(is1)                                    10d20s17
                   do i3=0,idoub(isao)-1                                 10d19s17
                    inot=i4p+noc(is1)*i3                                10d20s17
                    iadi=i4o+inot                                       10d20s17
                    do i5=0,iacto(isbo)-1                                10d19s17
                     jda=ida+i5+iacto(isbo)*i3                          10d19s17
                     inot=i4+iacto(is1)*i5                              10d20s17
                     jden=iden+inot                                     10d20s17
                     bc(jda)=bc(jda)-bc(iadi)*bc(jden)                  10d20s17
                    end do                                               10d19s17
                   end do                                                10d19s17
                  end do                                                 10d19s17
                 end if                                                 10d19s17
                 i4o=i4o+nrowi                                          10d19s17
                end do                                                  10d19s17
                i10=1                                                   10d19s17
               end do                                                   10d19s17
              end if                                                    10d20s17
             end do                                                     10d20s17
            end if                                                      10d19s17
           end do                                                       10d19s17
          end if                                                        10d19s17
          noo=(idoub(isao)*(idoub(isao)+1))/2                           4d27s17
          if(idwsdeb.gt.0)then
          write(6,*)('oop after i1=i3 sum ')
          write(6,*)(bc(ioop+i),i=0,noo-1)
          end if
          call dws_gsumf(bc(ioop),noo)
          call dws_gsumf(bc(ida),nda)                                   9d27s17
          if(idwsdeb.gt.10)then
          write(6,*)('2e part of oop: ')
          write(6,*)(bc(ioop+i),i=0,noo-1)
          end if
          joop=ioop
          nh=nbasdws(isao)*nvirtc(isao)/nvirt(isao)                     9d25s17
          if(idoit(1).ne.0)then
c
c     dddd vd vn contribution
c
           jda=ida                                                      9d27s17
           do i3=0,idoub(isao)-1                                        9d27s17
            ih0=ipt(isao)+idoub(isao)+1+nh*i3                           9d27s17
            do i4=0,iacto(isao)-1                                       9d27s17
             bc(jda+i4)=bc(jda+i4)-2d0*h0mo(ih0+i4)                     9d27s17
            end do                                                      9d27s17
            jda=jda+iacto(isao)                                         9d27s17
           end do                                                       9d27s17
          end if                                                        10d2s17
          if(idoit(2).ne.0)then                                         10d2s17
c
c     aa vd vn contribution
c
           do i1=0,idoub(isao)-1                                        10d2s17
            iamo=ipt(isao)+nh*i1+idoub(isao)+1                          10d2s17
            jda=ida+iacto(isao)*i1                                      10d2s17
            do i2=0,iacto(isao)-1                                       10d2s17
             iden=iden1(isao)+iacto(isao)*i2                            10d2s17
             fact=h0mo(iamo+i2)                                         10d2s17
             do i3=0,iacto(isao)-1                                      10d2s17
              bc(jda+i3)=bc(jda+i3)-bc(iden+i3)*fact                    10d2s17
             end do                                                     10d2s17
            end do                                                      10d2s17
           end do                                                       10d2s17
          end if                                                        10d2s17
          if(idoit(1).ne.0)then                                         9d25s17
          do i3=0,idoub(isao)-1                                         4d27s17
           do i4=0,i3
            iadh=ipt(isao)+i4+nh*i3+1
            bc(joop)=bc(joop)-4d0*h0mo(iadh)
            joop=joop+1
           end do
          end do
          end if                                                        9d25s17
          if(idwsdeb.gt.10)then
          write(6,*)('full ioop ')
          write(6,*)(bc(ioop+i),i=0,noo-1)
          end if
          if(idoit(3).ne.0)then                                         9d25s17
          do is=1,nsdlk
           if(isblk(3,is).eq.isav.and.isblk(4,is).eq.isav)then          3d29s16
            isc=isblk(1,is)
            call ilimts(nvirtc(isav),nvirtc(isav),mynprocg,mynowprog,il,3d29s16
     $           ih,i1s,i1e,i2s,i2e)
            ii=jmats(is)-1
            i10=i1s
            i1n=nvirtc(isav)                                            3d29s16
            nrowj=(noc(isc)*(noc(isc)+1))/2
            do i2=i2s,i2e
             if(i2.eq.i2e)i1n=i1e
             do i1=i10,i1n
              if(i1.ge.i2)then
               jvv=ivv+((i1*(i1-1))/2)+i2-1
               do i34=1,idoub(isc)                                      4d27s17
                iadj=ii+(i34*(i34+1))/2
                bc(jvv)=bc(jvv)+8d0*bc(iadj)
               end do
              end if
              ii=ii+nrowj
             end do
             i10=1
            end do
           end if
          end do
          if(idwsdeb.gt.10)then
          write(6,*)('vv after 34 sum ')
          write(6,*)(bc(ivv+i),i=0,nvv-1)
          end if
          do is=1,nsdlkk
           if(isblkk(3,is).eq.isav.and.isblkk(4,is).eq.isav)then        3d29s16
            isc=isblkk(1,is)
            call ilimts(nvirtc(isav),nvirtc(isav),mynprocg,mynowprog,il,3d29s16
     $           ih,i1s,i1e,i2s,i2e)
            ii=kmats(is)
            i10=i1s
            i1n=nvirtc(isav)                                            3d29s16
            nrowk=noc(isc)*noc(isc)
            do i2=i2s,i2e
             if(i2.eq.i2e)i1n=i1e
             do i1=i10,i1n
              if(i1.ge.i2)then
               jvv=ivv+((i1*(i1-1))/2)+i2-1
               do i34=0,idoub(isc)-1                                    4d27s17
                iadk=ii+i34*(noc(isc)+1)
                bc(jvv)=bc(jvv)-4d0*bc(iadk)
               end do
              end if
              ii=ii+nrowk
             end do
             i10=1
            end do
           end if
          end do
          end if                                                        9d25s17
          call dws_gsumf(bc(ivv),nvv)
          nh=noc(isav)+nvirtc(isav)                                     3d29s16
          i10=i1sh
          i1n=nvirtc(isbv)
          ii=ihessc(isa,isb)
          do i2=i2sh,i2eh
           if(i2.eq.i2eh)i1n=i1eh
           do i1=i10,i1n
            ix=max(i1,i2)
            in=min(i1,i2)
            jvv=ivv+((ix*(ix-1))/2)+in-1
            iadh0=ipt(isav)+i2+noc(isav)+nh*(i1+noc(isav)-1)
            if(idoit(1).ne.0)then                                       9d25s17
             h0add=4d0*h0mo(iadh0)
            else                                                        9d25s17
             h0add=0d0                                                  9d25s17
            end if                                                      9d25s17
            do i34=0,idoub(isao)-1                                      4d27s17
             iadh=ii+i34*(noc(isao)+1)
             bc(iadh)=bc(iadh)+bc(jvv)+h0add
            end do
            if(i1.eq.i2)then
             jda=ida                                                    9d27s17
             do i4=0,idoub(isao)-1                                      9d27s17
              do i3=0,iacto(isao)-1                                     9d27s17
               i3p=i3+idoub(isao)                                       9d27s17
               iad=ii+i3p+noc(isao)*i4                                  9d27s17
               bc(iad)=bc(iad)+bc(jda+i3)                               9d27s17
               iad=ii+i4+noc(isao)*i3p                                  9d27s17
               bc(iad)=bc(iad)+bc(jda+i3)                               9d27s17
              end do                                                    9d27s17
              jda=jda+iacto(isao)                                       9d27s17
             end do                                                     9d27s17
             do i4=0,idoub(isao)-1                                      4d27s17
              do i3=0,idoub(isao)-1                                     4d27s17
               ix=max(i3,i4)
               in=min(i3,i4)
               joop=ioop+((ix*(ix+1))/2)+in
               bc(ii)=bc(ii)+bc(joop)
               ii=ii+1
              end do
              ii=ii+iacto(isao)                                         7d19s17
             end do
             ii=ii+iacto(isao)*noc(isao)                                7d19s17
            else
             ii=ii+nrowh
            end if
           end do
           i10=1
          end do
          ioffh=ioffh+nh*nh
         end if
         if(isa.eq.isb)then                                             12d4s17
          i10=i1sh                                                      12d4s17
          i1n=nvirtc(isbv)                                              12d4s17
          ii=ihessc(isa,isb)                                            12d4s17
          ihessdo=ihessd(isa)+idoub(isao)*iacto(isa)                    12d4s17
          nocp=noc(isao)+1                                              12d4s17
          do i2=i2sh,i2eh                                               12d4s17
           if(i2.eq.i2eh)i1n=i1eh                                       12d4s17
           i2m=i2-1                                                     12d4s17
           do i1=i10,i1n                                                12d4s17
            if(i1.eq.i2)then                                            12d4s17
             do i34=0,noc(isao)-1                                       12d4s17
              iad=ihessdo+i2m+nvirtc(isa)*i34                           12d4s17
              iad2=ii+i34*nocp                                          12d4s17
              bc(iad)=bc(iad2)                                          12d4s17
             end do                                                     12d4s17
            end if                                                      12d4s17
            ii=ii+nrowh                                                 12d4s17
           end do                                                       12d4s17
           i10=1                                                        12d4s17
          end do                                                        12d4s17
          call dws_gsumf(bc(ihessd(isa)),nhessd)                        12d4s17
          idp=idoub(isao)+1                                             12d4s17
          iap=idoub(isao)*idoub(isao)*(iacto(isa)+1)                    12d4s17
          do ia=0,iacto(isa)-1                                          12d4s17
           do id=0,idoub(isao)-1                                        12d4s17
            iad1=ihessd(isa)+ia+iacto(isa)*id                           12d4s17
            iad2=ihessa(isa,isa)+id*idp+ia*iap                          12d4s17
            bc(iad1)=bc(iad2)                                           12d4s17
           end do                                                       12d4s17
          end do                                                        12d4s17
         end if                                                         12d4s17
         if(idwsdeb.gt.10.and.lsym)then                                 11d30s17
          rms=0d0
          do i=0,nrowh*nhere-1
           rms=rms+bc(ihessc(isa,isb)+i)**2
          end do
          if(rms.ne.0d0)then
           rms=sqrt(rms/dfloat(nrowh*nhere))
          end if
          write(6,*)('hessc rms: '),rms,isao,isbo,isbv,isav
          if(rms.gt.1d-5)then
         write(6,*)('hessian for symmetry block '),isao,isbo,isbv,isav
     $        ,ihessc(isa,isb)
         call prntm2(bc(ihessc(isa,isb)),nrowh,nhere,nrowh)
         if(isa.eq.isb.and.mynprocg.eq.1)then
          write(6,*)
     $         ('we have 1 proc ... average off diagonal elements ')
          iht=ibcoff
          ibcoff=iht+nrowh*nhere
          call enough('buildhesscas. 11',bc,ibc)
          do i=0,nrowh*nhere-1
           bc(iht+i)=bc(ihessc(isa,isb)+i)
          end do
          do i2=0,nvirt(isa)-1
           do i1=0,nvirt(isb)-1
            do i4=0,noc(isbo)-1
             do i3=0,noc(isao)-1
              ij=iht+i3+noc(isao)*(i4+noc(isbo)*(i2+nvirt(isa)*i1))
              ji=iht+i4+noc(isbo)*(i3+noc(isao)*(i1+nvirt(isb)*i2))
              avg=0.5d0*(bc(ij)+bc(ji))
              diff=bc(ij)-bc(ji)
              if(abs(diff).gt.1d-10)write(6,*)i3,i4,i1,i2,bc(ij),
     $             bc(ji),diff,avg,ij-iht+ihessc(isa,isb),
     $             ji-iht+ihessc(isa,isb)
              bc(ij)=avg
              bc(ji)=avg
             end do
            end do
           end do
          end do
          call prntm2(bc(iht),nrowh,nhere,nrowh)
          ibcoff=iht
         end if
         end if
         if(nsymb.eq.1)then
          iha=ibcoff
          nwds=nrowh*nvirtc(1)*nvirtc(1)                                4d11s17
          ibcoff=iha+nwds
          call enough('buildhesscas. 12',bc,ibc)
          do i1234=0,nwds-1
           bc(iha+i1234)=0d0
          end do
          do i12=0,nhere-1
           i12p=i12+ilh-1                                               4d11s17
           iad1=iha+nrowh*i12p                                          4d11s17
           iad2=ihessc(isa,isb)+nrowh*i12                                4d11s17
           do i34=0,nrowh-1
            bc(iad1+i34)=bc(iad2+i34)                                   4d11s17
           end do
          end do
          call dws_gsumf(bc(iha),nwds)
          if(isao.eq.isbo.and.isav.eq.isbv)then
          write(6,*)('let us average off diagonal elements ')
          do i2=0,nvirtc(isbv)-1
           do i1=0,nvirtc(isbv)-1
            ic12=i1+nvirtc(isbv)*i2
            ic21=i2+nvirtc(isbv)*i1
            do i4=0,noc(isbo)-1
             do i3=0,noc(isbo)-1
              ir34=i3+noc(isbo)*i4
              ir43=i4+noc(isbo)*i3
              iad=iha+ir34+nrowh*ic12
              iadt=iha+ir43+nrowh*ic21
              avg=0.5d0*(bc(iad)+bc(iadt))
              bc(iad)=avg
              bc(iadt)=avg
             end do
            end do
           end do
          end do
          end if
          call prntm2(bc(iha),nrowh,nvirtc(1)*nvirtc(1),nrowh)
          write(6,*)('full hess ')
          call printa(bc(iha),noc,0,noc,0,nvirtc,noc,nvirtc,noc,
     $         bc(ibcoff))
          ibcoff=iha
         end if
         end if
         ibcoff=ibtop                                                   3d29s16
        end if
       end do
      end do
      return
      end
