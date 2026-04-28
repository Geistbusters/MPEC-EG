c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine denmx12bodc(bc,ibc,nh0av,idoubo,irefo,                 3d26s25
     $     nbasdwsc,ncsf,iwave,mdon,mdoo,ixw1,ixw2,maxbx,maxbxd,srh,sr2,3d26s25
     $     multh,norb,nuniq,bodcp,oder1,noderw,ism,irel)                3d26s25
      implicit real*8 (a-h,o-z)
      external second                                                   10d17s23
c                                                                       3d3s21
c     compute two particle density matrices and compute contributions   3d17s25
c     to bodc.                                                          3d17s25
c
c     memory mapping ...
c     ioooo etc
c     vectors
c     compute densities (need to create space for 4v density) and
c     overwrite ioooo recompute all integrals, including 4v.
c     problem if number of roots is gt 1
c     compute orbital terms - store results for each der
c     clean up - just need density now.
c     compute densities in ao basis ...
c
      logical lwrite,ldebug                                             6d14s24
      integer*8 i18,i28,i48,ipack8                                      3d26s25
      character*10 dlabel
      include "common.hf"                                               7d28s22
      include "common.print"                                            6d14s24
      dimension ipack4(2)                                               3d26s25
      equivalence (ipack8,ipack4)                                       3d26s25
      dimension ncsf(*),irel(*),ism(*),irefo(*),multh(8,8),nh0av(*),
     $     iden(8),iorbno(8),id4o(idbk),jmden(idbk),iptoh(8,8,8),       3d19s25
     $     kmden(idbk),id1x(idbk),id3x(idbk),iorbf(8),nbasdwsc(*),      3d26s25
     $    iorb(8),iwave(*),idoubo(*),                                   3d26s25
     $     nl(4),ichoice(11),idatta(7,3),data(3),iamatu(8),xjt(512),    5d20s24
     $     iaddr(36,2,6),naddr(36,3),itt4v(8,8),isou(8),isoub(8),       8d23s24
     $     isoua1(8),i3xb(512),ionexb(512),jmatt(512),kmatt(512),       6d12s24
     $     ionexc(512),kmatd(64),i3x3(512),ionexbt(512),d4v(4),         7d11s24
     $     jmtsd(512),kmtsd(512),ioood(512),ionxd(512),i3xxd(512),      7d11s24
     $     i4xbd(512),idervcv(2,8),bodcp(*),oder1(noderw,*)             3d19s25
      include "common.store"
      include "common.basis"                                            4d29s24
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      common/drsigncm/drsign                                            8d20s24
      common/timerocm/tovrx,telapo(15)                                   4d26s18
      data loopx/200/
      integer*8 i,ibcoffo,id1x,id3x,id4o,iden,iorbf,iptrfbit,jmden,jptrf
     $kmden
      ldebug=iprtr(32).ne.0
      write(6,*)('Hi, my name is denmx12bodc'),ibcoff
      write(6,*)('mdoo '),mdoo,('mdon '),mdon
      mdoop=mdoo+1                                                      3d26s25
      icsf=1                                                            7d27s21
      icsf2=1                                                           7d27s21
      nff0k=iwave(4)+iwave(9)                                           7d27s21
      jff0k=iwave(4)+iwave(10)                                          7d27s21
      call countnfcn(ibc(nff0k),mdon,mdoo,nfcn)                         3d26s25
      nec=iwave(7)                                                      7d27s21
      ilook=iwave(4)+iwave(13)                                          7d27s21
      isymmrci=iwave(2)                                                 3d19s25
      nroot=iwave(3)                                                    3d19s25
      ipack8=ibc(ilook)                                                 5d12s21
      ncsftk=ipack4(1)
      ilook=ilook+1
      imff1=ilook+nroot*(ncsftk+1)                                      3d19s25
      write(6,*)('imff1 '),imff1,ibc(imff1)                                        3d19s25
      if(ldebug)write(6,*)('for ket, imff1 = '),imff1
      mff1=ibc(imff1)                                                   7d9s21
      inext=imff1+1                                                     7d21s21
      ihsdiag=1                                                         1d17s23
      nff1=1                                                            1d17s23
      iff1=1                                                            1d17s23
      if(mff1.gt.0)then                                                 7d9s21
       ihsdiag=imff1+1                                                  7d9s21
       nff1=ihsdiag+nsymb*mdoop*2                                       7d9s21
       iff1=nff1+nsymb*mdoop                                            7d9s21
       icsf=iff1+ibc(imff1)                                             7d9s21
       inext=icsf+mdoop-mdon                                            7d21s21
      end if
      mff2=ibc(inext)                                                   7d21s21
      ivdknon=1                                                         1d17s23
      if(ldebug)write(6,*)('for ket, mff2 = '),mff2
       ihddiag=1                                                        8d3s21
       nff2=1                                                           8d3s21
       iff2=1                                                           8d3s21
       icsf=1                                                           8d3s21
       icsf2=1                                                          8d3s21
       nfdat=1                                                          8d3s21
       ivdk=1                                                           8d3s21
      if(mff2.gt.0)then                                                 7d27s21
       ihddiag=inext+1                                                  7d21s21
       nff2=ihddiag+nsymb*mdoop*2                                       7d21s21
       iff2=nff2+nsymb*mdoop                                            7d21s21
       icsf=iff2+mff2                                                   7d21s21
       icsf2=icsf+mdoop-mdon                                            7d21s21
       nfdat=1                                                          8d3s21
       ivdk=1                                                           8d3s21
      else if(mff2.lt.0)then                                            8d3s21
       mdoubstore=inext+1
       ipack8=ibc(mdoubstore)                                           8d3s21
       ndoub=ipack4(1)                                                  8d12s21
       mdoub=ipack4(2)                                                  8d3s21
       mff2a=iabs(mff2)                                                 8d3s21
       iff2=mdoubstore+1                                                8d3s21
       nff2=iff2+mff2a                                                  8d3s21
       nfdat=nff2+mdoop*nsymb                                           8d3s21
       ivdk=nfdat+10*nsymb                                               8d3s21
       ivdknon=ivdk+nroot*ndoub                                         3d19s25
       icsf=ivdknon+nroot*mdoub                                         3d19s25
       icsf2=icsf+mdoop-mdon                                            8d3s21
      else                                                              8d3s21
      end if                                                            7d27s21
      loop=0
c     oo
      ichoice(1)=1
c     aa
      ichoice(2)=1
c     2e part of shift
      ichoice(3)=1
c     mixed oa
      ichoice(4)=1
c     d4o
      ichoice(5)=1
c     d1x
      ichoice(6)=1
c     dJ
      ichoice(7)=1
c     dK
      ichoice(8)=1
c     3x
      ichoice(9)=1
c     4x
      ichoice(10)=1
c     4v
      ichoice(11)=1
      if(ldebug)write(6,282)ichoice
  282 format('what we have for ichoice: ',11i2)
c
c     1: dcont
c     2: dorth
c     3: dorb
c     4: doortho
c     5: 1e ders
c     6: 2e ders
c
      ibcoffo=ibcoff                                                    3d3s21
      write(6,*)('stamp1 ')
      igoal=1962320
      write(6,*)('nsymb = '),nsymb
      write(6,*)('nsdlk '),nsdlk
      do isb=1,nsymb                                                    3d3s21
       iden(isb)=ibcoff                                                 3d3s21
       ibcoff=iden(isb)+nh0av(isb)*nh0av(isb)*nroot                     3d3s21
       iorbf(isb)=ibcoff                                                9d29s22
       ibcoff=iorbf(isb)+nh0av(isb)*nh0av(isb)*nroot                    9d29s22
      end do                                                            3d3s21
      do is=1,nsdlk                                                     7d28s22
       write(6,*)is,(isblk(j,is),j=1,4)
       if(isblk(1,is).eq.isblk(2,is))then                               7d28s22
        nrow=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2              7d28s22
        ncol=(irefo(isblk(3,is))*(irefo(isblk(3,is))+1))/2              7d28s22
       else                                                             7d28s22
        nrow=irefo(isblk(1,is))*irefo(isblk(2,is))                      7d28s22
        ncol=irefo(isblk(3,is))*irefo(isblk(4,is))                       7d28s22
       end if                                                           7d28s22
       id4o(is)=ibcoff                                                  7d28s22
       ibcoff=id4o(is)+nrow*ncol*nroot                                  7d28s22
      end do                                                            7d28s22
      ndstot=ibcoff-ibcoffo                                             5d2s23
      do is=1,nsdlk                                                     7d28s22
       if(isblk(1,is).eq.isblk(2,is))then                               7d28s22
        nrow=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2              7d28s22
       else                                                             7d28s22
        nrow=irefo(isblk(1,is))*irefo(isblk(2,is))                      7d28s22
       end if                                                           7d28s22
       call ilimts(nvirt(isblk(3,is)),nvirt(isblk(4,is)),mynprocg,      7d28s22
     $      mynowprog,il,ih,i1s,i1e,i2s,i2e)                            7d28s22
       ncol=ih+1-il                                                     7d28s22
       jmden(is)=ibcoff                                                 7d28s22
       ibcoff=jmden(is)+nrow*ncol*nroot                                 7d28s22
      end do                                                            7d28s22
      do is=1,nsdlk1                                                    7d28s22
       if(isblk1(1,is).eq.isblk1(2,is))then                             7d28s22
        nrow1=(irefo(isblk1(1,is))*(irefo(isblk1(1,is))+1))/2           7d28s22
       else                                                             7d28s22
        nrow1=irefo(isblk1(1,is))*irefo(isblk1(2,is))                   7d28s22
       end if                                                           7d28s22
       call ilimts(irefo(isblk1(3,is)),nvirt(isblk1(4,is)),mynprocg,    7d28s22
     $      mynowprog,il,ih,i1s,i1e,i2s,i2e)                            7d28s22
       ncol=ih+1-il                                                     7d28s22
       id1x(is)=ibcoff                                                  7d28s22
       ibcoff=id1x(is)+nrow1*irefo(isblk1(3,is))*nvirt(isblk1(4,is))    4d5s24
      end do                                                            7d28s22
      do is=1,nsdlk1                                                    7d28s22
       if(isblk1(1,is).eq.isblk1(2,is))then                             7d28s22
        nrow3=(nvirt(isblk1(1,is))*(nvirt(isblk1(1,is))+1))/2           7d28s22
       else                                                             7d28s22
        nrow3=nvirt(isblk1(1,is))*nvirt(isblk1(2,is))                   7d28s22
       end if                                                           7d28s22
       call ilimts(irefo(isblk1(3,is)),nvirt(isblk1(4,is)),mynprocg,    7d28s22
     $      mynowprog,il,ih,i1s,i1e,i2s,i2e)                            7d28s22
       ncol=ih+1-il                                                     7d28s22
       id3x(is)=ibcoff                                                  4d5s24
       ibcoff=id3x(is)+nrow3*ncol*nroot                                 7d28s22
      end do                                                            7d28s22
      do is=1,nsdlkk                                                    7d28s22
       nrow=irefo(isblkk(1,is))*irefo(isblkk(2,is))                     7d28s22
       call ilimts(nvirt(isblkk(3,is)),nvirt(isblkk(4,is)),mynprocg,    7d28s22
     $      mynowprog,il,ih,i1s,i1e,i2s,i2e)                            7d28s22
       ncol=ih+1-il                                                     7d28s22
       kmden(is)=ibcoff                                                 7d28s22
       ibcoff=kmden(is)+nrow*ncol*nroot                                 7d28s22
      end do                                                            7d28s22
      call enough('denmx12.  1',bc,ibc)
      do i=ibcoffo,ibcoff-1                                             3d3s21
       bc(i)=0d0                                                        3d3s21
      end do                                                            3d3s21
      write(6,*)('stamp2 ')
      if(lwrite)write(6,*)('computing densities')                       6d14s24
      jptrf=ibcoff                                                      3d26s25
      ibasis=ibcoff+jptrf+2*mdoop                                       3d26s25
      nww=nfcn*3                                                        3d26s25
      nww2=nww/2                                                        3d26s25
      if(nww2*2.ne.nww)nww2=nww2+1                                      3d26s25
      idorbf=ibasis+nww2                                                3d26s25
      nww=nfcn*norb                                                     3d26s25
      nww8=nww/8                                                        3d26s25
      if(nww8*8.ne.nww)nww8=nww8+1                                      3d26s25
      isorbf=idorbf+nww8                                                3d26s25
      ibcoff=isorbf+nww8                                                3d26s25
      write(6,*)('stamp3 ')
      call maptoold(mdon,mdoo,ibc(nff0k),ibc(jff0k),nfcn,ibc(ibasis),   3d26s25
     $     ibc(idorbf),ibc(isorbf),ibc(jptrf),norb,nec,bc,ibc)          3d26s25
      write(6,*)('stamp4 ')
      iptrfbit=ibcoff                                                   1d21s21
      ibcoff=iptrfbit+nsymb*mdoop                                       1d21s21
      call int1tobit(ibc(jptrf),ibc(iptrfbit),ibc(idorbf),ibc(isorbf),  1d21s21
     $     idorbfbit,isorbfbit,mdoop,1,nec,bc,ibc)                      11d10s22
      call second(time1)                                                10d17s23
      call hccsfd12(bc(ilook),ncsftk,ibc(ibasis),ncsf,ibc(iptrfbit),    3d26s25
     $     nfcn,iden,id4o,nroot,mdon,nec,isymmrci,ixw1,ixw2,mdoo,nh0av, 3d26s25
     $     iorbf,bc,ibc,1)                                              3d26s25
      write(6,*)('stamp5 ')
      call second(time2)                                                10d17s23
      telap=time2-time1-tovr                                            10d17s23
      if(lwrite)write(6,*)('time for hccsfd12 '),telap
      if(nsing.ne.0)then                                                3d4s21
       call second(time1)
       call hcsid12(ihsdiag,nff1,iff1,nff0,iff0,vint,nctf,ncsf,nec,     7d29s22
     $      mdon,mdoo,nsymb,multh,ixw1,ixw2,iden,nh0av,id1x,nvirt,      5d12s23
     $      .false.,maxbx,bc,ibc,igoal)                                       5d8s23
       call second(time2)
       telap=time2-time1-tovr
       if(lwrite)write(6,*)('time for hcsid12 '),telap
       call second(time1)
       call hcssd12(ihsdiag,nff1,iff1,ncsf,nec,mdon,mdoo,nsymb,multh,   5d12s23
     $      ixw1,ixw2,iden,nh0av,id4o,jmden,kmden,nvirt,nroot,ism,irel, 5d12s23
     $      irefo,isymmrci,norb,.false.,maxbx,dum,dum,bc,ibc,igoal)           5d12s23
       call second(time2)
       telap=time2-time1-tovr
       if(lwrite)write(6,*)('time for hcssd12'),telap                   6d14s24
       do is=1,nsdlk1
        if(isblk1(1,is).eq.isblk1(2,is))then                            10d16s23
         nrow=(irefo(isblk1(1,is))*(irefo(isblk1(1,is))+1))/2           10d16s23
        else                                                            10d16s23
         nrow=irefo(isblk1(1,is))*irefo(isblk1(2,is))                   10d16s23
        end if                                                          10d16s23
        ncol=irefo(isblk1(3,is))*nvirt(isblk1(4,is))                    10d16s23
        if(min(ncol,nrow).gt.0)then                                     10d16s23
         do iz=0,nrow*ncol-1                                            10d16s23
          bc(id1x(is)+iz)=bc(id1x(is)+iz)*2d0                           10d16s23
         end do                                                         10d16s23
        end if                                                          10d16s23
       end do
      end if                                                            3d4s21
      dot4v=0d0                                                         10d16s23
      do is=1,nsdlk                                                     7d10s23
       if(isblk(1,is).eq.isblk(2,is))then                               7d28s22
        nrow=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2              7d28s22
        ncol=(irefo(isblk(3,is))*(irefo(isblk(3,is))+1))/2              7d28s22
        write(6,*)('for id4o type '),(isblk(j,is),j=1,4)
        if(min(nrow,ncol).gt.0)then                                      7d10s23
         write(6,*)('density: ')
         call prntm2(bc(id4o(is)),nrow,ncol,nrow)
         call dws_gsumf(bc(id4o(is)),nrow*ncol)
         write(6,*)('global summed ')
         call prntm2(bc(id4o(is)),nrow,ncol,nrow)
         write(6,*)('oder1 '),loc(oder1)
         call prntm2(oder1,noderw,nuniq,noderw)
         do i=1,nuniq                                                    3d19s25
          trace=0d0                                                        7d10s23
          do i4=0,irefo(isblk(4,is))-1                                  3d26s25
           i4p=i4+idoubo(isblk(4,is))                                   3d26s25
           itri=((i4*(i4+1))/2)                                          3d19s25
           itri=id4o(is)+nrow*itri                                      3d19s25
           do i3=0,i4                                                   3d19s25
            i3p=i3+idoubo(isblk(4,is))                                  3d26s25
            iad34=i3p+1+nbasdwsc(isblk(3,is))*i4p                       3d26s25
            do j2=0,irefo(isblk(2,is))-1                                3d26s25
             jtri=((j2*(j2+1))/2)                                        3d19s25
             jtri=itri+jtri                                             3d19s25
             j2p=j2+idoubo(isblk(2,is))                                 3d26s25
             iad2=1+idoubo(isblk(2,is))+nbasdwsc(isblk(1,is))*j2p       3d26s25
             do j1=0,j2                                                 3d19s25
              orig=trace
              trace=trace+bc(jtri+j1)*oder1(iad2+j1,i)*oder1(iad34,i)   3d26s25
              if(abs(orig-trace).gt.1d-10)write(6,*)orig,
     $             bc(jtri+j1),oder1(iad2+j1,i),oder1(iad34,i),
     $             j1,j2,i3,i4,trace
             end do                                                     3d19s25
            end do                                                      3d19s25
           end do                                                       3d19s25
          end do                                                        3d19s25
          write(6,*)('for unique direction '),i                         3d19s25
          write(6,*)('2e- trace is '),trace,trace*0.5d0                 3d19s25
          bodcp(i)=bodcp(i)+trace*0.5d0                                 3d19s25
         end do                                                         3d19s25
        end if                                                          3d19s25
       end if                                                           7d28s22
      end do                                                            7d10s23
      ibcoff=ibcoffo                                                    3d17s25
      return
      end
      subroutine countnfcn(nff0k,mdon,mdoo,nfcn)                        3d26s25
      dimension nff0k(*)
      write(6,*)('Hi, my name is countnfcn '),mdon,mdoo
      nfcn=0                                                            3d26s25
      do nclo=mdon,mdoo                                                 3d26s25
       nclop=nclo+1                                                     3d26s25
       write(6,*)nclo,nff0k(nclop)                                      3d26s25
       nfcn=nfcn+nff0k(nclop)                                           3d26s25
      end do                                                            3d26s25
      return                                                            3d26s25
      end                                                               3d26s25
